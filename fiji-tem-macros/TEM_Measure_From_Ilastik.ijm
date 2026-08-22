// =============================================================================
//  TEM_Measure_From_Ilastik.ijm  = CODE 2 (automatic part)
//  Turns ilastik "Simple Segmentation" images into measured data - no tracing.
//  Run this AFTER ilastik has batch-segmented your images (see the steps I sent).
//
//  Per image it measures, in microns:
//    - ER, mitochondria, Golgi, vacuoles, lipid bodies: count, area, distance to
//      the nucleus, and near/far
//    - nucleus shape + count  -> multinucleation
//    - micronuclei (small separate nucleus-class bodies)
//
//  Outputs: Organelles.csv, Nuclei.csv, Micronuclei.csv, ImageSummary.csv
// =============================================================================

// --- CLASS ORDER: must match the order you added labels in ilastik ---
// (1=first label you added, 2=second, ...). Keep the nucleus class named "Nucleus".
var CLASS_VALUES = newArray( 1,          2,        3,    4,             5,      6,        7);
var CLASS_NAMES  = newArray("Background","Nucleus","ER","Mitochondria","Golgi","Vacuole","LipidBody");

// --- calibration: microns per pixel for each magnification (your microscope) ---
var CAL_MAG = newArray(2000,   2600,     11000,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.001355, 0.0006479);

// --- settings ---
var NEAR_FAR_UM      = 1.0;    // organelle centre <= this from nucleus = "near"
var MIN_OBJ_UM2      = 0.01;   // ignore organelle specks smaller than this
var MIN_NUC_UM2      = 2.0;    // Nucleus-class object >= this = a real nucleus
var MICRONUC_MIN_UM2 = 0.05;   // between this and MIN_NUC = a micronucleus
var PIXEL_SIZE_UM    = 0.005;  // last-resort fallback only

// --- state ---
var FIXED_PX = 0;
var nNucCount = 0;
var nMicroCount = 0;
var fallbackCount = 0;
var total = 0;
var idx = 0;
var _curImg = "";
var _curDir = "";
var DMAP_TITLE = "__nucDistMap";
var segRoot = "";

setBatchMode(true);
segRoot = getDirectory("Choose the folder of ilastik SEGMENTATION images");
outDir  = getDirectory("Choose an OUTPUT folder for the result CSVs");
loadCalibration(findCalib(segRoot, outDir));

// how is the pixel size decided?
Dialog.create("Scale for this run");
Dialog.addChoice("Pixel size comes from:",
    newArray("Magnification in each filename", "One magnification for the whole run"),
    "Magnification in each filename");
Dialog.addChoice("...if one magnification, which:", newArray("2000", "2600", "11000", "22000"), "2600");
Dialog.show();
scaleMode = Dialog.getChoice();
chosenMag = Dialog.getChoice();
if (scaleMode == "One magnification for the whole run") FIXED_PX = pxForMag(parseFloat(chosenMag));

freshTable("Organelles");
freshTable("Nuclei");
freshTable("Micronuclei");
freshTable("ImageSummary");
run("Set Measurements...", "area mean centroid perimeter shape redirect=None decimal=4");

total = countTiffs(segRoot);
idx = 0;
print("\\Clear");
print("Measuring " + total + " segmentation images...");
processTree(segRoot);

selectWindow("Organelles");   Table.save(outDir + "Organelles.csv");
selectWindow("Nuclei");       Table.save(outDir + "Nuclei.csv");
selectWindow("Micronuclei");  Table.save(outDir + "Micronuclei.csv");
selectWindow("ImageSummary"); Table.save(outDir + "ImageSummary.csv");
setBatchMode(false);
showMessage("Done", "Measured " + idx + "/" + total + " images.\n" +
    (fallbackCount > 0 ? ">> " + fallbackCount + " image(s) had no magnification match and used the fallback scale.\n" : "") +
    "CSVs saved in:\n" + outDir);

// =============================================================================
function processTree(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) processTree(dir + list[i]);
        else if (isTiff(list[i])) { idx = idx + 1; processSeg(dir + list[i], dir, list[i]); }
    }
}
function processSeg(path, dir, name) {
    _curImg = name; _curDir = dir;
    showProgress(idx, total);
    print("[" + idx + "/" + total + "] " + name);
    open(path);
    if (bitDepth() == 24) run("8-bit");
    calibrateSeg(name);
    seg = getTitle();

    nucVal = valueForName("Nucleus");
    nNucCount = 0;
    nMicroCount = 0;
    haveDmap = false;
    if (nucVal >= 0) {
        measureNucleiAndMicro(seg, nucVal, name);
        haveDmap = buildNucleusDistanceMap(seg, nucVal);
    }

    for (c = 0; c < CLASS_NAMES.length; c++) {
        nm = CLASS_NAMES[c];
        if (nm == "Background") continue;
        if (nm == "Nucleus") continue;
        measureOrganelleClass(seg, CLASS_VALUES[c], nm, name, haveDmap);
    }

    r = Table.size("ImageSummary");
    fld = relPath(dir);
    Table.set("Image",            r, name,        "ImageSummary");
    Table.set("Sample",           r, fld,         "ImageSummary");
    Table.set("NucleusCount",     r, nNucCount,   "ImageSummary");
    Table.set("Multinucleated",   r, (nNucCount > 1 ? 1 : 0), "ImageSummary");
    Table.set("MicronucleusCount",r, nMicroCount, "ImageSummary");
    Table.set("ER_count",         r, countOf("ER"),           "ImageSummary");
    Table.set("Mito_count",       r, countOf("Mitochondria"), "ImageSummary");
    Table.set("Golgi_count",      r, countOf("Golgi"),        "ImageSummary");
    Table.set("Vacuole_count",    r, countOf("Vacuole"),      "ImageSummary");
    Table.set("Lipid_count",      r, countOf("LipidBody"),    "ImageSummary");
    Table.update("ImageSummary");

    if (isOpen(DMAP_TITLE)) { selectWindow(DMAP_TITLE); close(); }
    selectWindow(seg); close();
    run("Clear Results");
}

function measureOrganelleClass(seg, val, oname, img, haveDmap) {
    selectWindow(seg);
    run("Duplicate...", "title=__cls");
    setThreshold(val, val);
    run("Convert to Mask");
    run("Set Measurements...", "area mean centroid perimeter shape redirect=None decimal=4");
    run("Clear Results");
    run("Analyze Particles...", "size=" + MIN_OBJ_UM2 + "-Infinity circularity=0.00-1.00 display clear");
    n = nResults;
    for (k = 0; k < n; k++) {
        cx = getResult("X", k);
        cy = getResult("Y", k);
        dist = NaN;
        nearFar = "NA";
        if (haveDmap) {
            dist = distanceFromEdge(cx, cy);
            if (dist <= NEAR_FAR_UM) nearFar = "near";
            else nearFar = "far";
        }
        r = Table.size("Organelles");
        fld = relPath(_curDir);
        Table.set("Image",   r, img,                    "Organelles");
        Table.set("Folder",  r, fld,                    "Organelles");
        Table.set("Type",    r, oname,                  "Organelles");
        Table.set("Index",   r, k + 1,                  "Organelles");
        Table.set("Area_um2",r, getResult("Area", k),   "Organelles");
        Table.set("Circularity", r, getResult("Circ.", k), "Organelles");
        Table.set("Dist_to_nucleus_um", r, dist,        "Organelles");
        Table.set("NearFar", r, nearFar,                "Organelles");
        Table.update("Organelles");
    }
    selectWindow("__cls"); close();
    selectWindow(seg);
}

// nucleus-class objects: big = nucleus, small = micronucleus
function measureNucleiAndMicro(seg, val, img) {
    selectWindow(seg);
    run("Duplicate...", "title=__nuc");
    setThreshold(val, val);
    run("Convert to Mask");
    run("Fill Holes");
    run("Clear Results");
    run("Analyze Particles...", "size=" + MICRONUC_MIN_UM2 + "-Infinity circularity=0.00-1.00 display clear");
    n = nResults;
    nNucCount = 0;
    nMicroCount = 0;
    for (k = 0; k < n; k++) {
        a = getResult("Area", k);
        fld = relPath(_curDir);
        if (a >= MIN_NUC_UM2) {
            nNucCount = nNucCount + 1;
            r = Table.size("Nuclei");
            Table.set("Image",       r, img,                   "Nuclei");
            Table.set("Folder",      r, fld,                   "Nuclei");
            Table.set("NucleusIndex",r, nNucCount,             "Nuclei");
            Table.set("Area_um2",    r, a,                     "Nuclei");
            Table.set("Perimeter_um",r, getResult("Perim.", k),"Nuclei");
            Table.set("Circularity", r, getResult("Circ.", k), "Nuclei");
            Table.set("AspectRatio", r, getResult("AR", k),    "Nuclei");
            Table.set("Roundness",   r, getResult("Round", k), "Nuclei");
            Table.set("Solidity",    r, getResult("Solidity", k),"Nuclei");
            Table.update("Nuclei");
        } else {
            nMicroCount = nMicroCount + 1;
            r = Table.size("Micronuclei");
            Table.set("Image",       r, img,                   "Micronuclei");
            Table.set("Folder",      r, fld,                   "Micronuclei");
            Table.set("Index",       r, nMicroCount,           "Micronuclei");
            Table.set("Area_um2",    r, a,                     "Micronuclei");
            Table.set("Circularity", r, getResult("Circ.", k), "Micronuclei");
            Table.update("Micronuclei");
        }
    }
    selectWindow("__nuc"); close();
    selectWindow(seg);
}

function buildNucleusDistanceMap(seg, val) {
    selectWindow(seg);
    run("Duplicate...", "title=" + DMAP_TITLE);
    setThreshold(val, val);
    run("Convert to Mask");
    run("Fill Holes");
    getStatistics(area, mean);
    if (mean == 0) { selectWindow(DMAP_TITLE); close(); selectWindow(seg); return false; }
    setOption("BlackBackground", true);
    run("Invert");
    run("Distance Map");
    selectWindow(seg);
    return true;
}
function distanceFromEdge(cx, cy) {
    if (!isOpen(DMAP_TITLE)) return NaN;
    cur = getTitle();
    selectWindow(DMAP_TITLE);
    getPixelSize(u, pw, ph);
    px = round(cx / pw);
    py = round(cy / ph);
    if (px < 0) px = 0;
    if (py < 0) py = 0;
    if (px >= getWidth())  px = getWidth() - 1;
    if (py >= getHeight()) py = getHeight() - 1;
    dpx = getPixel(px, py);
    selectWindow(cur);
    return dpx * pw;
}

// set the pixel size on the open segmentation image
function calibrateSeg(name) {
    pw = PIXEL_SIZE_UM;
    src = "fallback";
    if (FIXED_PX > 0) { pw = FIXED_PX; src = "fixed"; }
    else {
        mag = parseMag(name);
        if (mag > 0) {
            p = pxForMag(mag);
            if (p > 0) { pw = p; src = "mag"; }
        }
    }
    if (src == "fallback") fallbackCount = fallbackCount + 1;
    setVoxelSize(pw, pw, 1, "micron");
    return pw;
}

// ---- calibration table helpers ----
function findCalib(root, outDir) {
    cands = newArray(outDir + "magnification_calibration.csv", root + "magnification_calibration.csv");
    for (i = 0; i < cands.length; i++)
        if (File.exists(cands[i])) return cands[i];
    return "";
}
function loadCalibration(path) {
    if (path == "" || !File.exists(path)) return CAL_MAG.length > 0;   // keep built-in
    CAL_MAG = newArray(0);
    CAL_PX = newArray(0);
    lines = split(File.openAsString(path), "\n");
    for (i = 0; i < lines.length; i++) {
        ln = String.trim(lines[i]);
        if (ln == "") continue;
        c = split(ln, ",");
        if (c.length < 2) continue;
        a = String.trim(c[0]);
        if (!matches(a, "[0-9].*")) continue;
        CAL_MAG = Array.concat(CAL_MAG, parseFloat(a));
        CAL_PX = Array.concat(CAL_PX, parseFloat(String.trim(c[1])));
    }
    return CAL_MAG.length > 0;
}
function pxForMag(mag) {
    for (i = 0; i < CAL_MAG.length; i++)
        if (CAL_MAG[i] == mag) return CAL_PX[i];
    return -1;
}
function parseMag(name) {
    s = name;
    if (matches(s, ".*[0-9]+[kK]?[xX].*")) {
        num = replace(s, ".*?([0-9]+)([kK]?)[xX].*", "$1");
        kf = replace(s, ".*?([0-9]+)([kK]?)[xX].*", "$2");
        v = parseFloat(num);
        if (kf == "k" || kf == "K") v = v * 1000;
        return v;
    }
    if (matches(s, ".*[xX][0-9]+.*"))
        return parseFloat(replace(s, ".*[xX]([0-9]+).*", "$1"));
    return -1;
}

// ---- small utilities ----
function isTiff(name) {
    n = toLowerCase(name);
    return endsWith(n, ".tif") || endsWith(n, ".tiff");
}
function countTiffs(dir) {
    c = 0;
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) c = c + countTiffs(dir + list[i]);
        else if (isTiff(list[i])) c = c + 1;
    }
    return c;
}
function valueForName(nm) {
    for (i = 0; i < CLASS_NAMES.length; i++)
        if (CLASS_NAMES[i] == nm) return CLASS_VALUES[i];
    return -1;
}
function relPath(dir) {
    parts = split(dir, "/");
    n = parts.length;
    if (n >= 2) return parts[n-2] + "/" + parts[n-1];
    if (n == 1) return parts[0];
    return dir;
}
function countOf(otype) {
    n = Table.size("Organelles");
    c = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, "Organelles") == _curImg && Table.getString("Type", i, "Organelles") == otype)
            c = c + 1;
    return c;
}
function freshTable(name) {
    if (isOpen(name)) { selectWindow(name); run("Close"); }
    Table.create(name);
}
