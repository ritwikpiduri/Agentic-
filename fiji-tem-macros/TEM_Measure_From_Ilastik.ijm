// =============================================================================
//  TEM_Measure_From_Ilastik.ijm
//  Turns ilastik "Simple Segmentation" label images into fully-measured,
//  no-tracing organelle data. Run this AFTER ilastik has batch-segmented your
//  dataset (see ilastik/README_ilastik_pipeline.md and run_ilastik_batch.sh).
//
//  INPUT : a folder tree of ilastik segmentation TIFFs, where every pixel value
//          is a CLASS index (1..N) matching the classes you trained. One seg
//          image per original TEM image, same relative path.
//  OUTPUT: Organelles.csv + ImageSummary.csv (same schema as the interactive
//          macro), so tem_stats.py runs on them unchanged. Also Nuclei.csv with
//          the auto nucleus shape.
//
//  Per image it:
//    1. builds a distance map from the NUCLEUS class (for near/far),
//    2. for each organelle class, finds every object (connected component),
//       records count, area, circularity, centroid, distance to nucleus edge,
//       and near/far,
//    3. measures nucleus shape + multinucleation.
//
//  NO MANUAL TRACING. Spot-check a handful against the raw images; if a class
//  is over/under-segmented, add more training labels in ilastik and re-run.
// =============================================================================

// -------------------- MUST MATCH YOUR ILASTIK CLASS ORDER --------------------
// ilastik Simple Segmentation labels classes 1,2,3,... in the order you added
// them. List the SAME order here. Put background first. Give the nucleus class
// the exact name "Nucleus" so distance measurement finds it.
var CLASS_VALUES = newArray( 1,          2,        3,    4,             5,      6,        7);
var CLASS_NAMES  = newArray("Background","Nucleus","ER","Mitochondria","Golgi","Vacuole","LipidBody");
// -----------------------------------------------------------------------------

// -------------------- USER SETTINGS --------------------
var PIXEL_SIZE_UM = 0.005;    // <-- EDIT: microns/pixel (same as segmentation)
var UNIT          = "micron";
var NEAR_FAR_UM   = 1.0;      // organelle centroid <= this from nucleus = "near"
var MIN_OBJ_UM2   = 0.01;     // ignore specks smaller than this (noise)
var MIN_NUC_UM2   = 2.0;      // min area to count something as a nucleus
// -------------------------------------------------------

var DMAP_TITLE = "__nucDistMap";
var ORIG_ROOT  = "";         // set below; used to read the real scale per image
var fallbackCount = 0;

setBatchMode(true);
segRoot = getDirectory("Choose the folder of ilastik SEGMENTATION images");
// ilastik label images lose the original scale, so read it from the ORIGINAL
// TIFFs (which carry it in their metadata). Point to the raw dataset root; its
// folder structure must mirror the segmentation folder.
showMessage("Original images",
    "Next, choose the ROOT of your ORIGINAL TIFFs so the real scale can be read\n" +
    "from each image's metadata. (Cancel to use the PIXEL_SIZE_UM fallback.)");
ORIG_ROOT = getDirectory("Choose the ORIGINAL TIFF root (or Cancel for fallback)");
if (ORIG_ROOT == "0" || ORIG_ROOT == "") ORIG_ROOT = "";
outDir  = getDirectory("Choose an OUTPUT folder for CSVs");
loadCalibration(findCalib(segRoot, outDir));  // magnification -> pixel-size table

freshTable("Organelles");
freshTable("Nuclei");
freshTable("ImageSummary");
run("Set Measurements...", "area mean centroid perimeter shape redirect=None decimal=4");

var total = 0;
var idx = 0;
total = countTiffs(segRoot);
print("\\Clear");
print("Measuring " + total + " segmentation images from ilastik...");
processTree(segRoot);

selectWindow("Organelles");   Table.save(outDir + "Organelles.csv");
selectWindow("Nuclei");       Table.save(outDir + "Nuclei.csv");
selectWindow("ImageSummary"); Table.save(outDir + "ImageSummary.csv");
setBatchMode(false);
showMessage("Done", "Measured " + idx + "/" + total + " images.\n" +
    (fallbackCount > 0 ?
        ">> " + fallbackCount + " image(s) used the PIXEL_SIZE_UM fallback\n" +
        "   (no matching original / uncalibrated). Areas & distances for\n" +
        "   those may be off.\n" : "") +
    "CSVs in:\n" + outDir);

// =============================================================================
function processTree(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) processTree(dir + list[i]);
        else if (isTiff(list[i])) { idx++; processSeg(dir + list[i], dir, list[i]); }
    }
}
var _curImg = "";
var _curDir = "";
function processSeg(path, dir, name) {
    _curImg = name; _curDir = dir;
    showProgress(idx, total);
    print("[" + idx + "/" + total + "] " + relPath(dir) + " -> " + name);
    open(path);
    if (bitDepth() == 24) run("8-bit");
    calibrateSegFromOriginal(dir, name);   // sets real scale on the seg image
    seg = getTitle();

    // ---- nucleus: distance map + shape + count ----
    nucVal = valueForName("Nucleus");
    nNuclei = 0;
    haveDmap = false;
    if (nucVal >= 0) {
        nNuclei = measureNuclei(seg, nucVal, name);
        haveDmap = buildNucleusDistanceMap(seg, nucVal);
    }

    // ---- organelles ----
    counts = newArray(CLASS_NAMES.length);
    for (c = 0; c < CLASS_NAMES.length; c++) {
        nm = CLASS_NAMES[c];
        if (nm == "Background" || nm == "Nucleus") continue;
        counts[c] = measureOrganelleClass(seg, CLASS_VALUES[c], nm, name, haveDmap);
    }

    // ---- per-image summary row ----
    r = Table.size("ImageSummary");
    Table.set("Image",            r, name,              "ImageSummary");
    Table.set("Sample",           r, relPath(dir),      "ImageSummary");
    Table.set("NucleusCount",     r, nNuclei,           "ImageSummary");
    Table.set("Multinucleated",   r, (nNuclei > 1 ? 1 : 0), "ImageSummary");
    Table.set("MicronucleusCount",r, 0,                 "ImageSummary");  // see note
    Table.set("ER_count",         r, countOf("ER"),          "ImageSummary");
    Table.set("Mito_count",       r, countOf("Mitochondria"),"ImageSummary");
    Table.set("Golgi_count",      r, countOf("Golgi"),       "ImageSummary");
    Table.set("Vacuole_count",    r, countOf("Vacuole"),     "ImageSummary");
    Table.set("Lipid_count",      r, countOf("LipidBody"),   "ImageSummary");
    Table.update("ImageSummary");

    if (isOpen(DMAP_TITLE)) { selectWindow(DMAP_TITLE); close(); }
    selectWindow(seg); close();
    run("Clear Results");
}

// count organelle objects of one class in the CURRENT segmentation, measuring
// each and appending to Organelles table. Returns the count.
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
        cx = getResult("X", k); cy = getResult("Y", k);
        dist = NaN; nearFar = "NA";
        if (haveDmap) { dist = distanceFromEdge(cx, cy); nearFar = (dist <= NEAR_FAR_UM) ? "near" : "far"; }
        r = Table.size("Organelles");
        Table.set("Image",   r, img,                    "Organelles");
        Table.set("Folder",  r, relPath(_curDir),       "Organelles");
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
    return n;
}

function measureNuclei(seg, val, img) {
    selectWindow(seg);
    run("Duplicate...", "title=__nuc");
    setThreshold(val, val);
    run("Convert to Mask");
    run("Fill Holes");
    run("Clear Results");
    run("Analyze Particles...", "size=" + MIN_NUC_UM2 + "-Infinity circularity=0.00-1.00 display clear");
    n = nResults;
    for (k = 0; k < n; k++) {
        r = Table.size("Nuclei");
        Table.set("Image",       r, img,                   "Nuclei");
        Table.set("Folder",      r, relPath(_curDir),      "Nuclei");
        Table.set("NucleusIndex",r, k + 1,                 "Nuclei");
        Table.set("Area_um2",    r, getResult("Area", k),  "Nuclei");
        Table.set("Perimeter_um",r, getResult("Perim.", k),"Nuclei");
        Table.set("Circularity", r, getResult("Circ.", k), "Nuclei");
        Table.set("AspectRatio", r, getResult("AR", k),    "Nuclei");
        Table.set("Roundness",   r, getResult("Round", k), "Nuclei");
        Table.set("Solidity",    r, getResult("Solidity", k),"Nuclei");
        Table.update("Nuclei");
    }
    selectWindow("__nuc"); close();
    selectWindow(seg);
    return n;
}

// distance map: outside-nucleus pixels get calibrated distance to nucleus edge
function buildNucleusDistanceMap(seg, val) {
    selectWindow(seg);
    run("Duplicate...", "title=" + DMAP_TITLE);
    setThreshold(val, val);
    run("Convert to Mask");           // nucleus = 255
    run("Fill Holes");
    getStatistics(area, mean);
    if (mean == 0) { selectWindow(DMAP_TITLE); close(); selectWindow(seg); return false; }
    setOption("BlackBackground", true);
    run("Invert");                    // outside = 255
    run("Distance Map");              // outside pixels -> pixels-to-nucleus
    selectWindow(seg);
    return true;
}
function distanceFromEdge(cx, cy) {
    if (!isOpen(DMAP_TITLE)) return NaN;
    cur = getTitle();
    selectWindow(DMAP_TITLE);
    getPixelSize(u, pw, ph);                  // real scale carried on the map
    px = round(cx / pw); py = round(cy / ph);
    if (px < 0) px = 0; if (py < 0) py = 0;
    if (px >= getWidth())  px = getWidth() - 1;
    if (py >= getHeight()) py = getHeight() - 1;
    dpx = getPixel(px, py);
    selectWindow(cur);
    return dpx * pw;
}
// Read the true pixel size from the matching ORIGINAL TIFF and stamp it on the
// currently-open segmentation image. Falls back to PIXEL_SIZE_UM if no original
// is found or it is uncalibrated.
function calibrateSegFromOriginal(dir, name) {
    pw = PIXEL_SIZE_UM; src = "fallback";
    // (1) try the original image's embedded metadata scale
    if (ORIG_ROOT != "") {
        op = originalPath(dir, name);
        if (op != "") {
            seg = getTitle();
            open(op);
            if (normalizeToMicron()) { getPixelSize(u, p2, ph2); pw = p2; src = "embedded"; }
            close();                       // close the original
            selectWindow(seg);
        }
    }
    // (2) else read magnification from the filename -> calibration table
    if (src == "fallback") {
        mag = parseMag(name);
        if (mag > 0) { px = pxForMag(mag); if (px > 0) { pw = px; src = "mag=" + mag; } }
    }
    if (src == "fallback") fallbackCount++;
    setVoxelSize(pw, pw, 1, "micron");
    return pw;
}
// ---- magnification-from-filename calibration table ----
var CAL_MAG = newArray(0);
var CAL_PX  = newArray(0);
function findCalib(root, outDir) {
    cands = newArray(outDir + "magnification_calibration.csv",
                     root   + "magnification_calibration.csv");
    for (i = 0; i < cands.length; i++) if (File.exists(cands[i])) return cands[i];
    Dialog.create("Calibration table");
    Dialog.addString("Path to magnification_calibration.csv (blank = fallback):", "", 60);
    Dialog.show();
    return String.trim(Dialog.getString());
}
function loadCalibration(path) {
    CAL_MAG = newArray(0); CAL_PX = newArray(0);
    if (path == "" || !File.exists(path)) return false;
    lines = split(File.openAsString(path), "\n");
    for (i = 0; i < lines.length; i++) {
        ln = String.trim(lines[i]); if (ln == "") continue;
        c = split(ln, ","); if (c.length < 2) continue;
        a = String.trim(c[0]); if (!matches(a, "[0-9].*")) continue;
        CAL_MAG = Array.concat(CAL_MAG, parseFloat(a));
        CAL_PX  = Array.concat(CAL_PX, parseFloat(String.trim(c[1])));
    }
    return CAL_MAG.length > 0;
}
function pxForMag(mag) {
    for (i = 0; i < CAL_MAG.length; i++) if (CAL_MAG[i] == mag) return CAL_PX[i];
    return -1;
}
function parseMag(name) {
    s = name;
    if (matches(s, ".*[0-9]+[.]?[0-9]*[kK]?[xX].*")) {
        num = replace(s, ".*?([0-9]+[.]?[0-9]*)([kK]?)[xX].*", "$1");
        kfl = replace(s, ".*?([0-9]+[.]?[0-9]*)([kK]?)[xX].*", "$2");
        v = parseFloat(num); if (kfl == "k" || kfl == "K") v = v * 1000; return v;
    }
    if (matches(s, ".*[xX][0-9]+.*"))
        return parseFloat(replace(s, ".*[xX]([0-9]+).*", "$1"));
    return -1;
}
// map a segmentation path (…/base_seg.tif) back to its original TIFF under ORIG_ROOT
function originalPath(dir, name) {
    relDir = substring(dir, lengthOf(segRoot));    // path under segRoot
    b = stripExt(name);
    if (endsWith(b, "_seg")) b = substring(b, 0, lengthOf(b) - 4);
    cand = newArray(ORIG_ROOT + relDir + b + ".tif",
                    ORIG_ROOT + relDir + b + ".tiff",
                    ORIG_ROOT + relDir + b + ".TIF",
                    ORIG_ROOT + relDir + b + ".TIFF");
    for (i = 0; i < cand.length; i++) if (File.exists(cand[i])) return cand[i];
    return "";
}
function stripExt(nm) { d = lastIndexOf(nm, "."); return (d > 0) ? substring(nm, 0, d) : nm; }
// convert the open image's calibration to microns; false if uncalibrated
function normalizeToMicron() {
    getPixelSize(unit, pw, ph);
    u = toLowerCase(unit);
    if (u == "micron" || u == "microns" || u == "um" || u == "µm")
        return (pw != 1);
    factor = 0;
    if (u == "nm" || u == "nanometer" || u == "nanometre" || u == "nanometers")
        factor = 1.0/1000.0;
    else if (u == "a" || u == "angstrom" || u == "ang" || u == "å")
        factor = 1.0/10000.0;
    else if (u == "mm" || u == "millimeter")
        factor = 1000.0;
    else return false;
    setVoxelSize(pw * factor, ph * factor, 1, "micron");
    return true;
}

// ---- small utilities ----
function isTiff(name) { n = toLowerCase(name); return endsWith(n, ".tif") || endsWith(n, ".tiff"); }
function countTiffs(dir) {
    c = 0; list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) c += countTiffs(dir + list[i]);
        else if (isTiff(list[i])) c++;
    }
    return c;
}
function valueForName(nm) {
    for (i = 0; i < CLASS_NAMES.length; i++) if (CLASS_NAMES[i] == nm) return CLASS_VALUES[i];
    return -1;
}
function relPath(dir) {
    parts = split(dir, "/"); n = parts.length;
    if (n >= 2) return parts[n-2] + "/" + parts[n-1];
    if (n == 1) return parts[0];
    return dir;
}
// count Organelles rows of one type for the image currently being processed
function countOf(otype) {
    n = Table.size("Organelles"); c = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, "Organelles") == _curImg &&
            Table.getString("Type",  i, "Organelles") == otype) c++;
    return c;
}
function freshTable(name) {
    if (isOpen(name)) { selectWindow(name); run("Close"); }
    Table.create(name);
}
