// =============================================================================
//  TEM_Manual_Code2.ijm  = CODE 2 (manual marking, everything else automatic)
//  It OPENS each image in a folder for you, auto-calibrates from the filename
//  magnification, and gives a per-image menu. You only identify/trace/click the
//  structures; it measures, counts, and logs everything.
//
//  Run it (no image needs to be open). Pick your image folder + an output folder.
//  For each image: trace nucleus, mark organelles/micronuclei, measure membrane,
//  click pores, then "Next image". Skip images you don't want. Quit saves all.
//
//  Output CSVs: Nuclei, Organelles, Micronuclei, Membrane, Pores, ImageSummary
//  (all in um/nm; feed into the stats macro to compare your 4 samples).
// =============================================================================

var NEAR_FAR_UM = 1.0;          // organelle centre <= this (um) from nucleus = "near"
var DMAP = "__nucDistMap";
var HAVE_DMAP = false;
var _curImg = "";
var _curDir = "";

// --- calibration ---
// Exact pixel sizes (um/pixel) you gave for these magnifications:
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
// For ANY other magnification (e.g. 4500), pixel size is estimated from your
// camera constant: magnification x pixel-size(nm) ~ 14400  ->  um/px = 14.4 / mag
var CAM_CONST = 14.4;

// ---- pick folders ----
inDir  = getDirectory("Choose the FOLDER of images to analyse");
outDir = getDirectory("Choose an OUTPUT folder for the CSVs");
initTables();

// build a flat list of all tiffs (including subfolders)
var files = newArray(0);
files = listTiffs(inDir);
if (files.length == 0) { showMessage("No .tif/.tiff images found in that folder."); exit; }

setBatchMode(false);
quit = false;
done = 0;
for (f = 0; f < files.length; f++) {
    if (quit) f = files.length;
    else {
        path = files[f];
        open(path);
        if (bitDepth() == 24 || bitDepth() == 16) run("8-bit");   // grayscale for measuring
        _curImg = getTitle();
        _curDir = File.getParent(path);
        HAVE_DMAP = false;
        autoCalibrate();                       // sets scale from filename magnification

        // per-image action menu
        imgDone = false;
        while (!imgDone) {
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg);
            Dialog.addChoice("Action:", newArray(
                "Trace nucleus (shape + chromatin + multinucleation)",
                "Mark organelles (type, size, near/far)",
                "Mark micronuclei",
                "Membrane / perinuclear thickness",
                "Nuclear pores (count + density)",
                "NEXT image (save this image's summary)",
                "SKIP this image",
                "QUIT and save all CSVs"), "Trace nucleus (shape + chromatin + multinucleation)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Trace nucleus"))   actNucleus();
            else if (startsWith(a, "Mark organelles")) actOrganelles();
            else if (startsWith(a, "Mark micronuclei"))actMicronuclei();
            else if (startsWith(a, "Membrane"))        actMembrane();
            else if (startsWith(a, "Nuclear pores"))   actPores();
            else if (startsWith(a, "NEXT"))            { actSummary(); imgDone = true; }
            else if (startsWith(a, "SKIP"))            { imgDone = true; }
            else if (startsWith(a, "QUIT"))            { imgDone = true; quit = true; }
        }
        if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
saveAll(outDir);
showMessage("Done", "Processed " + done + " image(s).\nCSVs saved in:\n" + outDir);

// =============================================================================
function autoCalibrate() {
    mag = parseMag(getTitle());
    if (mag <= 0) return false;
    px = pxForMag(mag);
    setVoxelSize(px, px, 1, "micron");
    return true;
}
function parseMag(name) {
    s = name;
    if (matches(s, ".*[0-9]+[ ]?[kK]?[ ]?[xX].*")) {
        num = replace(s, ".*?([0-9]+)[ ]?([kK]?)[ ]?[xX].*", "$1");
        kf  = replace(s, ".*?([0-9]+)[ ]?([kK]?)[ ]?[xX].*", "$2");
        v = parseFloat(num);
        if (kf == "k" || kf == "K") v = v * 1000;
        return v;
    }
    return -1;
}
function pxForMag(mag) {
    for (i = 0; i < CAL_MAG.length; i++) if (CAL_MAG[i] == mag) return CAL_PX[i];
    if (mag > 0) return CAM_CONST / mag;   // estimate for magnifications not listed
    return -1;
}

function actNucleus() {
    setMeas(); setTool("freehand");
    waitForUser("Trace the NUCLEUS",
        "Outline the nucleus (freehand), then OK.\n(Repeat for each extra nucleus in this image.)");
    if (selectionType() < 0) { showMessage("No selection."); return; }
    run("Measure"); row = nResults - 1;
    n = countImg("Nuclei", _curImg) + 1; r = Table.size("Nuclei");
    Table.set("Image", r, _curImg, "Nuclei");
    Table.set("Folder", r, folderOf(), "Nuclei");
    Table.set("NucleusIndex", r, n, "Nuclei");
    Table.set("Area_um2", r, getResult("Area", row), "Nuclei");
    Table.set("Perimeter_um", r, getResult("Perim.", row), "Nuclei");
    Table.set("Circularity", r, getResult("Circ.", row), "Nuclei");
    Table.set("AspectRatio", r, getResult("AR", row), "Nuclei");
    Table.set("Roundness", r, getResult("Round", row), "Nuclei");
    Table.set("Solidity", r, getResult("Solidity", row), "Nuclei");
    getStatistics(nArea, nMean, nMin, nMax, nStd);
    thr = round(nMean - 0.5 * nStd);
    Roi.getBounds(bx, by, bw, bh);
    dark = 0; tot = 0;
    for (yy = by; yy < by + bh; yy++) {
        for (xx = bx; xx < bx + bw; xx++) {
            if (Roi.contains(xx, yy)) { tot = tot + 1; if (getPixel(xx, yy) <= thr) dark = dark + 1; }
        }
    }
    het = NaN; if (tot > 0) het = 100.0 * dark / tot;
    Table.set("Heterochrom_pct", r, het, "Nuclei");
    Table.set("Euchrom_pct", r, 100 - het, "Nuclei");
    Table.update("Nuclei");
    buildDistMap(); run("Select None");
    showMessage("Nucleus " + n + " recorded." + (n > 1 ? "\n>> Multinucleated." : ""));
}

function actOrganelles() {
    types = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody");
    Dialog.create("Organelle type");
    Dialog.addChoice("Marking which organelle?", types, types[0]);
    Dialog.addNumber("Near/far cutoff (um):", NEAR_FAR_UM);
    Dialog.show();
    otype = Dialog.getChoice(); NEAR_FAR_UM = Dialog.getNumber();
    if (!HAVE_DMAP) showMessage("Tip: trace the nucleus first for distances.");
    setMeas(); setTool("freehand");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser(otype + " #" + i, "Trace this " + otype + " (freehand), then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); row = nResults - 1;
            cx = getResult("X", row); cy = getResult("Y", row);
            dist = NaN; nf = "NA";
            if (HAVE_DMAP) { dist = distEdge(cx, cy); if (dist <= NEAR_FAR_UM) nf = "near"; else nf = "far"; }
            r = Table.size("Organelles");
            Table.set("Image", r, _curImg, "Organelles");
            Table.set("Folder", r, folderOf(), "Organelles");
            Table.set("Type", r, otype, "Organelles");
            Table.set("Index", r, i, "Organelles");
            Table.set("Area_um2", r, getResult("Area", row), "Organelles");
            Table.set("Dist_to_nucleus_um", r, dist, "Organelles");
            Table.set("NearFar", r, nf, "Organelles");
            Table.update("Organelles");
        }
        more = getBoolean("Mark another " + otype + "?");
    }
    run("Select None"); showMessage(i + " " + otype + "(s) recorded.");
}

function actMicronuclei() {
    setMeas(); setTool("freehand");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser("Micronucleus #" + i, "Trace a small body separate from the nucleus, then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); row = nResults - 1; r = Table.size("Micronuclei");
            Table.set("Image", r, _curImg, "Micronuclei");
            Table.set("Folder", r, folderOf(), "Micronuclei");
            Table.set("Index", r, i, "Micronuclei");
            Table.set("Area_um2", r, getResult("Area", row), "Micronuclei");
            Table.set("Circularity", r, getResult("Circ.", row), "Micronuclei");
            Table.update("Micronuclei");
        }
        more = getBoolean("Mark another micronucleus?");
    }
    run("Select None"); showMessage(i + " micronucleus/nuclei recorded.");
}

function actMembrane() {
    setTool("line");
    n = getNumber("How many thickness measurements on this image?", 5);
    i = 0;
    while (i < n) {
        waitForUser("Thickness " + (i+1) + " of " + n, "Draw a short line ACROSS the membrane, then OK.");
        if (selectionType() == 5) {
            lenUm = getValue("Length"); r = Table.size("Membrane");
            Table.set("Image", r, _curImg, "Membrane");
            Table.set("Folder", r, folderOf(), "Membrane");
            Table.set("Measurement", r, i+1, "Membrane");
            Table.set("Thickness_nm", r, lenUm*1000.0, "Membrane");
            Table.set("Thickness_um", r, lenUm, "Membrane");
            Table.update("Membrane");
            i = i + 1;
        }
    }
    run("Select None"); showMessage(n + " thickness measurements recorded.");
}

function actPores() {
    setTool("multipoint");
    waitForUser("Nuclear pores", "Multi-point tool: CLICK each nuclear pore, then OK.");
    if (selectionType() != 10) { showMessage("Use the multi-point tool and click the pores."); return; }
    getSelectionCoordinates(xs, ys); count = xs.length;
    perim = nucPerim(_curImg); dens = NaN; if (perim > 0) dens = count / perim;
    r = Table.size("Pores");
    Table.set("Image", r, _curImg, "Pores");
    Table.set("Folder", r, folderOf(), "Pores");
    Table.set("PoreCount", r, count, "Pores");
    Table.set("NucPerimeter_um", r, perim, "Pores");
    Table.set("Pores_per_um", r, dens, "Pores");
    Table.update("Pores");
    run("Select None");
    showMessage(count + " pores." + (perim > 0 ? "\nDensity = " + d2s(dens,3) + " pores/um." : "\n(Trace nucleus first for density.)"));
}

function actSummary() {
    nNuc = countImg("Nuclei", _curImg); r = Table.size("ImageSummary");
    Table.set("Image", r, _curImg, "ImageSummary");
    Table.set("Sample", r, folderOf(), "ImageSummary");
    Table.set("NucleusCount", r, nNuc, "ImageSummary");
    Table.set("Multinucleated", r, (nNuc > 1 ? 1 : 0), "ImageSummary");
    Table.set("MicronucleusCount", r, countImg("Micronuclei", _curImg), "ImageSummary");
    Table.set("ER_count", r, countType(_curImg, "ER"), "ImageSummary");
    Table.set("Mito_count", r, countType(_curImg, "Mitochondria"), "ImageSummary");
    Table.set("Golgi_count", r, countType(_curImg, "Golgi"), "ImageSummary");
    Table.set("Vacuole_count", r, countType(_curImg, "Vacuole"), "ImageSummary");
    Table.set("Lipid_count", r, countType(_curImg, "LipidBody"), "ImageSummary");
    Table.set("Pore_count", r, poreCount(_curImg), "ImageSummary");
    Table.update("ImageSummary");
}

function saveAll(dir) {
    saveTab("Nuclei", dir); saveTab("Organelles", dir); saveTab("Micronuclei", dir);
    saveTab("Membrane", dir); saveTab("Pores", dir); saveTab("ImageSummary", dir);
}

// ---- distance map from the nucleus ROI ----
function buildDistMap() {
    if (selectionType() < 0) return;
    orig = getTitle(); w = getWidth(); h = getHeight();
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    newImage(DMAP, "8-bit black", w, h, 1);
    selectWindow(orig); roiManager("reset"); roiManager("add");
    selectWindow(DMAP); roiManager("select", 0); setColor(255); fill(); run("Select None");
    setOption("BlackBackground", true); run("Invert"); run("Distance Map");
    HAVE_DMAP = true; selectWindow(orig); roiManager("reset");
}
function distEdge(cx, cy) {
    if (!HAVE_DMAP || !isOpen(DMAP)) return NaN;
    cur = getTitle(); selectWindow(DMAP); getPixelSize(u, pw, ph);
    px = round(cx / pw); py = round(cy / ph);
    if (px < 0) px = 0; if (py < 0) py = 0;
    if (px >= getWidth()) px = getWidth() - 1;
    if (py >= getHeight()) py = getHeight() - 1;
    d = getPixel(px, py) * pw; selectWindow(cur); return d;
}

// ---- file walking ----
function listTiffs(dir) {
    out = newArray(0);
    l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        if (endsWith(l[i], "/")) out = Array.concat(out, listTiffs(dir + l[i]));
        else { n = toLowerCase(l[i]); if (endsWith(n, ".tif") || endsWith(n, ".tiff")) out = Array.concat(out, dir + l[i]); }
    }
    return out;
}

// ---- helpers ----
function initTables() {
    mk("Nuclei"); mk("Organelles"); mk("Micronuclei"); mk("Membrane"); mk("Pores"); mk("ImageSummary");
}
function mk(name) { if (!isOpen(name)) Table.create(name); }
function setMeas() { run("Set Measurements...", "area mean centroid perimeter shape redirect=None decimal=3"); }
function countImg(tbl, img) {
    n = Table.size(tbl); c = 0;
    for (i = 0; i < n; i++) if (Table.getString("Image", i, tbl) == img) c = c + 1;
    return c;
}
function countType(img, otype) {
    n = Table.size("Organelles"); c = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, "Organelles") == img && Table.getString("Type", i, "Organelles") == otype) c = c + 1;
    return c;
}
function nucPerim(img) {
    n = Table.size("Nuclei"); s = 0;
    for (i = 0; i < n; i++) if (Table.getString("Image", i, "Nuclei") == img) s = s + Table.get("Perimeter_um", i, "Nuclei");
    return s;
}
function poreCount(img) {
    n = Table.size("Pores"); c = 0;
    for (i = 0; i < n; i++) if (Table.getString("Image", i, "Pores") == img) c = c + Table.get("PoreCount", i, "Pores");
    return c;
}
function folderOf() {
    p = split(_curDir, "/\\"); n = p.length;
    if (n >= 1) return p[n-1];
    return _curDir;
}
function saveTab(name, dir) {
    if (isOpen(name) && Table.size(name) > 0) { selectWindow(name); Table.save(dir + name + ".csv"); }
}
