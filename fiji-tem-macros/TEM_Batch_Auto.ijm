// =============================================================================
//  TEM_Batch_Auto.ijm
//  Unattended batch pass over a whole TEM dataset tree for the metrics that can
//  be measured WITHOUT human identification:
//     - nuclear shape (circularity, area, AR, roundness, solidity)
//     - heterochromatin vs euchromatin fraction inside the nucleus
//  It walks every .tif/.tiff under a chosen root (e.g. 4 folders x 20 cells x
//  243 images), auto-segments the nucleus, and writes one CSV row per image.
//
//  IMPORTANT / HONEST LIMITS
//  * Auto nucleus segmentation on TEM is approximate. It thresholds the large,
//    homogeneous, mid-density nuclear region and keeps the biggest object.
//    ALWAYS spot-check with the QC overlays this macro can save (SAVE_QC=true)
//    and drop rows where segmentation clearly failed. Treat this as a first
//    pass; refine borderline cells in the interactive macro.
//  * ER / Golgi / mitochondria / lipid / vacuole are NOT auto-detected here -
//    they have no reliable intensity signature in TEM. Do those in
//    TEM_Analysis_Interactive.ijm, or train a classifier (Weka/Labkit/ilastik)
//    and adapt the "measure from mask" section below.
//
//  CALIBRATION: set PIXEL_SIZE_UM to the real pixel size for this dataset. If
//  different magnifications are mixed, split them into separate runs, or encode
//  pixel size per subfolder and edit pixelSizeFor().
// =============================================================================

// -------------------- USER SETTINGS --------------------
// Your TIFFs carry the scale in their metadata, so leave USE_EMBEDDED_SCALE
// true: each image is measured with its OWN embedded pixel size (this also
// handles mixed magnifications automatically). PIXEL_SIZE_UM is only a FALLBACK
// used for any image that opens uncalibrated (unit = pixel). nm/Angstrom units
// are converted to microns so all outputs are in um / um^2.
var USE_EMBEDDED_SCALE = true;
var PIXEL_SIZE_UM   = 0.005;   // FALLBACK only (used if an image has no scale)
var UNIT            = "micron";
var MIN_NUC_UM2     = 2.0;     // ignore objects smaller than this as nuclei
var MAX_NUC_FRAC    = 0.95;    // ignore objects bigger than 95% of frame (junk)
var HET_STD_FACTOR  = 0.5;     // heterochromatin threshold = mean - factor*std
var SAVE_QC         = true;    // save an outline overlay per image for checking
var QC_SUBDIR       = "QC_overlays";
// -------------------------------------------------------

setBatchMode(true);
root = getDirectory("Choose the ROOT folder of the TEM dataset");
outDir = getDirectory("Choose an OUTPUT folder for results");
if (SAVE_QC) File.makeDirectory(outDir + QC_SUBDIR);
loadCalibration(findCalib(root, outDir));   // magnification -> pixel-size table

if (isOpen("BatchResults")) { selectWindow("BatchResults"); run("Close"); }
Table.create("BatchResults");

run("Set Measurements...",
    "area mean standard min centroid perimeter shape redirect=None decimal=4");

var fileCount = 0;
var doneCount = 0;
var totalTiffs = 0;
var startTime  = 0;
totalTiffs = countTiffs(root);   // pre-count so we can show "X of TOTAL"
startTime  = getTime();
print("\\Clear");
print("=== TEM batch: " + totalTiffs + " TIFF images found under ===");
print(root);
processTree(root);

selectWindow("BatchResults");
Table.save(outDir + "TEM_Batch_NucleusChromatin.csv");
setBatchMode(false);
showMessage("Batch complete",
    "Processed " + doneCount + " / " + fileCount + " images.\n" +
    (fallbackCount > 0 ?
        ">> WARNING: " + fallbackCount + " image(s) had NO embedded scale and\n" +
        "   used the fallback PIXEL_SIZE_UM. Check the ScaleSource column;\n" +
        "   'fallback' rows may have wrong areas/distances.\n" : "") +
    "Results: " + outDir + "TEM_Batch_NucleusChromatin.csv" +
    (SAVE_QC ? "\nQC overlays: " + outDir + QC_SUBDIR : ""));

// =============================================================================
function processTree(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        p = dir + list[i];
        if (endsWith(list[i], "/")) {
            processTree(p);                       // recurse into subfolder
        } else if (isTiff(list[i])) {
            fileCount++;
            processImage(p, dir, list[i]);
        }
    }
}
function isTiff(name) {
    n = toLowerCase(name);
    return endsWith(n, ".tif") || endsWith(n, ".tiff");
}
// Calibration priority: (1) embedded metadata scale, (2) magnification read
// from the filename -> pixel size from the calibration table, (3) fallback.
var fallbackCount = 0;
function applyCalibration() { return applyCalibrationFor(getTitle()); }
function applyCalibrationFor(name) {
    if (USE_EMBEDDED_SCALE && normalizeToMicron()) return "embedded";
    mag = parseMag(name);
    if (mag > 0) {
        px = pxForMag(mag);
        if (px > 0) { setVoxelSize(px, px, 1, "micron"); return "mag=" + mag; }
    }
    setVoxelSize(PIXEL_SIZE_UM, PIXEL_SIZE_UM, 1, "micron");
    fallbackCount++;
    return "fallback";
}
// ---- magnification-from-filename calibration table ----
var CAL_MAG = newArray(0);
var CAL_PX  = newArray(0);
function findCalib(root, outDir) {
    cands = newArray(outDir + "magnification_calibration.csv",
                     root   + "magnification_calibration.csv");
    for (i = 0; i < cands.length; i++) if (File.exists(cands[i])) return cands[i];
    Dialog.create("Calibration table");
    Dialog.addMessage("magnification_calibration.csv not found in the dataset or\n" +
                      "output folder. Build it with TEM_Build_Calibration.ijm.");
    Dialog.addString("Paste full path to it (blank = use fallback pixel size):", "", 60);
    Dialog.show();
    return String.trim(Dialog.getString());
}
function loadCalibration(path) {
    CAL_MAG = newArray(0); CAL_PX = newArray(0);
    if (path == "" || !File.exists(path)) { print("No calibration table loaded."); return false; }
    lines = split(File.openAsString(path), "\n");
    for (i = 0; i < lines.length; i++) {
        ln = String.trim(lines[i]);
        if (ln == "") continue;
        c = split(ln, ",");
        if (c.length < 2) continue;
        a = String.trim(c[0]);
        if (!matches(a, "[0-9].*")) continue;               // skip header
        CAL_MAG = Array.concat(CAL_MAG, parseFloat(a));
        CAL_PX  = Array.concat(CAL_PX, parseFloat(String.trim(c[1])));
    }
    print("Loaded " + CAL_MAG.length + " magnification calibrations from " + path);
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
        v = parseFloat(num);
        if (kfl == "k" || kfl == "K") v = v * 1000;
        return v;
    }
    if (matches(s, ".*[xX][0-9]+.*"))
        return parseFloat(replace(s, ".*[xX]([0-9]+).*", "$1"));
    return -1;
}
// If the image is calibrated in a length unit, convert its scale to microns
// and return true. If it is uncalibrated (unit = pixel) or unknown, return false.
function normalizeToMicron() {
    getPixelSize(unit, pw, ph);
    u = toLowerCase(unit);
    if (u == "micron" || u == "microns" || u == "um" || u == "µm" || u == "microns/pixel")
        return (pw != 1 || u != "pixel");                 // already microns
    factor = 0;
    if (u == "nm" || u == "nanometer" || u == "nanometre" || u == "nanometers")
        factor = 1.0/1000.0;
    else if (u == "a" || u == "angstrom" || u == "ang" || u == "å" || u == "Å")
        factor = 1.0/10000.0;
    else if (u == "mm" || u == "millimeter")
        factor = 1000.0;
    else
        return false;                                     // pixel / unknown
    setVoxelSize(pw * factor, ph * factor, 1, "micron");
    return true;
}
// pre-count all TIFFs under a folder tree (for progress display)
function countTiffs(dir) {
    c = 0; list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) c += countTiffs(dir + list[i]);
        else if (isTiff(list[i])) c++;
    }
    return c;
}

function processImage(path, dir, name) {
    // --- progress line: which folder/cell/image, count, and ETA ---
    elapsed = (getTime() - startTime) / 1000.0;                 // seconds so far
    rate    = (doneCount > 0) ? elapsed / doneCount : 0;        // s per image
    etaMin  = (rate * (totalTiffs - fileCount)) / 60.0;         // minutes left
    showStatus("TEM " + fileCount + "/" + totalTiffs + "  " + name);
    showProgress(fileCount, totalTiffs);
    print("[" + fileCount + "/" + totalTiffs + "]  " + relPath(dir) + "  ->  " +
          name + (rate > 0 ? "   (~" + d2s(etaMin,1) + " min left)" : ""));

    open(path);
    if (bitDepth() == 24) run("8-bit");           // TEM should be grayscale
    calFrom = applyCalibration();                 // "embedded" or "fallback"
    getPixelSize(u, pw, ph);
    w = getWidth(); h = getHeight();
    frameArea = w * h * pw * ph;

    // ---- auto-segment the nucleus ----
    // work on a smoothed duplicate so texture doesn't fragment the region
    orig = getTitle();
    run("Duplicate...", "title=work");
    run("Gaussian Blur...", "sigma=3");
    setAutoThreshold("Otsu dark");                // nucleus tends mid/dark & large
    run("Convert to Mask");
    run("Options...", "iterations=3 count=1 black");
    run("Close-");                                // fill small gaps
    run("Fill Holes");

    // keep the largest plausible particle as the nucleus
    run("Set Measurements...",
        "area mean standard min centroid perimeter shape redirect=None decimal=4");
    roiManager("reset");
    run("Analyze Particles...",
        "size=" + MIN_NUC_UM2 + "-Infinity circularity=0.1-1.00 " +
        "show=Nothing add");
    nRoi = roiManager("count");
    bestRoi = -1; bestArea = 0;
    for (r = 0; r < nRoi; r++) {
        roiManager("select", r);
        getStatistics(a);
        if (a > bestArea && a < MAX_NUC_FRAC * frameArea) { bestArea = a; bestRoi = r; }
    }

    row = Table.size("BatchResults");
    Table.set("Image",  row, name,          "BatchResults");
    Table.set("Folder", row, relPath(dir),  "BatchResults");
    Table.set("PixelSize_um", row, pw,      "BatchResults");
    Table.set("ScaleSource",  row, calFrom, "BatchResults");

    if (bestRoi < 0) {
        Table.set("Status", row, "no_nucleus_found", "BatchResults");
        Table.update("BatchResults");
        cleanup();
        return;
    }

    // measure shape on the ORIGINAL (mask ROI applied to real pixels)
    selectWindow(orig);
    roiManager("select", bestRoi);
    run("Measure");
    m = nResults - 1;
    circ = getResult("Circ.", m);
    area = getResult("Area", m);
    peri = getResult("Perim.", m);
    ar   = getResult("AR", m);
    rnd  = getResult("Round", m);
    sol  = getResult("Solidity", m);
    getStatistics(nucArea, nucMean, nucMin, nucMax, nucStd);

    // ---- chromatin fraction inside the nucleus ----
    thr = round(nucMean - HET_STD_FACTOR * nucStd);
    Roi.getBounds(bx, by, bw, bh);
    hetero = 0; total = 0;
    for (y = by; y < by + bh; y++) {
        for (x = bx; x < bx + bw; x++) {
            if (Roi.contains(x, y)) {
                total++;
                if (getPixel(x, y) <= thr) hetero++;
            }
        }
    }
    hetPct = (total > 0) ? 100.0 * hetero / total : NaN;

    Table.set("Status",          row, "ok",                "BatchResults");
    Table.set("Nucleus_Area_um2",row, area,                "BatchResults");
    Table.set("Perimeter_um",    row, peri,                "BatchResults");
    Table.set("Circularity",     row, circ,                "BatchResults");
    Table.set("AspectRatio",     row, ar,                  "BatchResults");
    Table.set("Roundness",       row, rnd,                 "BatchResults");
    Table.set("Solidity",        row, sol,                 "BatchResults");
    Table.set("Nucleus_MeanGray",row, nucMean,             "BatchResults");
    Table.set("Chromatin_thr",   row, thr,                 "BatchResults");
    Table.set("Heterochrom_pct", row, hetPct,              "BatchResults");
    Table.set("Euchrom_pct",     row, 100 - hetPct,        "BatchResults");
    Table.set("Hetero_dominant", row, (hetPct > 50 ? 1 : 0),"BatchResults");
    Table.update("BatchResults");
    doneCount++;

    // ---- QC overlay ----
    if (SAVE_QC) {
        selectWindow(orig);
        roiManager("select", bestRoi);
        run("Flatten");
        saveAs("PNG", outDir + QC_SUBDIR + "/" + stripExt(name) + "_QC.png");
        close();
    }
    cleanup();
}

function cleanup() {
    roiManager("reset");
    while (nImages > 0) { selectImage(nImages); close(); }
    run("Clear Results");
}
function relPath(dir) {
    // last two path components, for grouping (folder / cell)
    parts = split(dir, "/");
    n = parts.length;
    if (n >= 2) return parts[n-2] + "/" + parts[n-1];
    if (n == 1) return parts[0];
    return dir;
}
function stripExt(name) {
    dot = lastIndexOf(name, ".");
    return (dot > 0) ? substring(name, 0, dot) : name;
}
