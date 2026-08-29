// ============================================================================
//  TEM_PeripheralChromatin.ijm  --  peripheral vs interior heterochromatin.
//  NO re-tracing: it reuses the yellow nucleus outlines in QC_nucleus.
//
//  For each nucleus it splits the area into:
//    PERIPHERY = a band of set width just inside the nuclear envelope
//    INTERIOR  = everything deeper in
//  and measures heterochromatin (dark pixels) in each zone. This tests whether
//  heterochromatin sits at the edge (lamin-tethered) or has moved inward.
//
//  Reports, per nucleus:
//    Periph_Het_pct     % dark (heterochromatin) in the peripheral band
//    Interior_Het_pct   % dark in the interior
//    Periph_MeanGrey / Interior_MeanGrey   (lower = darker = more heterochromatin)
//    PeriphEnrichment   Periph_Het_pct / Interior_Het_pct
//                       ( >1 = heterochromatin concentrated at the edge;
//                         ~1  = evenly dispersed;  <1 = more in the middle )
//
//  A lower PeriphEnrichment in the mutant = heterochromatin has left the edge.
//  Output: TEM_PeripheralChromatin.csv  (+ a QC picture: cyan band = periphery).
//
//  NOTE: uses a FIXED grey cutoff -> only fair if imaging contrast was
//  consistent across samples. MeanGrey columns are the cutoff-free backup.
// ============================================================================

var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;

var TBL = "TEM_PeripheralChromatin";
var BAND_UM   = 0.5;      // width of the peripheral band (microns)
var FIXED_THR = 100;      // grey cutoff: pixels <= this = heterochromatin
var _origDir = "";
var _sample  = "";
var _outDir  = "";
var _qcDir   = "";

_origDir = getDirectory("Choose the folder of ORIGINAL images (.tif) for ONE sample");
qcDir    = getDirectory("Choose that sample's QC_nucleus folder (yellow outlines)");
_outDir  = getDirectory("Choose an OUTPUT folder");
_qcDir   = _outDir + "QC_periphery/";
File.makeDirectory(_qcDir);
_sample = folderOf(_origDir);

Dialog.create("Peripheral chromatin settings");
Dialog.addNumber("Peripheral band width (microns):", BAND_UM);
Dialog.addNumber("Heterochromatin cutoff (grey 0-255):", FIXED_THR);
Dialog.show();
BAND_UM = Dialog.getNumber();
FIXED_THR = Dialog.getNumber();

if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
Table.create(TBL);
run("Set Measurements...", "area mean redirect=None decimal=4");
setBatchMode(true);

list = getFileList(qcDir);
n = 0;
for (i = 0; i < list.length; i++) {
    nm = list[i];
    if (endsWith(toLowerCase(nm), "_qc.png")) {
        base = replace(nm, "_n[0-9]+_[qQ][cC]\\.png$", "");
        if (base == nm) base = replace(nm, "_[qQ][cC]\\.png$", "");
        if (process(qcDir + nm, base, nm)) n = n + 1;
    }
}
setBatchMode(false);
selectWindow(TBL);
Table.save(_outDir + "TEM_PeripheralChromatin.csv");
showMessage("Done", "Measured " + n + " nuclei.\nSaved: " + _outDir + "TEM_PeripheralChromatin.csv\nBand QC: " + _qcDir);

// ============================================================================
function process(qcPath, base, qcName) {
    open(qcPath); qcTitle = getTitle();
    w = getWidth(); h = getHeight();
    ok = yellowToRoi(qcTitle);           // nucleus ROI -> roiManager[0]
    closeAllImages();
    if (!ok) { print("no outline: " + base); return false; }

    op = _origDir + base + ".tif";
    if (!File.exists(op)) op = _origDir + base + ".tiff";
    if (!File.exists(op)) { print("no original: " + base); roiManager("reset"); return false; }

    open(op); orig = getTitle();
    if (bitDepth() == 24 || bitDepth() == 16) run("8-bit");
    if (getWidth() != w || getHeight() != h) { print("size mismatch: " + base); close(); roiManager("reset"); return false; }
    autoCalibrate();
    getPixelSize(u, pw, ph);
    bandPx = round(BAND_UM / pw);
    if (bandPx < 1) bandPx = 1;
    if (bandPx > 254) bandPx = 254;

    // distance-from-edge map of the nucleus interior
    newImage("nz", "8-bit black", w, h, 1);
    selectWindow("nz"); roiManager("select", 0); setColor(255); fill(); run("Select None");
    setOption("BlackBackground", true); run("Distance Map");   // nz now = distance to edge

    // peripheral band (1..bandPx) and interior (>bandPx) as ROIs
    selectWindow("nz"); setThreshold(1, bandPx); run("Create Selection");
    roiManager("add");                    // index 1 = periphery
    run("Select None"); resetThreshold();
    selectWindow("nz"); setThreshold(bandPx + 1, 255); run("Create Selection");
    roiManager("add");                    // index 2 = interior
    run("Select None"); resetThreshold();
    if (isOpen("nz")) { selectWindow("nz"); close(); }

    // measure on the ORIGINAL image
    selectWindow(orig);
    roiManager("select", 1); getStatistics(pArea, pMean); pHet = darkFrac(FIXED_THR);
    roiManager("select", 2); getStatistics(iArea, iMean); iHet = darkFrac(FIXED_THR);

    enrich = NaN; if (!isNaN(iHet) && iHet > 0) enrich = pHet / iHet;

    r = Table.size(TBL);
    Table.set("Image", r, base, TBL);
    Table.set("Sample", r, _sample, TBL);
    Table.set("BandWidth_um", r, BAND_UM, TBL);
    Table.set("Cutoff", r, FIXED_THR, TBL);
    Table.set("Periph_Het_pct", r, pHet, TBL);
    Table.set("Interior_Het_pct", r, iHet, TBL);
    Table.set("Periph_MeanGrey", r, pMean, TBL);
    Table.set("Interior_MeanGrey", r, iMean, TBL);
    Table.set("PeriphEnrichment", r, enrich, TBL);
    Table.update(TBL);
    selectWindow(TBL); Table.save(_outDir + "TEM_PeripheralChromatin.csv");

    // QC: cyan = peripheral band on the original
    selectWindow(orig); run("Remove Overlay");
    roiManager("select", 1); Roi.setStrokeColor("cyan"); Roi.setStrokeWidth(3); run("Add Selection...");
    run("Flatten"); saveAs("PNG", _qcDir + noExt(qcName) + "_band.png"); close();

    selectWindow(orig); close(); roiManager("reset");
    return true;
}

// % of current selection darker than thr
function darkFrac(thr) {
    getHistogram(hv, hc, 256);
    d = 0; t = 0;
    for (k = 0; k < hv.length; k++) { t = t + hc[k]; if (hv[k] <= thr) d = d + hc[k]; }
    if (t > 0) return 100.0 * d / t;
    return NaN;
}

// extract yellow outline -> filled region -> ROI in roiManager[0]
function yellowToRoi(t) {
    selectWindow(t);
    if (bitDepth() != 24) run("RGB Color");
    run("Split Channels");
    rID = t + " (red)"; gID = t + " (green)"; bID = t + " (blue)";
    selectWindow(rID); setThreshold(180, 255); setOption("BlackBackground", true); run("Convert to Mask");
    selectWindow(bID); setThreshold(0, 110); run("Convert to Mask");
    imageCalculator("AND create", rID, bID);
    res = "Result of " + rID;
    if (!isOpen(res)) return false;
    selectWindow(res); run("Fill Holes");
    getStatistics(a2, m2);
    if (m2 == 0) return false;
    run("Create Selection");
    roiManager("reset"); roiManager("add");
    return true;
}

function autoCalibrate() {
    mag = readMag(getTitle());
    if (mag <= 0) return false;
    px = pxForMag(mag);
    setVoxelSize(px, px, 1, "micron");
    return true;
}
function readMag(name) {
    s = name;
    if (matches(s, ".*[0-9]+[ ]?[kK]?[ ]?[xX].*")) {
        num = replace(s, ".*?([0-9]+)[ ]?([kK]?)[ ]?[xX].*", "$1");
        kf  = replace(s, ".*?([0-9]+)[ ]?([kK]?)[ ]?[xX].*", "$2");
        v = parseFloat(num); if (kf == "k" || kf == "K") v = v * 1000; return v;
    }
    return -1;
}
function pxForMag(mag) {
    for (i = 0; i < CAL_MAG.length; i++) if (CAL_MAG[i] == mag) return CAL_PX[i];
    if (mag > 0) return CAM_CONST / mag;
    return -1;
}
function closeAllImages() { while (nImages > 0) { selectImage(nImages); close(); } }
function folderOf(dir) { p = split(dir, "/\\"); k = p.length; if (k >= 1) return p[k-1]; return dir; }
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
