// ============================================================================
//  TEM_Chromatin_FromQC.ijm  --  re-measure chromatin WITHOUT re-tracing.
//
//  The nuclear macro already saved a yellow outline for every nucleus in the
//  QC_nucleus folder. This macro re-uses those outlines: it extracts the
//  yellow region, lays it back on the ORIGINAL image, and measures chromatin
//  a better way -- so it can actually detect a real heterochromatin difference
//  (the old per-nucleus threshold self-normalised to ~31% for everyone).
//
//  It reports, per nucleus:
//    MeanGrey            average brightness inside the nucleus (LOWER = darker =
//                        MORE heterochromatin). Threshold-free -> the cleanest.
//    Heterochrom_pct     % of nucleus darker than a FIXED grey cutoff (same for
//                        every nucleus). A darker/more-condensed nucleus reads
//                        HIGHER here.
//    Euchrom_pct         100 - Heterochrom_pct
//
//  Output: TEM_Chromatin.csv (one row per nucleus). Re-running with a different
//  cutoff is instant (no tracing), so you can try a couple of values.
//
//  NOTE: a FIXED cutoff is only fair if imaging contrast was consistent across
//  samples. If unsure, rely on MeanGrey, or tell me and we normalise per image.
// ============================================================================

var TBL = "TEM_Chromatin";
var FIXED_THR = 100;      // grey cutoff (0-255): pixels <= this = heterochromatin
var _origDir = "";
var _sample  = "";
var _outDir  = "";

_origDir = getDirectory("Choose the folder of ORIGINAL images (.tif) for ONE sample");
qcDir    = getDirectory("Choose that sample's QC_nucleus folder (yellow outlines)");
_outDir  = getDirectory("Choose an OUTPUT folder");

_sample = folderOf(_origDir);

Dialog.create("Chromatin threshold");
Dialog.addNumber("Fixed heterochromatin cutoff (grey 0-255):", FIXED_THR);
Dialog.addMessage("Pixels darker than this = heterochromatin, for EVERY nucleus.\n" +
    "Tip: open one nucleus, Image > Adjust > Threshold, slide until the dark\n" +
    "condensed chromatin is just covered, read the upper value, use it here.\n" +
    "MeanGrey is reported too and needs no cutoff.");
Dialog.show();
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
        base = replace(nm, "_n[0-9]+_[qQ][cC]\\.png$", "");   // strip _n1_QC.png
        if (base == nm) base = replace(nm, "_[qQ][cC]\\.png$", "");
        if (process(qcDir + nm, base)) n = n + 1;
    }
}
setBatchMode(false);
selectWindow(TBL);
Table.save(_outDir + "TEM_Chromatin.csv");
showMessage("Done", "Re-measured " + n + " nuclei.\nSaved: " + _outDir + "TEM_Chromatin.csv");

// ============================================================================
function process(qcPath, base) {
    open(qcPath); qcTitle = getTitle();
    w = getWidth(); h = getHeight();
    ok = yellowToRoi(qcTitle);
    closeAllImages();
    if (!ok) { print("no outline found: " + base); return false; }

    op = _origDir + base + ".tif";
    if (!File.exists(op)) op = _origDir + base + ".tiff";
    if (!File.exists(op)) { print("no original for: " + base); roiManager("reset"); return false; }

    open(op);
    if (bitDepth() == 24 || bitDepth() == 16) run("8-bit");
    if (getWidth() != w || getHeight() != h) { print("size mismatch: " + base); close(); roiManager("reset"); return false; }

    roiManager("select", 0);
    getStatistics(area, meanGrey, gmin, gmax, gsd);
    getHistogram(hv, hc, 256);
    dark = 0; tot = 0;
    for (k = 0; k < hv.length; k++) { tot = tot + hc[k]; if (hv[k] <= FIXED_THR) dark = dark + hc[k]; }
    het = NaN; if (tot > 0) het = 100.0 * dark / tot;

    r = Table.size(TBL);
    Table.set("Image", r, base, TBL);
    Table.set("Sample", r, _sample, TBL);
    Table.set("MeanGrey", r, meanGrey, TBL);
    Table.set("Heterochrom_pct", r, het, TBL);
    Table.set("Euchrom_pct", r, 100 - het, TBL);
    Table.set("Cutoff", r, FIXED_THR, TBL);
    Table.update(TBL);
    selectWindow(TBL); Table.save(_outDir + "TEM_Chromatin.csv");

    close(); roiManager("reset");
    return true;
}

// extract the yellow nucleus outline -> filled region -> ROI (in roiManager[0])
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
    selectWindow(res);
    run("Fill Holes");
    getStatistics(a2, m2);
    if (m2 == 0) return false;
    run("Create Selection");
    roiManager("reset"); roiManager("add");
    return true;
}

function closeAllImages() { while (nImages > 0) { selectImage(nImages); close(); } }
function folderOf(dir) {
    p = split(dir, "/\\"); k = p.length;
    if (k >= 1) return p[k-1];
    return dir;
}
