// ============================================================================
//  TEM_Nuclear.ijm  --  SESSION 1: everything about the NUCLEUS.
//  Run this on your LOW-MAG whole-cell images (the ~20 cells).
//
//  Per image, a menu lets you:
//    Nucleus       -> TRACE the nucleus -> circularity, aspect ratio, roundness,
//                     solidity, perimeter, area, heterochromatin % / euchromatin %.
//                     Trace a 2nd/3rd nucleus in the same image = multinucleated.
//    Micronuclei   -> trace each small body -> area, circularity
//    Membrane      -> draw a line across the envelope -> perinuclear thickness (nm)
//    NEXT / SKIP / QUIT
//
//  OUTPUT: one file TEM_Nuclear.csv (+ a QC picture of every nucleus).
//  Saved after EVERY measurement, so Esc / QUIT never loses data.
//  (Organelles + nuclear pores are the SECOND macro: TEM_Organelle.ijm.)
// ============================================================================

// ---- calibration (your pixel sizes, um/pixel) ----
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;

var HET_FACTOR = 0.5;      // heterochromatin = pixels darker than mean-0.5*SD
var TBL   = "TEM_Nuclear";
var _curImg = "";
var _curDir = "";
var _outDir = "";
var _qcDir  = "";

inDir   = getDirectory("Choose the FOLDER of images (low-mag / whole nucleus)");
_outDir = getDirectory("Choose an OUTPUT folder for TEM_Nuclear.csv");
_qcDir  = _outDir + "QC_nucleus/";
File.makeDirectory(_qcDir);

if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
Table.create(TBL);

files = listTiffs(inDir);
if (files.length == 0) { showMessage("No .tif/.tiff images found."); exit; }

setBatchMode(false);
quit = false; done = 0;
for (f = 0; f < files.length; f++) {
    if (quit) f = files.length;
    else {
        open(files[f]);
        if (bitDepth() == 24 || bitDepth() == 16) run("8-bit");
        _curImg = getTitle();
        _curDir = File.getParent(files[f]);
        autoCalibrate();
        run("Enhance Contrast", "saturated=0.35");   // view only

        imgDone = false;
        while (!imgDone) {
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg);
            Dialog.addChoice("Action:", newArray(
                "Nucleus  (trace - shape + chromatin + area + multinucleation)",
                "Micronuclei  (trace each)",
                "Membrane / perinuclear thickness  (line)",
                "NEXT image",
                "SKIP this image",
                "QUIT and save"), "Nucleus  (trace - shape + chromatin + area + multinucleation)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Nucleus"))     actNucleus();
            else if (startsWith(a, "Micronuclei")) actMicronuclei();
            else if (startsWith(a, "Membrane"))    actMembrane();
            else if (startsWith(a, "NEXT"))        imgDone = true;
            else if (startsWith(a, "SKIP"))        imgDone = true;
            else if (startsWith(a, "QUIT"))        { imgDone = true; quit = true; }
        }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(_outDir + "TEM_Nuclear.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + _outDir + "TEM_Nuclear.csv");

// ============================================================================
function writeRow(rec, typ, idx, area, perim, circ, ar, rnd, sol, het, eu, dist, nf, thick, pcnt, pdens) {
    r = Table.size(TBL);
    smp = folderOf();
    Table.set("Image", r, _curImg, TBL);
    Table.set("Sample", r, smp, TBL);
    Table.set("RecordType", r, rec, TBL);
    Table.set("Type", r, typ, TBL);
    Table.set("Index", r, idx, TBL);
    Table.set("Area_um2", r, area, TBL);
    Table.set("Perimeter_um", r, perim, TBL);
    Table.set("Circularity", r, circ, TBL);
    Table.set("AspectRatio", r, ar, TBL);
    Table.set("Roundness", r, rnd, TBL);
    Table.set("Solidity", r, sol, TBL);
    Table.set("Heterochrom_pct", r, het, TBL);
    Table.set("Euchrom_pct", r, eu, TBL);
    Table.set("Dist_to_nucleus_um", r, dist, TBL);
    Table.set("NearFar", r, nf, TBL);
    Table.set("Thickness_nm", r, thick, TBL);
    Table.set("PoreCount", r, pcnt, TBL);
    Table.set("Pores_per_um", r, pdens, TBL);
    Table.update(TBL);
    selectWindow(TBL); Table.save(_outDir + "TEM_Nuclear.csv");
}

function actNucleus() {
    setMeas(); setTool("freehand");
    run("Select None");
    waitForUser("Trace the NUCLEUS",
        "Trace around the nuclear envelope (freehand), then OK.\n(No whole nucleus here? Click OK without tracing.)");
    if (selectionType() < 0) { showMessage("Nothing traced - skipped."); return; }
    roiManager("reset"); roiManager("add");
    selectWindow(_curImg); roiManager("select", 0);
    run("Measure"); m = nResults - 1;
    n = countRec("Nucleus", _curImg) + 1;
    getStatistics(dummyA, nMean, nMin, nMax, nSD);
    thr = round(nMean - HET_FACTOR * nSD);
    getHistogram(hv, hc, 256);
    dark = 0; tot = 0;
    for (k = 0; k < hv.length; k++) { tot = tot + hc[k]; if (hv[k] <= thr) dark = dark + hc[k]; }
    het = NaN; if (tot > 0) het = 100.0 * dark / tot;
    writeRow("Nucleus", "", n,
        getResult("Area", m), getResult("Perim.", m), getResult("Circ.", m),
        getResult("AR", m), getResult("Round", m), getResult("Solidity", m),
        het, 100 - het, NaN, "", NaN, NaN, NaN);
    // QC picture
    selectWindow(_curImg); run("Remove Overlay"); roiManager("select", 0);
    Roi.setStrokeColor("yellow"); Roi.setStrokeWidth(4); run("Add Selection...");
    run("Flatten"); saveAs("PNG", _qcDir + noExt(_curImg) + "_n" + n + "_QC.png"); close();
    selectWindow(_curImg); run("Select None"); roiManager("reset");
    msg = "Nucleus " + n + " recorded.";
    if (n > 1) msg = msg + "\n>> Multinucleated (" + n + " nuclei).";
    showMessage(msg);
}

function actMicronuclei() {
    setMeas(); setTool("freehand");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser("Micronucleus #" + i, "Trace a small body separate from the nucleus, then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); m = nResults - 1;
            writeRow("Micronucleus", "", i, getResult("Area", m), NaN, getResult("Circ.", m),
                NaN, NaN, NaN, NaN, NaN, NaN, "", NaN, NaN, NaN);
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
        waitForUser("Thickness " + (i+1) + " of " + n, "Draw a short line ACROSS the membrane / perinuclear space, then OK.");
        if (selectionType() == 5) {
            lenUm = getValue("Length");
            writeRow("Membrane", "", i+1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "", lenUm*1000.0, NaN, NaN);
            i = i + 1;
        }
    }
    run("Select None"); showMessage(n + " thickness measurements recorded.");
}

// ---- calibration + helpers ----
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
function listTiffs(dir) {
    out = newArray(0); l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        if (endsWith(l[i], "/")) out = Array.concat(out, listTiffs(dir + l[i]));
        else { n = toLowerCase(l[i]); if (endsWith(n, ".tif") || endsWith(n, ".tiff")) out = Array.concat(out, dir + l[i]); }
    }
    return out;
}
function setMeas() { run("Set Measurements...", "area mean standard centroid perimeter shape redirect=None decimal=4"); }
function countRec(rec, img) {
    n = Table.size(TBL); c = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, TBL) == img && Table.getString("RecordType", i, TBL) == rec) c = c + 1;
    return c;
}
function folderOf() {
    p = split(_curDir, "/\\"); n = p.length;
    if (n >= 1) return p[n-1];
    return _curDir;
}
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
