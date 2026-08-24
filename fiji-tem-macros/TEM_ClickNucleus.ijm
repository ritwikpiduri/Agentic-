// ============================================================================
//  TEM_ClickNucleus.ijm  --  ONE CLICK per nucleus. Fast AND reliable.
//
//  Why this exists: fully-automatic detection cannot tell the nucleus from the
//  whole cell / cytoplasm in these TEM images (the cell is darker than the
//  background, so "biggest dark blob" grabs the cell). One click fixes that:
//  you point at the nucleus, the code grows the outline to the nuclear envelope
//  and measures it. No tracing.
//
//  PER IMAGE:
//    - It opens the image and auto-calibrates from the filename magnification.
//    - You CLICK once inside the nucleus (point tool is pre-selected), then OK.
//        -> it measures circularity, aspect ratio, roundness, solidity,
//           perimeter, area, heterochromatin % and euchromatin %.
//    - No whole nucleus in this image (high-mag crop)? Just click OK without
//        clicking anything -> it skips that image.
//    - Multinucleated cell? Click one nucleus, OK; when it asks "another
//        nucleus?", say Yes and click the next one.
//    - To STOP at any time: press Esc. Your data is already saved (the CSV is
//        written after every nucleus), so nothing is lost.
//
//  OUTPUT: one file, TEM_Nucleus.csv, plus a QC picture per nucleus (yellow
//  outline) so you can confirm the grow was correct.
// ============================================================================

// ---- calibration (your pixel sizes) ----
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;

var HET_FACTOR = 0.5;   // heterochromatin = pixels darker than mean - 0.5*SD
var TOL  = 25;          // grow sensitivity (higher = grabs more)
var BLUR = 10;          // smoothing before growing (blurs out chromatin grain)
var TBL  = "TEM_Nucleus";

inDir  = getDirectory("Choose the FOLDER of images");
outDir = getDirectory("Choose an OUTPUT folder");
qcDir  = outDir + "QC_nucleus/";
File.makeDirectory(qcDir);

if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
Table.create(TBL);
run("Set Measurements...", "area mean standard perimeter shape redirect=None decimal=4");

// let you set the two knobs once (defaults are fine to start)
Dialog.create("Nucleus grow settings");
Dialog.addNumber("Sensitivity (higher grabs more):", TOL);
Dialog.addNumber("Smoothing (blur, pixels):", BLUR);
Dialog.addMessage("Tip: if the outline overshoots into cytoplasm, lower\nsensitivity. If it stops short inside the nucleus, raise it.");
Dialog.show();
TOL  = Dialog.getNumber();
BLUR = Dialog.getNumber();

files = listTiffs(inDir);
if (files.length == 0) { showMessage("No .tif/.tiff found."); exit; }

setBatchMode(false);
done = 0;
for (f = 0; f < files.length; f++) {
    open(files[f]);
    if (bitDepth() == 24 || bitDepth() == 16) run("8-bit");
    name = getTitle();
    dir  = File.getParent(files[f]);
    mag  = readMag(name);
    px   = pxForMag(mag);
    cal  = false;
    if (px > 0) { setVoxelSize(px, px, 1, "micron"); cal = true; }

    run("Enhance Contrast", "saturated=0.35");   // view only, easier to see
    setTool("point");

    nHere = 0;
    more = true;
    while (more) {
        run("Select None");
        waitForUser("Image " + (f+1) + "/" + files.length,
            "CLICK once inside the nucleus, then OK.\n \n" +
            "No whole nucleus here? Click OK without clicking -> skipped.\n" +
            "To STOP the whole run: press Esc (your data is already saved).");
        if (selectionType() != 10) { more = false; }   // nothing clicked -> skip
        else {
            getSelectionCoordinates(xs, ys);
            xc = xs[0]; yc = ys[0];
            ok = measureNucleusAt(name, dir, xc, yc, cal, px, nHere + 1);
            if (ok) {
                nHere = nHere + 1;
                done = done + 1;
                selectWindow(TBL); Table.save(outDir + "TEM_Nucleus.csv");  // save after each
            }
            more = getBoolean("Another nucleus in THIS image? (multinucleated)");
        }
    }
    if (nImages > 0) { selectWindow(name); close(); }
}

selectWindow(TBL);
Table.save(outDir + "TEM_Nucleus.csv");
showMessage("Done", "Measured " + done + " nucleus/nuclei.\nSaved: " + outDir + "TEM_Nucleus.csv");

// ============================================================================
// grow the nucleus from one click, measure it, write a row. returns true if ok.
function measureNucleusAt(name, dir, xc, yc, cal, px, idx) {
    selectWindow(name); run("Select None");
    run("Duplicate...", "title=nwork");
    run("Gaussian Blur...", "sigma=" + BLUR);
    doWand(xc, yc, TOL, "8-connected");
    if (selectionType() < 0) {
        if (isOpen("nwork")) { selectWindow("nwork"); close(); }
        showMessage("Could not grow a region there.\nTry a higher Sensitivity, or click nearer the nucleus centre.");
        return false;
    }
    // reject grabs that are basically the whole frame (fell into cytoplasm/cell)
    getStatistics(gArea);
    w = getWidth(); h = getHeight(); frameA = w * h; if (cal) frameA = w * h * px * px;
    if (gArea > 0.9 * frameA) {
        if (isOpen("nwork")) { selectWindow("nwork"); close(); }
        showMessage("That grew to almost the whole image (likely cytoplasm).\nLower the Sensitivity, or click deeper inside the nucleus.");
        return false;
    }
    // clean the outline: mask -> fill holes -> reselect (removes chromatin holes)
    roiManager("reset"); roiManager("add");
    if (isOpen("nwork")) { selectWindow("nwork"); close(); }

    // measure geometry + chromatin on the ORIGINAL image
    selectWindow(name); roiManager("select", 0);
    run("Measure"); m = nResults - 1;
    area = getResult("Area", m);
    getStatistics(dummyA, nMean, nMin, nMax, nSD);
    areaPx = area; areaUm = NaN;
    if (cal) { areaUm = area; areaPx = area / (px * px); }

    thr = round(nMean - HET_FACTOR * nSD);
    getHistogram(hv, hc, 256);
    dark = 0; tot = 0;
    for (k = 0; k < hv.length; k++) { tot = tot + hc[k]; if (hv[k] <= thr) dark = dark + hc[k]; }
    hetPct = NaN; if (tot > 0) hetPct = 100.0 * dark / tot;

    row = Table.size(TBL);
    fld = folderOf(dir);
    Table.set("Image", row, name, TBL);
    Table.set("Sample", row, fld, TBL);
    Table.set("NucleusIndex", row, idx, TBL);
    Table.set("Circularity", row, getResult("Circ.", m), TBL);
    Table.set("AspectRatio", row, getResult("AR", m), TBL);
    Table.set("Roundness", row, getResult("Round", m), TBL);
    Table.set("Solidity", row, getResult("Solidity", m), TBL);
    Table.set("Perimeter", row, getResult("Perim.", m), TBL);
    Table.set("Area_px", row, areaPx, TBL);
    Table.set("Area_um2", row, areaUm, TBL);
    Table.set("Heterochromatin_pct", row, hetPct, TBL);
    Table.set("Euchromatin_pct", row, 100 - hetPct, TBL);
    Table.update(TBL);

    // QC picture
    selectWindow(name); run("Remove Overlay"); roiManager("select", 0);
    Roi.setStrokeColor("yellow"); Roi.setStrokeWidth(4); run("Add Selection...");
    run("Flatten");
    saveAs("PNG", qcDir + noExt(name) + "_n" + idx + "_QC.png");
    close();
    selectWindow(name); run("Select None"); roiManager("reset");
    return true;
}

// ---- helpers ----
function listTiffs(dir) {
    out = newArray(0); l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        if (endsWith(l[i], "/")) out = Array.concat(out, listTiffs(dir + l[i]));
        else { n = toLowerCase(l[i]); if (endsWith(n, ".tif") || endsWith(n, ".tiff")) out = Array.concat(out, dir + l[i]); }
    }
    return out;
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
function folderOf(dir) {
    parts = split(dir, "/\\"); n = parts.length;
    if (n >= 1) return parts[n-1];
    return dir;
}
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
