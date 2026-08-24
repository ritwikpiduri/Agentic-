// ============================================================================
//  TEM_Complete.ijm  --  ONE code, everything measured, ONE CSV.
//
//  It opens each image in a folder, auto-calibrates from the filename
//  magnification, and shows a small menu. Pick what you want to do on that
//  image; you can do several things on one image before moving on.
//
//  RELIABLE + FAST (1 click):
//    Nucleus       -> click ONCE inside the nucleus. It grows the outline and
//                     measures circularity, AR, roundness, solidity, perimeter,
//                     area, heterochromatin % and euchromatin %.
//                     Click a 2nd/3rd nucleus in the same image = multinucleated.
//
//  YOU MARK, CODE MEASURES (only on images that need it):
//    Micronuclei   -> trace each small body -> area, circularity
//    Organelles    -> pick type, trace each -> area, circularity, distance to
//                     the nucleus, and near/far
//    Membrane      -> draw a line across it -> thickness in nm
//    Nuclear pores -> click each pore -> count + pores per micron of envelope
//
//  Do the Nucleus first: it builds the distance map that gives organelle
//  near/far and the envelope length for pore density.
//
//  OUTPUT: one file TEM_Results.csv (one row per measurement, tagged by
//  RecordType) + a QC picture of each nucleus. The CSV is saved after EVERY
//  measurement, so you can press Esc / QUIT any time and lose nothing.
// ============================================================================

// ---- calibration (your pixel sizes, um/pixel) ----
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;      // other mags: um/px = 14.4 / magnification

var HET_FACTOR  = 0.5;     // heterochromatin = pixels darker than mean-0.5*SD
var TOL         = 25;      // nucleus grow sensitivity (higher grabs more)
var BLUR        = 10;      // smoothing before growing
var NEAR_FAR_UM = 1.0;     // organelle <= this from nucleus edge = "near"

var TBL   = "TEM_Results";
var DMAP  = "__nucDist";
var HAVE_DMAP = false;
var _curImg = "";
var _curDir = "";
var _outDir = "";
var _qcDir  = "";

inDir   = getDirectory("Choose the FOLDER of images to analyse");
_outDir = getDirectory("Choose an OUTPUT folder for TEM_Results.csv");
_qcDir  = _outDir + "QC_nucleus/";
File.makeDirectory(_qcDir);

if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
Table.create(TBL);

Dialog.create("Nucleus grow settings");
Dialog.addNumber("Sensitivity (higher grabs more):", TOL);
Dialog.addNumber("Smoothing (blur, pixels):", BLUR);
Dialog.addNumber("Organelle near/far cutoff (um):", NEAR_FAR_UM);
Dialog.show();
TOL = Dialog.getNumber(); BLUR = Dialog.getNumber(); NEAR_FAR_UM = Dialog.getNumber();

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
        HAVE_DMAP = false;
        autoCalibrate();
        run("Enhance Contrast", "saturated=0.35");   // view only

        imgDone = false;
        while (!imgDone) {
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg);
            Dialog.addChoice("Action:", newArray(
                "Nucleus  (1 click - shape + chromatin + multinucleation)",
                "Nucleus edge for distance only  (partial nucleus, high-mag)",
                "Micronuclei  (trace each)",
                "Organelles  (type, size, near/far)",
                "Membrane thickness  (line)",
                "Nuclear pores  (click each)",
                "NEXT image",
                "SKIP this image",
                "QUIT and save"), "Nucleus  (1 click - shape + chromatin + multinucleation)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Nucleus edge")) actNucleusEdge();
            else if (startsWith(a, "Nucleus"))      actNucleus();
            else if (startsWith(a, "Micronuclei"))  actMicronuclei();
            else if (startsWith(a, "Organelles"))   actOrganelles();
            else if (startsWith(a, "Membrane"))     actMembrane();
            else if (startsWith(a, "Nuclear pores")) actPores();
            else if (startsWith(a, "NEXT"))         imgDone = true;
            else if (startsWith(a, "SKIP"))         imgDone = true;
            else if (startsWith(a, "QUIT"))         { imgDone = true; quit = true; }
        }
        if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(_outDir + "TEM_Results.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + _outDir + "TEM_Results.csv");

// ============================================================================
// one row per measurement, all columns (NaN / "" where not applicable)
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
    selectWindow(TBL); Table.save(_outDir + "TEM_Results.csv");   // save after every row
}

// ---- Nucleus: one click -> grow -> measure ----
function actNucleus() {
    setMeas(); setTool("point");
    run("Select None");
    waitForUser("Nucleus",
        "CLICK once inside the nucleus, then OK.\n(No whole nucleus here? Click OK without clicking.)");
    if (selectionType() != 10) { showMessage("Nothing clicked - skipped."); return; }
    getSelectionCoordinates(xs, ys); xc = xs[0]; yc = ys[0];
    run("Select None");

    run("Duplicate...", "title=nwork");
    run("Gaussian Blur...", "sigma=" + BLUR);
    doWand(xc, yc, TOL, "8-connected");
    if (selectionType() < 0) {
        if (isOpen("nwork")) { selectWindow("nwork"); close(); }
        showMessage("Could not grow there. Try higher Sensitivity or click nearer the centre.");
        return;
    }
    getStatistics(gArea);
    w = getWidth(); h = getHeight(); getPixelSize(u2, pw2, ph2);
    frameA = w * h * pw2 * ph2;
    if (gArea > 0.9 * frameA) {
        if (isOpen("nwork")) { selectWindow("nwork"); close(); }
        showMessage("That grew to almost the whole image. Lower Sensitivity or click deeper inside the nucleus.");
        return;
    }
    roiManager("reset"); roiManager("add");
    if (isOpen("nwork")) { selectWindow("nwork"); close(); }

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

    // distance map for organelle near/far + pore density
    selectWindow(_curImg); roiManager("select", 0); buildDistMap();
    run("Select None");
    msg = "Nucleus " + n + " recorded.";
    if (n > 1) msg = msg + "\n>> Multinucleated (" + n + " nuclei).";
    showMessage(msg);
}

// ---- Nucleus edge: trace the partial nucleus for near/far ONLY (no shape row) ----
function actNucleusEdge() {
    setMeas(); setTool("freehand");
    run("Select None");
    waitForUser("Nucleus edge (distance only)",
        "Trace the part of the NUCLEUS shown in this image (the nuclear side),\n" +
        "then OK. This sets up organelle near/far.\nNo nucleus shape is recorded.");
    if (selectionType() < 0) { showMessage("Nothing traced - skipped."); return; }
    buildDistMap(); run("Select None");
    showMessage("Nucleus edge set for this image.\nNow trace organelles to get their near/far.");
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

function actOrganelles() {
    types = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody");
    Dialog.create("Organelle type");
    Dialog.addChoice("Marking which organelle?", types, types[0]);
    Dialog.show();
    otype = Dialog.getChoice();
    if (!HAVE_DMAP) showMessage("Tip: do the Nucleus first so distances (near/far) are recorded.");
    setMeas(); setTool("freehand");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser(otype + " #" + i, "Trace this " + otype + " (freehand), then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); m = nResults - 1;
            cx = getResult("X", m); cy = getResult("Y", m);
            dist = NaN; nf = "";
            if (HAVE_DMAP) { dist = distEdge(cx, cy); if (dist <= NEAR_FAR_UM) nf = "near"; else nf = "far"; }
            writeRow("Organelle", otype, i, getResult("Area", m), NaN, getResult("Circ.", m),
                NaN, NaN, NaN, NaN, NaN, dist, nf, NaN, NaN, NaN);
        }
        more = getBoolean("Mark another " + otype + "?");
    }
    run("Select None"); showMessage(i + " " + otype + "(s) recorded.");
}

function actMembrane() {
    setTool("line");
    n = getNumber("How many thickness measurements on this image?", 5);
    i = 0;
    while (i < n) {
        waitForUser("Thickness " + (i+1) + " of " + n, "Draw a short line ACROSS the membrane, then OK.");
        if (selectionType() == 5) {
            lenUm = getValue("Length");
            writeRow("Membrane", "", i+1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "", lenUm*1000.0, NaN, NaN);
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
    writeRow("Pore", "", 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "", NaN, count, dens);
    run("Select None");
    pmsg = count + " pores.";
    if (perim > 0) pmsg = pmsg + "\nDensity = " + d2s(dens,3) + " pores/um.";
    else pmsg = pmsg + "\n(Do the Nucleus first for density.)";
    showMessage(pmsg);
}

// ---- distance map from the current nucleus selection ----
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

// ---- calibration ----
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

// ---- helpers ----
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
function nucPerim(img) {
    n = Table.size(TBL); s = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, TBL) == img && Table.getString("RecordType", i, TBL) == "Nucleus")
            s = s + Table.get("Perimeter_um", i, TBL);
    return s;
}
function folderOf() {
    p = split(_curDir, "/\\"); n = p.length;
    if (n >= 1) return p[n-1];
    return _curDir;
}
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
