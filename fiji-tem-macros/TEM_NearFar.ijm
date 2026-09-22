// ============================================================================
//  TEM_NearFar.ijm  --  organelle NEAR/FAR + NUCLEAR PORES, in one pass.
//
//  Redo of the two readouts that need the nucleus as a reference. Your organelle
//  AREA, count and circularity from the first pass are fine and are NOT redone.
//
//  Per image you get a menu:
//    Near/far   -> TRACE the nucleus (freehand, follow the envelope), then just
//                  CLICK each organelle (multipoint). Each click's distance to
//                  the nuclear edge -> near / far.  (No oval needed = fast.)
//    Pores      -> TRACE ALONG a stretch of nuclear envelope (freehand LINE),
//                  then CLICK each pore on it. -> count + pores per um.
//    NEXT       -> go to the next image.
//    QUIT       -> save and stop.
//  You can do near/far AND pores on the same image before NEXT.
//  If no nucleus edge is in the frame, just hit NEXT (both need a nucleus).
//
//  Output: TEM_NearFar.csv, saved after every entry (Esc / QUIT safe).
//    RecordType "OrganelleNF": Type, Dist_to_nucleus_um, NearFar
//    RecordType "Pore":        Perimeter_um (envelope length), PoreCount, Pores_per_um
//  TEM_Analyze.ijm reads both. Also saves each traced nucleus in ROI_nucleus/.
//
//  Run all sample folders (or one at a time) -- the Sample column is your folder
//  name, so samples stay separated in the output.
// ============================================================================

var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;

var TYPES = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody", "Lysosome");
var NEAR_FAR_UM = 1.0;     // organelle <= this from the nucleus edge = "near"
var TBL   = "TEM_NearFar";
var DMAP  = "__nucDistNF";
var REF_PW = 0;
var _curImg = "";
var _curDir = "";
var _outDir = "";

inDir   = getDirectory("Choose the FOLDER of high-mag images (organelles + nucleus edge)");
_outDir = getDirectory("Choose an OUTPUT folder for TEM_NearFar.csv");

csvPath = _outDir + "TEM_NearFar.csv";
if (File.exists(csvPath)) {
    if (getBoolean("Found an existing TEM_NearFar.csv here.\n \nYES = RESUME (keep it, add to it)\nNO = start NEW (old one renamed to _prev)")) {
        open(csvPath); TBL = "TEM_NearFar.csv";
    } else {
        File.rename(csvPath, _outDir + "TEM_NearFar_prev.csv");
        if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
        Table.create(TBL);
    }
} else {
    if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
    Table.create(TBL);
}

Dialog.create("Settings");
Dialog.addNumber("Near/far cutoff (um from nuclear edge):", NEAR_FAR_UM);
Dialog.show();
NEAR_FAR_UM = Dialog.getNumber();

files = listTiffs(inDir);
if (files.length == 0) { showMessage("No .tif/.tiff images found."); exit; }

startAt = getNumber("Start from image number (1 = beginning).\nThe title bar shows the number you were on.", 1);
if (startAt < 1) startAt = 1;

setBatchMode(false);
quit = false; done = 0;
for (f = startAt - 1; f < files.length; f++) {
    if (quit) f = files.length;
    else {
        open(files[f]);
        if (bitDepth() == 24 || bitDepth() == 16) run("8-bit");
        _curImg = getTitle();
        _curDir = File.getParent(files[f]);
        REF_PW = 0;
        autoCalibrate();
        run("Enhance Contrast", "saturated=0.35");

        imgDone = false;
        while (!imgDone) {
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg);
            Dialog.addChoice("Action:", newArray(
                "Near/far  (trace nucleus, click organelles)",
                "Nuclear pores  (trace envelope, click pores)",
                "NEXT image  /  SKIP",
                "QUIT and save"), "Near/far  (trace nucleus, click organelles)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Near/far"))      measureNearFar();
            else if (startsWith(a, "Nuclear pores")) measurePores();
            else if (startsWith(a, "NEXT"))          imgDone = true;
            else if (startsWith(a, "QUIT"))          { imgDone = true; quit = true; }
        }
        if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(_outDir + "TEM_NearFar.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + _outDir + "TEM_NearFar.csv");

// ============================================================================
// ---- organelle near/far: trace nucleus once, then click each organelle ----
function measureNearFar() {
    setMeas(); setTool("freehand");
    run("Select None");
    waitForUser("Trace the nucleus",
        "Trace AROUND the nucleus edge in this frame (freehand, follow the\n" +
        "envelope), then OK.  Partial nucleus is fine.");
    if (selectionType() < 0) { showMessage("Nothing traced - near/far skipped."); return; }
    buildDistMap();
    saveNucRoi();
    run("Select None");

    more = true;
    while (more) {
        Dialog.create("Which organelle?");
        Dialog.addChoice("Type:", TYPES, TYPES[1]);
        Dialog.show();
        ty = Dialog.getChoice();
        setTool("multipoint");
        waitForUser("Click each " + ty, "CLICK once on each " + ty + " in this field, then OK.\n(No clicks = none of this type.)");
        cNear = 0; cFar = 0;
        if (selectionType() == 10) {
            getSelectionCoordinates(xs, ys);
            for (k = 0; k < xs.length; k++) {
                d = distEdgePx(xs[k], ys[k]);
                nf = "far"; if (!isNaN(d) && d <= NEAR_FAR_UM) nf = "near";
                if (nf == "near") cNear = cNear + 1; else cFar = cFar + 1;
                writeNF(ty, k + 1, d, nf);
            }
        }
        run("Select None");
        showMessage(ty + ": " + cNear + " near, " + cFar + " far recorded.");
        more = getBoolean("Mark another organelle type in this image?");
    }
}

// ---- nuclear pores: trace the envelope stretch, then click the pores on it ----
function measurePores() {
    setMeas(); setTool("freeline");
    run("Select None");
    waitForUser("Envelope for pores",
        "With the FREEHAND-LINE tool, trace ALONG the stretch of nuclear\n" +
        "envelope where you will count pores (follow the membrane), then OK.\n" +
        "(Partial nucleus is fine - just trace the visible stretch.)");
    envLen = NaN; st = selectionType();
    if (st == 5 || st == 6 || st == 7) envLen = getValue("Length");
    run("Select None");

    setTool("multipoint");
    waitForUser("Nuclear pores", "CLICK each nuclear pore ALONG that stretch, then OK.");
    if (selectionType() != 10) { showMessage("Use the multi-point tool and click the pores."); return; }
    getSelectionCoordinates(xs, ys); count = xs.length;
    dens = NaN; if (!isNaN(envLen) && envLen > 0) dens = count / envLen;
    writePore(envLen, count, dens);
    run("Select None");
    pmsg = count + " pores.";
    if (!isNaN(dens)) pmsg = pmsg + "\nDensity = " + d2s(dens,3) + " pores/um of envelope (over " + d2s(envLen,2) + " um).";
    else pmsg = pmsg + "\nCount saved. (Trace the envelope line first to also get density.)";
    showMessage(pmsg);
}

function writeNF(typ, idx, dist, nf) {
    r = Table.size(TBL);
    smp = folderOf();
    Table.set("Image", r, _curImg, TBL);
    Table.set("Sample", r, smp, TBL);
    Table.set("RecordType", r, "OrganelleNF", TBL);
    Table.set("Type", r, typ, TBL);
    Table.set("Index", r, idx, TBL);
    Table.set("Dist_to_nucleus_um", r, dist, TBL);
    Table.set("NearFar", r, nf, TBL);
    Table.set("Cutoff_um", r, NEAR_FAR_UM, TBL);
    Table.update(TBL);
    selectWindow(TBL); Table.save(_outDir + "TEM_NearFar.csv");
}

function writePore(envLen, count, dens) {
    r = Table.size(TBL);
    smp = folderOf();
    Table.set("Image", r, _curImg, TBL);
    Table.set("Sample", r, smp, TBL);
    Table.set("RecordType", r, "Pore", TBL);
    Table.set("Type", r, "", TBL);
    Table.set("Perimeter_um", r, envLen, TBL);
    Table.set("PoreCount", r, count, TBL);
    Table.set("Pores_per_um", r, dens, TBL);
    Table.update(TBL);
    selectWindow(TBL); Table.save(_outDir + "TEM_NearFar.csv");
}

// ---- 32-bit distance map from the traced nucleus (distance to nuclear edge) --
function buildDistMap() {
    if (selectionType() < 0) return;
    orig = getTitle(); w = getWidth(); h = getHeight();
    getPixelSize(u, pw, ph); REF_PW = pw;
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    if (isOpen("__nucBinNF")) { selectWindow("__nucBinNF"); close(); }
    newImage("__nucBinNF", "8-bit black", w, h, 1);
    selectWindow(orig); roiManager("reset"); roiManager("add");
    selectWindow("__nucBinNF"); roiManager("select", 0); setColor(255); fill(); run("Select None");
    setOption("BlackBackground", true);
    run("Options...", "iterations=1 count=1 black edm=32-bit");
    run("Invert");
    run("Distance Map");
    mapT = getTitle();
    if (mapT != "__nucBinNF" && isOpen("__nucBinNF")) { selectWindow("__nucBinNF"); close(); }
    selectWindow(mapT); rename(DMAP);
    selectWindow(orig); roiManager("reset");
}
// distance (um) from a pixel coordinate to the nuclear edge
function distEdgePx(px, py) {
    if (!isOpen(DMAP) || REF_PW <= 0) return NaN;
    cur = getTitle();
    px = round(px); py = round(py);
    selectWindow(DMAP);
    if (px < 0) px = 0; if (py < 0) py = 0;
    if (px >= getWidth()) px = getWidth() - 1;
    if (py >= getHeight()) py = getHeight() - 1;
    d = getPixel(px, py) * REF_PW;
    selectWindow(cur);
    return d;
}
function saveNucRoi() {
    if (selectionType() < 0) return;
    roiDir = _outDir + "ROI_nucleus/";
    File.makeDirectory(roiDir);
    roiManager("reset"); roiManager("add");
    roiManager("save", roiDir + noExt(_curImg) + "_nuc.roi");
    roiManager("reset");
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
function setMeas() { run("Set Measurements...", "area mean centroid perimeter shape redirect=None decimal=4"); }
function folderOf() { p = split(_curDir, "/\\"); n = p.length; if (n >= 1) return p[n-1]; return _curDir; }
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
