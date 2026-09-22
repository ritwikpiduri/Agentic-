// ============================================================================
//  TEM_NearFar.ijm  --  organelle COUNT + DISTANCE-to-nucleus + NUCLEAR PORES,
//                       in ONE pass, built for HIGH-MAG (partial nucleus) images.
//
//  Your organelle AREA / circularity from the first pass are fine and NOT redone.
//  This fixes the two readouts that need the nucleus as a reference.
//
//  Because the nucleus is usually only PARTLY in the frame at high mag, the
//  reference is the VISIBLE MEMBRANE LINE (not a filled shape). Distance is then
//  the true shortest distance from each organelle to that membrane -- correct
//  even when most of the nucleus is off-frame.
//
//  PER IMAGE (one trace does everything):
//    1) With the FREEHAND-LINE tool, trace ALONG the visible nuclear envelope.
//    2) Organelles: pick a type, CLICK each one -> count + distance -> near/far.
//       (Repeat for each type.)
//    3) Pores: CLICK each pore along that envelope -> count + pores per um.
//  If no nucleus edge is visible in a frame, SKIP it (distance is impossible
//  without a membrane to measure to).
//
//  NOTE (partial nucleus): distance is to the membrane you TRACED. An organelle
//  could sit near a piece of membrane that is off-frame; you can only measure to
//  what is visible. That is a limitation of partial images, not a bug.
//
//  Output: TEM_NearFar.csv, saved after every entry (Esc / QUIT safe).
//    RecordType "OrganelleNF": Type, Dist_to_nucleus_um, NearFar
//    RecordType "Pore":        Perimeter_um (envelope length), PoreCount, Pores_per_um
//  TEM_Analyze.ijm reads both. Each traced envelope is saved in ROI_nucleus/.
// ============================================================================

var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;

var TYPES = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody", "Lysosome");
var NEAR_FAR_UM = 1.0;     // organelle <= this from the membrane = "near"
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
Dialog.addNumber("Near/far cutoff (um from the membrane):", NEAR_FAR_UM);
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

        Dialog.create("Image " + (f+1) + " / " + files.length);
        Dialog.addMessage(_curImg);
        Dialog.addChoice("Action:", newArray(
            "Measure  (trace membrane, then organelles + pores)",
            "SKIP  (no nucleus edge in this frame)",
            "QUIT and save"), "Measure  (trace membrane, then organelles + pores)");
        Dialog.show();
        a = Dialog.getChoice();
        if      (startsWith(a, "Measure")) measure();
        else if (startsWith(a, "QUIT"))    quit = true;

        if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(_outDir + "TEM_NearFar.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + _outDir + "TEM_NearFar.csv");

// ============================================================================
// one trace of the visible membrane -> organelle count+distance, then pores
function measure() {
    setMeas(); setTool("freeline");
    run("Select None");
    waitForUser("Trace the nuclear membrane",
        "With the FREEHAND-LINE tool, trace ALONG the visible nuclear\n" +
        "envelope (the membrane) in this frame, then OK.\n" +
        "Partial nucleus is fine - trace whatever membrane you can see.\n" +
        "This one trace is the reference for BOTH distance AND pores.");
    st = selectionType();
    if (!(st == 5 || st == 6 || st == 7)) { showMessage("No line traced - image skipped."); return; }
    envLen = getValue("Length");
    buildDistMap();        // distance-to-membrane-line map (sets REF_PW)
    saveNucRoi();
    run("Select None");

    // --- organelles: count + distance to membrane ---
    if (getBoolean("Mark organelles (count + distance to nucleus) in this image?")) {
        more = true;
        while (more) {
            Dialog.create("Which organelle?");
            Dialog.addChoice("Type:", TYPES, TYPES[1]);
            Dialog.show();
            ty = Dialog.getChoice();
            setTool("multipoint");
            waitForUser("Click each " + ty, "CLICK once on each " + ty + " in this field, then OK.\n(No clicks = none of this type.)");
            cNear = 0; cFar = 0; cnt = 0;
            if (selectionType() == 10) {
                getSelectionCoordinates(xs, ys); cnt = xs.length;
                for (k = 0; k < xs.length; k++) {
                    d = distEdgePx(xs[k], ys[k]);
                    nf = "far"; if (!isNaN(d) && d <= NEAR_FAR_UM) nf = "near";
                    if (nf == "near") cNear = cNear + 1; else cFar = cFar + 1;
                    writeNF(ty, k + 1, d, nf);
                }
            }
            run("Select None");
            showMessage(ty + ": " + cnt + " total  (" + cNear + " near, " + cFar + " far).");
            more = getBoolean("Mark another organelle type in this image?");
        }
    }

    // --- nuclear pores along the same traced envelope ---
    if (getBoolean("Count nuclear pores along this envelope?  (length " + d2s(envLen,2) + " um)")) {
        setTool("multipoint");
        waitForUser("Nuclear pores", "CLICK each nuclear pore along the traced envelope, then OK.");
        count = 0;
        if (selectionType() == 10) { getSelectionCoordinates(pxs, pys); count = pxs.length; }
        dens = NaN; if (!isNaN(envLen) && envLen > 0) dens = count / envLen;
        writePore(envLen, count, dens);
        run("Select None");
        if (!isNaN(dens)) showMessage(count + " pores.  Density = " + d2s(dens,3) + " pores/um of envelope.");
        else showMessage(count + " pores recorded.");
    }
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

// ---- 32-bit distance map from the traced membrane LINE (distance to membrane) --
function buildDistMap() {
    st = selectionType();
    if (!(st == 5 || st == 6 || st == 7)) return;
    orig = getTitle(); w = getWidth(); h = getHeight();
    getPixelSize(u, pw, ph); REF_PW = pw;
    roiManager("reset"); roiManager("add");     // store the membrane line
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    if (isOpen("__nucBinNF")) { selectWindow("__nucBinNF"); close(); }
    newImage("__nucBinNF", "8-bit black", w, h, 1);
    selectWindow("__nucBinNF");
    roiManager("select", 0);                     // restore the line here
    setForegroundColor(255, 255, 255);
    run("Line Width...", "line=1");
    run("Draw", "slice");                        // draw the membrane as 255
    run("Select None");
    setOption("BlackBackground", true);
    run("Options...", "iterations=1 count=1 black edm=32-bit");
    run("Invert");                               // membrane -> 0, rest -> 255
    run("Distance Map");                         // each pixel = distance (px) to membrane
    mapT = getTitle();
    if (mapT != "__nucBinNF" && isOpen("__nucBinNF")) { selectWindow("__nucBinNF"); close(); }
    selectWindow(mapT); rename(DMAP);
    selectWindow(orig); roiManager("reset");
}
// distance (um) from a pixel coordinate to the traced membrane
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
