// ============================================================================
//  TEM_Organelle.ijm  --  SESSION 2: organelles + nuclear pores.
//  Run this on your HIGH-MAG images (do this after the nuclear session).
//
//  Per image, a menu lets you:
//    Nucleus reference -> TRACE the nuclear envelope shown in the frame
//                         (freehand, following the membrane / around the nucleus).
//                         This is the reference for organelle near/far + pores.
//    Organelles        -> pick a type, then drag a quick OVAL over each one ->
//                         area, circularity, distance to nucleus, near/far.
//    Nuclear pores     -> click each pore -> count + pores per um of envelope.
//    NEXT / SKIP / QUIT
//
//  OUTPUT: one file TEM_Organelle.csv. Saved after EVERY measurement, so
//  Esc / QUIT never loses data. Same columns as TEM_Nuclear.csv, so the two
//  can be merged for the final analysis.
// ============================================================================

// ---- calibration (your pixel sizes, um/pixel) ----
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;

var NEAR_FAR_UM = 1.0;     // organelle <= this from the nucleus edge = "near"
var TBL   = "TEM_Organelle";
var DMAP  = "__nucDist";
var _curImg = "";
var _curDir = "";
var _outDir = "";

var HAVE_REF = false;
var REF_PERIM = 0;         // traced nucleus perimeter (um), for pore density

inDir   = getDirectory("Choose the FOLDER of images (high-mag / organelles)");
_outDir = getDirectory("Choose an OUTPUT folder for TEM_Organelle.csv");

csvPath = _outDir + "TEM_Organelle.csv";
if (File.exists(csvPath)) {
    if (getBoolean("Found an existing TEM_Organelle.csv here.\n \nYES = RESUME (keep it and add to it)\nNO = start a NEW file (old one renamed to _prev)")) {
        open(csvPath);
        TBL = "TEM_Organelle.csv";
    } else {
        File.rename(csvPath, _outDir + "TEM_Organelle_prev.csv");
        if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
        Table.create(TBL);
    }
} else {
    if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
    Table.create(TBL);
}

Dialog.create("Settings");
Dialog.addNumber("Organelle near/far cutoff (um):", NEAR_FAR_UM);
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
        HAVE_REF = false; REF_PERIM = 0;
        autoCalibrate();
        run("Enhance Contrast", "saturated=0.35");   // view only

        imgDone = false;
        while (!imgDone) {
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg);
            Dialog.addChoice("Action:", newArray(
                "Nucleus reference  (trace the nucleus - do first)",
                "Organelles  (oval each - size + near/far)",
                "Nuclear pores  (click each)",
                "NEXT image",
                "SKIP this image",
                "QUIT and save"), "Nucleus reference  (trace the nucleus - do first)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Nucleus reference")) actNucleusRef();
            else if (startsWith(a, "Organelles"))        actOrganelles();
            else if (startsWith(a, "Nuclear pores"))     actPores();
            else if (startsWith(a, "NEXT"))              imgDone = true;
            else if (startsWith(a, "SKIP"))              imgDone = true;
            else if (startsWith(a, "QUIT"))              { imgDone = true; quit = true; }
        }
        if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(_outDir + "TEM_Organelle.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + _outDir + "TEM_Organelle.csv");

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
    selectWindow(TBL); Table.save(_outDir + "TEM_Organelle.csv");
}

// ---- nucleus reference: trace the nucleus -> distance map + perimeter ----
function actNucleusRef() {
    setMeas(); setTool("freehand");
    run("Select None");
    waitForUser("Nucleus reference (trace)",
        "Trace AROUND the nucleus (freehand), following the envelope,\n" +
        "then OK. This is the reference for organelle near/far and pore density.");
    if (selectionType() < 0) { showMessage("Nothing traced - skipped."); return; }
    run("Measure"); m = nResults - 1;
    REF_PERIM = getResult("Perim.", m);
    buildDistMap();
    run("Select None");
    HAVE_REF = true;
    showMessage("Nucleus reference set (perimeter " + d2s(REF_PERIM,2) + " um).\nNow mark organelles / pores.");
}

// ---- organelles: quick oval each ----
function actOrganelles() {
    types = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody", "Lysosome");
    Dialog.create("Organelle type");
    Dialog.addChoice("Marking which organelle?", types, types[0]);
    Dialog.show();
    otype = Dialog.getChoice();
    if (!HAVE_REF) showMessage("Tip: trace the Nucleus reference first so near/far is recorded.");
    setMeas(); setTool("oval");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser(otype + " #" + i, "Drag a quick OVAL over this " + otype + ", then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); m = nResults - 1;
            cx = getResult("X", m); cy = getResult("Y", m);
            dist = NaN; nf = "";
            if (HAVE_REF) { dist = distEdge(cx, cy); if (dist <= NEAR_FAR_UM) nf = "near"; else nf = "far"; }
            writeRow("Organelle", otype, i, getResult("Area", m), NaN, getResult("Circ.", m),
                NaN, NaN, NaN, NaN, NaN, dist, nf, NaN, NaN, NaN);
        }
        more = getBoolean("Mark another " + otype + "?");
    }
    run("Select None"); showMessage(i + " " + otype + "(s) recorded.");
}

// ---- nuclear pores ----
function actPores() {
    setTool("multipoint");
    waitForUser("Nuclear pores", "Multi-point tool: CLICK each nuclear pore, then OK.");
    if (selectionType() != 10) { showMessage("Use the multi-point tool and click the pores."); return; }
    getSelectionCoordinates(xs, ys); count = xs.length;
    dens = NaN; if (HAVE_REF && REF_PERIM > 0) dens = count / REF_PERIM;
    writeRow("Pore", "", 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "", NaN, count, dens);
    run("Select None");
    pmsg = count + " pores.";
    if (HAVE_REF && REF_PERIM > 0) pmsg = pmsg + "\nDensity = " + d2s(dens,3) + " pores/um of envelope.";
    else pmsg = pmsg + "\n(Trace the Nucleus reference first for density.)";
    showMessage(pmsg);
}

// ---- distance map from the traced nucleus (distance to the nuclear edge) ----
function buildDistMap() {
    if (selectionType() < 0) return;
    orig = getTitle(); w = getWidth(); h = getHeight();
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    newImage(DMAP, "8-bit black", w, h, 1);
    selectWindow(orig); roiManager("reset"); roiManager("add");
    selectWindow(DMAP); roiManager("select", 0); setColor(255); fill(); run("Select None");
    setOption("BlackBackground", true); run("Invert"); run("Distance Map");
    selectWindow(orig); roiManager("reset");
}
function distEdge(cx, cy) {
    if (!isOpen(DMAP)) return NaN;
    cur = getTitle(); selectWindow(DMAP); getPixelSize(u, pw, ph);
    px = round(cx / pw); py = round(cy / ph);
    if (px < 0) px = 0; if (py < 0) py = 0;
    if (px >= getWidth()) px = getWidth() - 1;
    if (py >= getHeight()) py = getHeight() - 1;
    d = getPixel(px, py) * pw; selectWindow(cur); return d;
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
function folderOf() {
    p = split(_curDir, "/\\"); n = p.length;
    if (n >= 1) return p[n-1];
    return _curDir;
}
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
