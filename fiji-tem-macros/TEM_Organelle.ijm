// ============================================================================
//  TEM_Organelle.ijm  --  SESSION 2: organelles + nuclear pores.
//  Run this on your HIGH-MAG images (do this after the nuclear session).
//
//  Per image, a menu lets you:
//    Nucleus reference -> draw ONE line along the nuclear envelope in the frame.
//                         This is the reference for organelle near/far.
//    Organelles        -> pick a type, then drag a quick OVAL over each one ->
//                         area, circularity, distance to nucleus, near/far.
//    Nuclear pores     -> click each pore -> count (+ pores per um of the
//                         reference line, if you drew one).
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
var _curImg = "";
var _curDir = "";
var _outDir = "";

var HAVE_REF = false;
var LX1 = 0; var LY1 = 0; var LX2 = 0; var LY2 = 0;   // envelope line (microns)
var LINE_LEN = 0;          // length of the reference line (microns)

inDir   = getDirectory("Choose the FOLDER of images (high-mag / organelles)");
_outDir = getDirectory("Choose an OUTPUT folder for TEM_Organelle.csv");

if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
Table.create(TBL);

Dialog.create("Settings");
Dialog.addNumber("Organelle near/far cutoff (um):", NEAR_FAR_UM);
Dialog.show();
NEAR_FAR_UM = Dialog.getNumber();

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
        HAVE_REF = false; LINE_LEN = 0;
        autoCalibrate();
        run("Enhance Contrast", "saturated=0.35");   // view only

        imgDone = false;
        while (!imgDone) {
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg);
            Dialog.addChoice("Action:", newArray(
                "Nucleus reference  (draw line on envelope - do first)",
                "Organelles  (oval each - size + near/far)",
                "Nuclear pores  (click each)",
                "NEXT image",
                "SKIP this image",
                "QUIT and save"), "Nucleus reference  (draw line on envelope - do first)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Nucleus reference")) actNucleusRef();
            else if (startsWith(a, "Organelles"))        actOrganelles();
            else if (startsWith(a, "Nuclear pores"))     actPores();
            else if (startsWith(a, "NEXT"))              imgDone = true;
            else if (startsWith(a, "SKIP"))              imgDone = true;
            else if (startsWith(a, "QUIT"))              { imgDone = true; quit = true; }
        }
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

// ---- nucleus reference line ----
function actNucleusRef() {
    setTool("line");
    run("Select None");
    waitForUser("Nucleus reference (distance)",
        "Draw ONE line ALONG the nuclear envelope shown in this image,\n" +
        "then OK. This is the reference for organelle near/far\n(and for pore density).");
    if (selectionType() != 5) { showMessage("Draw a straight LINE along the envelope."); return; }
    getLine(x1, y1, x2, y2, lw);
    getPixelSize(u, pw, ph);
    LX1 = x1 * pw; LY1 = y1 * ph; LX2 = x2 * pw; LY2 = y2 * ph;
    LINE_LEN = sqrt((LX2-LX1)*(LX2-LX1) + (LY2-LY1)*(LY2-LY1));
    HAVE_REF = true;
    run("Select None");
    showMessage("Nucleus reference set (length " + d2s(LINE_LEN,3) + " um).\nNow mark organelles.");
}

// ---- organelles: quick oval each ----
function actOrganelles() {
    types = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody");
    Dialog.create("Organelle type");
    Dialog.addChoice("Marking which organelle?", types, types[0]);
    Dialog.show();
    otype = Dialog.getChoice();
    if (!HAVE_REF) showMessage("Tip: draw the Nucleus reference line first so near/far is recorded.");
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
            if (HAVE_REF) { dist = pointSegDist(cx, cy, LX1, LY1, LX2, LY2); if (dist <= NEAR_FAR_UM) nf = "near"; else nf = "far"; }
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
    dens = NaN; if (HAVE_REF && LINE_LEN > 0) dens = count / LINE_LEN;
    writeRow("Pore", "", 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "", NaN, count, dens);
    run("Select None");
    pmsg = count + " pores.";
    if (HAVE_REF && LINE_LEN > 0) pmsg = pmsg + "\nDensity = " + d2s(dens,3) + " pores/um (along the reference line).";
    else pmsg = pmsg + "\n(Draw the Nucleus reference line first for density.)";
    showMessage(pmsg);
}

// ---- distance: point to line segment (microns) ----
function pointSegDist(px, py, x1, y1, x2, y2) {
    dx = x2 - x1; dy = y2 - y1;
    if (dx == 0 && dy == 0) return sqrt((px-x1)*(px-x1) + (py-y1)*(py-y1));
    t = ((px-x1)*dx + (py-y1)*dy) / (dx*dx + dy*dy);
    if (t < 0) t = 0; if (t > 1) t = 1;
    qx = x1 + t*dx; qy = y1 + t*dy;
    return sqrt((px-qx)*(px-qx) + (py-qy)*(py-qy));
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
