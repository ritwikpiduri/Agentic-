// ============================================================================
//  LabPack_Nuclear.ijm  --  automated NUCLEAR TEM measurements for any lab/scope.
//
//  Same measurements as a manual Fiji session, but driven by a menu so you just
//  trace and click. Portable: you calibrate from YOUR scale bar (or pixel size),
//  not a fixed camera table, so it works on any microscope / cell model.
//
//  Per image, a menu lets you:
//    Set scale          -> calibrate THIS image from its scale bar (if needed)
//    Nucleus            -> trace the nucleus -> area, perimeter, circularity,
//                          aspect ratio, roundness, solidity, and heterochromatin
//                          % using a FIXED grey cutoff (set once; see note below).
//                          Tracing a 2nd/3rd nucleus in one image = multinucleated.
//    Micronuclei        -> trace each small separate body -> area, circularity
//    Perinuclear space  -> draw a line across the envelope gap -> nm
//    NEXT / SKIP / QUIT
//
//  OUTPUT: LabPack_Nuclear.csv (one row per measurement) + a QC picture of every
//  traced nucleus. Saved after EVERY entry, so Esc / QUIT never loses data.
//
//  CHROMATIN NOTE (important): heterochromatin % is the fraction of the nucleus
//  DARKER than a FIXED grey cutoff you set once at the start. Do NOT use a
//  per-nucleus/auto threshold -- it self-normalises and makes every sample look
//  equal. Keep imaging (staining + exposure) identical across all conditions,
//  and use the same cutoff for the whole dataset. MeanGrey is saved as a
//  cutoff-free backup (lower mean = darker = more heterochromatin).
// ============================================================================

var TBL     = "LabPack_Nuclear";
var _curImg = "";
var _curDir = "";
var _outDir = "";
var _qcDir  = "";

var CAL_MODE = "scalebar";   // "scalebar" | "fixed" | "embedded"
var FIXED_PX = 0;            // um/pixel, when CAL_MODE = "fixed"
var CUTOFF   = 100;          // grey 0-255: pixels <= this = heterochromatin
var _calOK   = false;        // is the current image calibrated?
var _calPx   = 0;            // current image um/pixel (for display)

inDir   = getDirectory("Choose the FOLDER of images for ONE sample/condition");
_outDir = getDirectory("Choose an OUTPUT folder for LabPack_Nuclear.csv");
_qcDir  = _outDir + "QC_nucleus/";
File.makeDirectory(_qcDir);

// ---- resume or start fresh ----
csvPath = _outDir + "LabPack_Nuclear.csv";
if (File.exists(csvPath)) {
    if (getBoolean("Found an existing LabPack_Nuclear.csv here.\n \nYES = RESUME (keep it, add to it)\nNO = start NEW (old one renamed to _prev)")) {
        open(csvPath); TBL = "LabPack_Nuclear.csv";
    } else {
        File.rename(csvPath, _outDir + "LabPack_Nuclear_prev.csv");
        if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
        Table.create(TBL);
    }
} else {
    if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
    Table.create(TBL);
}

// ---- calibration + chromatin settings (once) ----
Dialog.create("Settings");
Dialog.addChoice("Calibration:", newArray(
    "Set scale from the scale bar (recommended)",
    "Fixed pixel size for all images (same magnification)",
    "Images are already calibrated (embedded scale)"), "Set scale from the scale bar (recommended)");
Dialog.addNumber("If fixed: pixel size (um per pixel):", 0.001);
Dialog.addMessage("Heterochromatin cutoff: pixels darker than this (0-255) count as\nheterochromatin. Pick ONE value and use it for every condition.");
Dialog.addNumber("Heterochromatin grey cutoff (0-255):", CUTOFF);
Dialog.show();
cm = Dialog.getChoice();
FIXED_PX = Dialog.getNumber();
CUTOFF   = Dialog.getNumber();
if      (startsWith(cm, "Fixed")) CAL_MODE = "fixed";
else if (startsWith(cm, "Images")) CAL_MODE = "embedded";
else CAL_MODE = "scalebar";

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
        run("Enhance Contrast", "saturated=0.35");   // view only
        startCal();                                   // set/detect calibration

        imgDone = false;
        while (!imgDone) {
            calStr = "NOT calibrated - choose 'Set scale' first";
            if (_calOK) calStr = "calibrated: " + d2s(_calPx,5) + " um/px";
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg + "\n[" + calStr + "]");
            Dialog.addChoice("Action:", newArray(
                "Set scale  (draw on the scale bar)",
                "Nucleus  (trace - shape + chromatin + area)",
                "Micronuclei  (trace each)",
                "Perinuclear space  (line across the gap)",
                "NEXT image",
                "SKIP this image",
                "QUIT and save"), "Nucleus  (trace - shape + chromatin + area)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Set scale"))         actSetScale();
            else if (startsWith(a, "Nucleus"))           actNucleus();
            else if (startsWith(a, "Micronuclei"))       actMicronuclei();
            else if (startsWith(a, "Perinuclear"))       actMembrane();
            else if (startsWith(a, "NEXT"))              imgDone = true;
            else if (startsWith(a, "SKIP"))              imgDone = true;
            else if (startsWith(a, "QUIT"))              { imgDone = true; quit = true; }
        }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(_outDir + "LabPack_Nuclear.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + _outDir + "LabPack_Nuclear.csv");

// ============================================================================
// ---- calibration ----
function startCal() {
    _calOK = false; _calPx = 0;
    if (CAL_MODE == "fixed") {
        if (FIXED_PX > 0) { setVoxelSize(FIXED_PX, FIXED_PX, 1, "micron"); _calOK = true; _calPx = FIXED_PX; }
    } else if (CAL_MODE == "embedded") {
        getPixelSize(u, pw, ph);
        if (pw > 0 && pw != 1) { _calOK = true; _calPx = pw; }
    }
    // "scalebar" mode: stays uncalibrated until the user sets it
}
function actSetScale() {
    setTool("line");
    run("Select None");
    waitForUser("Scale bar", "Draw a straight line EXACTLY along the scale bar\n(match its length), then OK.");
    if (selectionType() != 5) { showMessage("Draw a straight line along the scale bar."); return; }
    getPixelSize(u, pw, ph);
    lenNow = getValue("Length");
    lenPx = lenNow; if (pw > 0) lenPx = lenNow / pw;   // length back in pixels
    if (lenPx <= 0) { showMessage("Could not read the line length."); return; }
    Dialog.create("Scale bar value");
    Dialog.addNumber("Scale bar length (the number printed on it):", 1);
    Dialog.addChoice("Its unit:", newArray("um (micron)", "nm"), "um (micron)");
    Dialog.show();
    known = Dialog.getNumber();
    unit  = Dialog.getChoice();
    knownUm = known; if (startsWith(unit, "nm")) knownUm = known / 1000.0;
    pxUm = knownUm / lenPx;
    setVoxelSize(pxUm, pxUm, 1, "micron");
    _calOK = true; _calPx = pxUm;
    run("Select None");
    showMessage("Scale set: " + d2s(pxUm,5) + " um/pixel.");
}

// ---- nucleus: shape + fixed-cutoff chromatin ----
function actNucleus() {
    if (!requireCal()) return;
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
    getHistogram(hv, hc, 256);
    dark = 0; tot = 0;
    for (k = 0; k < hv.length; k++) { tot = tot + hc[k]; if (hv[k] <= CUTOFF) dark = dark + hc[k]; }
    het = NaN; if (tot > 0) het = 100.0 * dark / tot;
    eu = NaN;  if (!isNaN(het)) eu = 100 - het;
    writeNucleus(n, getResult("Area", m), getResult("Perim.", m), getResult("Circ.", m),
        getResult("AR", m), getResult("Round", m), getResult("Solidity", m), het, eu, nMean);
    // QC picture
    selectWindow(_curImg); run("Remove Overlay"); roiManager("select", 0);
    Roi.setStrokeColor("yellow"); Roi.setStrokeWidth(4); run("Add Selection...");
    run("Flatten"); saveAs("PNG", _qcDir + noExt(_curImg) + "_n" + n + "_QC.png"); close();
    selectWindow(_curImg); run("Select None"); roiManager("reset");
    msg = "Nucleus " + n + " recorded.";
    if (n > 1) msg = msg + "\n>> Multinucleated (" + n + " nuclei in this image).";
    showMessage(msg);
}

function actMicronuclei() {
    if (!requireCal()) return;
    setMeas(); setTool("freehand");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser("Micronucleus #" + i, "Trace a small body SEPARATE from the nucleus, then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); m = nResults - 1;
            writeMicro(i, getResult("Area", m), getResult("Circ.", m));
        }
        more = getBoolean("Mark another micronucleus?");
    }
    run("Select None"); showMessage(i + " micronucleus/nuclei recorded.");
}

function actMembrane() {
    if (!requireCal()) return;
    setTool("line");
    n = getNumber("How many measurements on this image?", 5);
    i = 0;
    while (i < n) {
        waitForUser("Perinuclear space " + (i+1) + " of " + n,
            "Draw a short line ACROSS the gap between the inner and outer\n" +
            "nuclear membranes (perpendicular to the envelope), then OK.");
        if (selectionType() == 5) {
            lenUm = getValue("Length");
            writeMembrane(i + 1, lenUm * 1000.0);   // um -> nm
            i = i + 1;
        }
    }
    run("Select None"); showMessage(n + " perinuclear-space measurements recorded.");
}

function requireCal() {
    if (_calOK) return true;
    showMessage("Set the scale first (choose 'Set scale' and draw on the scale bar).");
    return false;
}

// ---- writers (one row each) ----
function writeNucleus(idx, area, perim, circ, ar, rnd, sol, het, eu, mg) {
    r = Table.size(TBL); smp = folderOf();
    base(r, smp, "Nucleus", idx);
    Table.set("Area_um2", r, area, TBL);
    Table.set("Perimeter_um", r, perim, TBL);
    Table.set("Circularity", r, circ, TBL);
    Table.set("AspectRatio", r, ar, TBL);
    Table.set("Roundness", r, rnd, TBL);
    Table.set("Solidity", r, sol, TBL);
    Table.set("Heterochrom_pct", r, het, TBL);
    Table.set("Euchrom_pct", r, eu, TBL);
    Table.set("MeanGrey", r, mg, TBL);
    Table.set("Cutoff", r, CUTOFF, TBL);
    flush();
}
function writeMicro(idx, area, circ) {
    r = Table.size(TBL); smp = folderOf();
    base(r, smp, "Micronucleus", idx);
    Table.set("Area_um2", r, area, TBL);
    Table.set("Circularity", r, circ, TBL);
    flush();
}
function writeMembrane(idx, nm) {
    r = Table.size(TBL); smp = folderOf();
    base(r, smp, "Perinuclear", idx);
    Table.set("PerinuclearSpace_nm", r, nm, TBL);
    flush();
}
function base(r, smp, rec, idx) {
    Table.set("Sample", r, smp, TBL);
    Table.set("Image", r, _curImg, TBL);
    Table.set("RecordType", r, rec, TBL);
    Table.set("Index", r, idx, TBL);
}
function flush() { Table.update(TBL); selectWindow(TBL); Table.save(_outDir + "LabPack_Nuclear.csv"); }

// ---- helpers ----
function listTiffs(dir) {
    out = newArray(0); l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        if (endsWith(l[i], "/")) out = Array.concat(out, listTiffs(dir + l[i]));
        else { nm = toLowerCase(l[i]); if (endsWith(nm, ".tif") || endsWith(nm, ".tiff")) out = Array.concat(out, dir + l[i]); }
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
function folderOf() { p = split(_curDir, "/\\"); k = p.length; if (k >= 1) return p[k-1]; return _curDir; }
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
