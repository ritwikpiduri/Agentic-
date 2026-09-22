// ============================================================================
//  TEM_OrganelleCount.ijm  --  count organelles per field (abundance / density).
//  A separate, simple macro just for "are there MORE or FEWER organelles?".
//  Works on ANY image (nucleus not needed), including pure organelle zoom-ins.
//
//  Per field you:
//    1) trace the region you are counting within (the cell/cytoplasm area),
//    2) for each organelle type, CLICK each one (multipoint) -> a count.
//  It records, per type:  Count  and  FieldArea_um2  and  Density_per_um2
//  ( Density = Count / area ), so the comparison across samples is FAIR
//  (raw counts are not comparable; density per um^2 is).
//
//  OUTPUT: TEM_OrganelleCount.csv, saved after every entry (Esc/QUIT safe).
//  Use ONE output folder per sample so files don't overwrite.
// ============================================================================

// ---- calibration (your pixel sizes, um/pixel) ----
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;

var TYPES = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody", "Lysosome");
var TBL = "TEM_OrganelleCount";
var _curImg = "";
var _curDir = "";
var _outDir = "";

inDir   = getDirectory("Choose the FOLDER of images for ONE sample");
_outDir = getDirectory("Choose an OUTPUT folder for TEM_OrganelleCount.csv");

csvPath = _outDir + "TEM_OrganelleCount.csv";
if (File.exists(csvPath)) {
    if (getBoolean("Found an existing TEM_OrganelleCount.csv here.\n \nYES = RESUME (keep it, add to it)\nNO = start NEW (old one renamed to _prev)")) {
        open(csvPath); TBL = "TEM_OrganelleCount.csv";
    } else {
        File.rename(csvPath, _outDir + "TEM_OrganelleCount_prev.csv");
        if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
        Table.create(TBL);
    }
} else {
    if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
    Table.create(TBL);
}

files = listTiffs(inDir);
if (files.length == 0) { showMessage("No .tif/.tiff images found."); exit; }

startAt = getNumber("Start from image number (1 = beginning).", 1);
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
        autoCalibrate();
        run("Enhance Contrast", "saturated=0.35");

        Dialog.create("Image " + (f+1) + " / " + files.length);
        Dialog.addMessage(_curImg);
        Dialog.addChoice("Action:", newArray("Count this field", "SKIP this image", "QUIT and save"), "Count this field");
        Dialog.show();
        a = Dialog.getChoice();
        if (startsWith(a, "Count")) countField();
        else if (startsWith(a, "QUIT")) quit = true;

        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(_outDir + "TEM_OrganelleCount.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + _outDir + "TEM_OrganelleCount.csv");

// ============================================================================
function countField() {
    setMeasAreaOnly(); setTool("freehand");
    run("Select None");
    waitForUser("Counting area",
        "Trace the region you are counting within (the cell/cytoplasm area in\n" +
        "this field), then OK.  (No selection = the whole image is used.)");
    getPixelSize(u, pw, ph);
    if (selectionType() >= 0) { run("Measure"); m = nResults - 1; fieldA = getResult("Area", m); }
    else { fieldA = getWidth() * getHeight() * pw * ph; }
    run("Select None");

    more = true;
    while (more) {
        Dialog.create("Count which organelle?");
        Dialog.addChoice("Type:", TYPES, TYPES[1]);
        Dialog.show();
        ct = Dialog.getChoice();
        setTool("multipoint");
        waitForUser("Count " + ct, "CLICK each " + ct + " in the counting area, then OK.\n(No clicks = zero of this type.)");
        cnt = 0;
        if (selectionType() == 10) { getSelectionCoordinates(xs, ys); cnt = xs.length; }
        writeCount(ct, cnt, fieldA);
        run("Select None");
        more = getBoolean("Count another organelle type in this field?");
    }
    showMessage("Counts recorded (field area " + d2s(fieldA,2) + " um2).");
}

function writeCount(type, cnt, fieldA) {
    r = Table.size(TBL);
    dens = NaN; if (fieldA > 0) dens = cnt / fieldA;
    Table.set("Image", r, _curImg, TBL);
    Table.set("Sample", r, folderOf(), TBL);
    Table.set("Type", r, type, TBL);
    Table.set("Count", r, cnt, TBL);
    Table.set("FieldArea_um2", r, fieldA, TBL);
    Table.set("Density_per_um2", r, dens, TBL);
    Table.update(TBL);
    selectWindow(TBL); Table.save(_outDir + "TEM_OrganelleCount.csv");
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
function setMeasAreaOnly() { run("Set Measurements...", "area redirect=None decimal=4"); }
function folderOf() {
    p = split(_curDir, "/\\"); n = p.length;
    if (n >= 1) return p[n-1];
    return _curDir;
}
