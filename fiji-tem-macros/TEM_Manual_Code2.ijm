// =============================================================================
//  TEM_Manual_Code2.ijm  = CODE 2  (manual marking, everything else automatic)
//  Opens each image in a folder for you, auto-calibrates from the filename
//  magnification, and lets you mark structures. Everything is measured, counted,
//  and saved into ONE file: TEM_Results.csv
//
//  Each row is one measurement, tagged by RecordType:
//    Nucleus      -> area, perimeter, circularity, AR, roundness, solidity,
//                    heterochromatin %, euchromatin %   (2+ per image = multinucleated)
//    Organelle    -> Type (ER/Mito/Golgi/Vacuole/LipidBody), area, circularity,
//                    distance to nucleus, near/far
//    Micronucleus -> area, circularity
//    Membrane     -> thickness (nm)
//    Pore         -> pore count + pores per um of nuclear envelope
//  Sample = the image's folder name (for grouping the 4 conditions).
// =============================================================================

var NEAR_FAR_UM = 1.0;
var DMAP = "__nucDistMap";
var HAVE_DMAP = false;
var _curImg = "";
var _curDir = "";
var TBL = "TEM_Results";

// --- calibration ---
// Exact pixel sizes (um/pixel) you gave for these magnifications:
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
// Any other magnification (e.g. 4500) is estimated from the camera constant
// (magnification x pixel-size(nm) ~ 14400):  um/px = 14.4 / magnification
var CAM_CONST = 14.4;

inDir  = getDirectory("Choose the FOLDER of images to analyse");
outDir = getDirectory("Choose an OUTPUT folder for TEM_Results.csv");
if (!isOpen(TBL)) Table.create(TBL);

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

        imgDone = false;
        while (!imgDone) {
            Dialog.create("Image " + (f+1) + " / " + files.length);
            Dialog.addMessage(_curImg);
            Dialog.addChoice("Action:", newArray(
                "Trace nucleus (shape + chromatin + multinucleation)",
                "Mark organelles (type, size, near/far)",
                "Mark micronuclei",
                "Membrane / perinuclear thickness",
                "Nuclear pores (count + density)",
                "NEXT image",
                "SKIP this image",
                "QUIT and save"), "Trace nucleus (shape + chromatin + multinucleation)");
            Dialog.show();
            a = Dialog.getChoice();
            if      (startsWith(a, "Trace nucleus"))    actNucleus();
            else if (startsWith(a, "Mark organelles"))  actOrganelles();
            else if (startsWith(a, "Mark micronuclei")) actMicronuclei();
            else if (startsWith(a, "Membrane"))         actMembrane();
            else if (startsWith(a, "Nuclear pores"))    actPores();
            else if (startsWith(a, "NEXT"))             imgDone = true;
            else if (startsWith(a, "SKIP"))             imgDone = true;
            else if (startsWith(a, "QUIT"))             { imgDone = true; quit = true; }
        }
        if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
        if (nImages > 0) { selectWindow(_curImg); close(); }
        done = done + 1;
    }
}
selectWindow(TBL);
Table.save(outDir + "TEM_Results.csv");
showMessage("Done", "Processed " + done + " image(s).\nSaved: " + outDir + "TEM_Results.csv");

// =============================================================================
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
}

function actNucleus() {
    setMeas(); setTool("freehand");
    waitForUser("Trace the NUCLEUS",
        "Outline the nucleus (freehand), then OK.\n(Repeat for each extra nucleus in this image.)");
    if (selectionType() < 0) { showMessage("No selection."); return; }
    run("Measure"); row = nResults - 1;
    n = countRec("Nucleus", _curImg) + 1;
    getStatistics(nArea, nMean, nMin, nMax, nStd);
    thr = round(nMean - 0.5 * nStd);
    // fast chromatin: histogram of the nucleus selection (instant, no pixel loop)
    getHistogram(hvals, hcounts, 256);
    dark = 0; tot = 0;
    for (k = 0; k < hvals.length; k++) { tot = tot + hcounts[k]; if (hvals[k] <= thr) dark = dark + hcounts[k]; }
    het = NaN; if (tot > 0) het = 100.0 * dark / tot;
    writeRow("Nucleus", "", n,
        getResult("Area", row), getResult("Perim.", row), getResult("Circ.", row),
        getResult("AR", row), getResult("Round", row), getResult("Solidity", row),
        het, 100 - het, NaN, "", NaN, NaN, NaN);
    buildDistMap(); run("Select None");
    showMessage("Nucleus " + n + " recorded." + (n > 1 ? "\n>> Multinucleated." : ""));
}

function actOrganelles() {
    types = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody");
    Dialog.create("Organelle type");
    Dialog.addChoice("Marking which organelle?", types, types[0]);
    Dialog.addNumber("Near/far cutoff (um):", NEAR_FAR_UM);
    Dialog.show();
    otype = Dialog.getChoice(); NEAR_FAR_UM = Dialog.getNumber();
    if (!HAVE_DMAP) showMessage("Tip: trace the nucleus first for distances.");
    setMeas(); setTool("freehand");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser(otype + " #" + i, "Trace this " + otype + " (freehand), then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); row = nResults - 1;
            cx = getResult("X", row); cy = getResult("Y", row);
            dist = NaN; nf = "";
            if (HAVE_DMAP) { dist = distEdge(cx, cy); if (dist <= NEAR_FAR_UM) nf = "near"; else nf = "far"; }
            writeRow("Organelle", otype, i, getResult("Area", row), NaN, getResult("Circ.", row),
                NaN, NaN, NaN, NaN, NaN, dist, nf, NaN, NaN, NaN);
        }
        more = getBoolean("Mark another " + otype + "?");
    }
    run("Select None"); showMessage(i + " " + otype + "(s) recorded.");
}

function actMicronuclei() {
    setMeas(); setTool("freehand");
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser("Micronucleus #" + i, "Trace a small body separate from the nucleus, then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure"); row = nResults - 1;
            writeRow("Micronucleus", "", i, getResult("Area", row), NaN, getResult("Circ.", row),
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
    showMessage(count + " pores." + (perim > 0 ? "\nDensity = " + d2s(dens,3) + " pores/um." : "\n(Trace nucleus first for density.)"));
}

// ---- distance map from the nucleus ROI ----
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
    mag = parseMag(getTitle());
    if (mag <= 0) return false;
    px = pxForMag(mag);
    setVoxelSize(px, px, 1, "micron");
    return true;
}
function parseMag(name) {
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

// ---- file walking + helpers ----
function listTiffs(dir) {
    out = newArray(0); l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        if (endsWith(l[i], "/")) out = Array.concat(out, listTiffs(dir + l[i]));
        else { n = toLowerCase(l[i]); if (endsWith(n, ".tif") || endsWith(n, ".tiff")) out = Array.concat(out, dir + l[i]); }
    }
    return out;
}
function setMeas() { run("Set Measurements...", "area mean centroid perimeter shape redirect=None decimal=3"); }
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
