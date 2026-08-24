// ============================================================================
//  TEM_AutoAll.ijm  --  FULLY AUTOMATIC. No tracing, no clicking.
//  Run it -> pick the image folder -> pick an output folder -> walk away.
//
//  For EVERY .tif/.tiff (subfolders included) it auto-measures, per image:
//    NUCLEUS      circularity, aspect ratio, roundness, solidity, perimeter,
//                 area, heterochromatin %, euchromatin %
//    MULTINUCLEATION  how many nucleus-sized bodies were found (>=2 = flagged)
//    MICRONUCLEI  count of small round bodies separate from the main nucleus
//    DENSE BODIES count of electron-dense cytoplasmic structures (an organelle
//                 screen: mitochondria / lipid droplets / dense vacuoles), with
//                 how many are NEAR vs FAR from the nucleus (distance map)
//
//  Everything lands in ONE file:  TEM_AutoResults.csv  (one row per image).
//  A QC picture per image is saved so you can eyeball the auto-detections:
//  yellow = nucleus, cyan = micronucleus candidates, magenta = dense bodies.
//
//  Calibration is automatic from the magnification in the filename (same table
//  you already use). Images where the whole nucleus is not visible (very high
//  mag crops) are simply marked "no_nucleus_found" and skipped -- that is
//  normal and expected.
// ============================================================================

// ---- calibration: magnification -> microns per pixel (your values) ----
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;   // for any other mag: um/px = 14.4 / magnification

// ---- tuning (µm² if calibrated, else px²). Sensible TEM defaults ----
var HET_FACTOR   = 0.5;    // heterochromatin = pixels darker than mean-0.5*SD
var NUC_MIN      = 3.0;    // a body this big or bigger counts as a nucleus (µm²)
var MICRO_MIN    = 0.15;   // micronucleus size window
var MICRO_MAX    = 4.0;
var MICRO_CIRC   = 0.55;   // micronuclei are fairly round
var DENSE_MIN    = 0.02;   // dense-body (organelle) size window
var DENSE_MAX    = 3.0;
var NEAR_FAR_UM  = 1.0;    // dense body <= this from nucleus edge = "near"
var SAVE_QC      = true;

var TBL = "TEM_AutoResults";
var nTotal = 0;
var nDone  = 0;

setBatchMode(true);
inDir  = getDirectory("Choose the FOLDER with your TEM images");
outDir = getDirectory("Choose a folder to SAVE the results");
qcDir  = outDir + "QC_overlays/";
if (SAVE_QC) File.makeDirectory(qcDir);

if (isOpen(TBL)) { selectWindow(TBL); run("Close"); }
Table.create(TBL);
run("Set Measurements...", "area mean standard min centroid perimeter shape redirect=None decimal=4");

nTotal = countTiffs(inDir);
print("\\Clear");
print("Found " + nTotal + " images. Working through them automatically...");
walk(inDir);

selectWindow(TBL);
Table.save(outDir + "TEM_AutoResults.csv");
setBatchMode(false);
showMessage("Done!",
    "Auto-measured " + nDone + " of " + nTotal + " images.\n" +
    "Saved: " + outDir + "TEM_AutoResults.csv\n" +
    "QC pictures: " + qcDir + "\n\n" +
    "Yellow = nucleus, Cyan = micronucleus, Magenta = dense body.");

// ============================================================================
function walk(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) walk(dir + list[i]);
        else if (isTiff(list[i])) processOne(dir + list[i], dir, list[i]);
    }
}

function processOne(path, dir, name) {
    print("[" + (nDone + 1) + "/" + nTotal + "] " + name);
    showProgress(nDone, nTotal);
    open(path);
    if (bitDepth() == 24 || bitDepth() == 16) run("8-bit");

    // --- calibrate from filename magnification ---
    scaleSource = "pixels"; px = 1;
    mag = readMag(name);
    p = pxForMag(mag);
    if (p > 0) { px = p; setVoxelSize(px, px, 1, "micron"); scaleSource = "mag=" + mag; }
    calibrated = (scaleSource != "pixels");

    orig = getTitle(); w = getWidth(); h = getHeight();

    // --- find nucleus/nuclei (largest smooth dark regions) ---
    run("Duplicate...", "title=work");
    run("Gaussian Blur...", "sigma=3");
    setAutoThreshold("Otsu dark");
    run("Convert to Mask");
    run("Options...", "iterations=3 count=1 black");
    run("Fill Holes");
    roiManager("reset");
    run("Analyze Particles...", "size=" + NUC_MIN + "-Infinity circularity=0.1-1.00 show=Nothing add");

    frameA = w * h; if (calibrated) frameA = w * h * px * px;

    // rank nucleus-like ROIs by area (exclude near-frame-filling blobs)
    nRoi = roiManager("count");
    nucIdx = newArray(0); nucArea = newArray(0);
    for (r = 0; r < nRoi; r++) {
        selectWindow(orig); roiManager("select", r); getStatistics(a);
        if (a < 0.95 * frameA && a >= NUC_MIN * areaUnit(calibrated, px)) {
            nucIdx  = Array.concat(nucIdx, r);
            nucArea = Array.concat(nucArea, a);
        }
    }

    row = Table.size(TBL);
    fld = folderOf(dir);
    Table.set("Image", row, name, TBL);
    Table.set("Folder", row, fld, TBL);
    Table.set("ScaleSource", row, scaleSource, TBL);

    if (nucIdx.length == 0) {
        Table.set("Status", row, "no_nucleus_found", TBL);
        Table.update(TBL);
        closeAll(); return;
    }

    // biggest = the main nucleus
    best = nucIdx[0]; bestA = nucArea[0];
    for (i = 1; i < nucIdx.length; i++) if (nucArea[i] > bestA) { bestA = nucArea[i]; best = nucIdx[i]; }

    // count nucleus-sized bodies for multinucleation (>= half the main one)
    nucCount = 0;
    for (i = 0; i < nucArea.length; i++) if (nucArea[i] >= 0.5 * bestA) nucCount = nucCount + 1;

    // --- measure the main nucleus ---
    selectWindow(orig); roiManager("select", best); run("Measure"); m = nResults - 1;
    area = getResult("Area", m);
    getStatistics(dummyA, nucMean, nucMin, nucMax, nucSD);
    areaPx = area; areaUm = NaN;
    if (calibrated) { areaUm = area; areaPx = area / (px * px); }

    // --- chromatin: fast histogram inside the nucleus selection ---
    thr = round(nucMean - HET_FACTOR * nucSD);
    getHistogram(hv, hc, 256);
    dark = 0; tot = 0;
    for (k = 0; k < hv.length; k++) { tot = tot + hc[k]; if (hv[k] <= thr) dark = dark + hc[k]; }
    hetPct = NaN; if (tot > 0) hetPct = 100.0 * dark / tot;

    // --- micronucleus candidates: small round bodies, not the nucleus ---
    micro = countMicronuclei(orig, best, calibrated, px);

    // --- dense cytoplasmic bodies + near/far via distance map ---
    buildNucDistMap(orig, best, w, h);
    dCount = 0; dNear = 0; dFar = 0; dArea = 0;
    denseScan(orig, best, calibrated, px);   // fills global dg_* below
    dCount = dg_count; dNear = dg_near; dFar = dg_far; dArea = dg_area;

    // --- write the row ---
    Table.set("Status", row, "ok", TBL);
    Table.set("NucleusCount", row, nucCount, TBL);
    Table.set("Multinucleated", row, (nucCount >= 2) + 0, TBL);
    Table.set("Circularity", row, getResult("Circ.", m), TBL);
    Table.set("AspectRatio", row, getResult("AR", m), TBL);
    Table.set("Roundness", row, getResult("Round", m), TBL);
    Table.set("Solidity", row, getResult("Solidity", m), TBL);
    Table.set("Perimeter", row, getResult("Perim.", m), TBL);
    Table.set("Area_px", row, areaPx, TBL);
    Table.set("Area_um2", row, areaUm, TBL);
    Table.set("Heterochromatin_pct", row, hetPct, TBL);
    Table.set("Euchromatin_pct", row, 100 - hetPct, TBL);
    Table.set("Micronuclei_count", row, micro, TBL);
    Table.set("DenseBodies_count", row, dCount, TBL);
    Table.set("DenseBodies_near", row, dNear, TBL);
    Table.set("DenseBodies_far", row, dFar, TBL);
    Table.set("DenseBodies_area", row, dArea, TBL);
    Table.update(TBL);
    nDone = nDone + 1;

    if (SAVE_QC) saveQC(orig, best, name);
    closeAll();
}

// ---- micronucleus candidates: round bodies in the size window, away from main nucleus
function countMicronuclei(orig, bestRoi, calibrated, px) {
    selectWindow(orig);
    roiManager("select", bestRoi); Roi.getBounds(nx, ny, nw, nh);
    ncx = nx + nw/2; ncy = ny + nh/2; nrad = 0.5 * (nw + nh) / 2.0;

    selectWindow(orig); run("Select None");
    run("Duplicate...", "title=mwork");
    run("Gaussian Blur...", "sigma=2");
    setAutoThreshold("Otsu dark");
    run("Convert to Mask");
    run("Fill Holes");
    run("Set Measurements...", "area centroid shape redirect=None decimal=4");
    run("Analyze Particles...", "size=" + MICRO_MIN + "-" + MICRO_MAX +
        " circularity=" + MICRO_CIRC + "-1.00 show=Nothing clear");
    cnt = 0;
    micX = newArray(0); micY = newArray(0);
    for (r = 0; r < nResults; r++) {
        cx = getResult("X", r); cy = getResult("Y", r);
        cxp = cx; cyp = cy;
        if (calibrated) { cxp = cx / px; cyp = cy / px; }
        d = sqrt((cxp - ncx)*(cxp - ncx) + (cyp - ncy)*(cyp - ncy));
        if (d > 1.2 * nrad) {           // clearly outside the main nucleus
            cnt = cnt + 1;
            micX = Array.concat(micX, cxp); micY = Array.concat(micY, cyp);
        }
    }
    close();  // mwork
    // stash for QC overlay
    qc_micX = micX; qc_micY = micY;
    run("Set Measurements...", "area mean standard min centroid perimeter shape redirect=None decimal=4");
    return cnt;
}

// ---- distance map from the nucleus edge ----
var DMAP = "__nucDist";
function buildNucDistMap(orig, bestRoi, w, h) {
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    newImage(DMAP, "8-bit black", w, h, 1);
    selectWindow(orig); roiManager("select", bestRoi);
    selectWindow(DMAP);
    roiManager("select", bestRoi); setColor(255); fill(); run("Select None");
    setOption("BlackBackground", true); run("Invert"); run("Distance Map");
    selectWindow(orig);
}

// ---- dense cytoplasmic bodies (organelle screen) ----
var dg_count = 0; var dg_near = 0; var dg_far = 0; var dg_area = 0;
var qc_dX; var qc_dY; var qc_micX; var qc_micY;
function denseScan(orig, bestRoi, calibrated, px) {
    dg_count = 0; dg_near = 0; dg_far = 0; dg_area = 0;
    dX = newArray(0); dY = newArray(0);

    selectWindow(orig); run("Select None");
    run("Duplicate...", "title=dwork");
    // very dark structures stand out; slight blur removes speckle
    run("Gaussian Blur...", "sigma=1");
    setAutoThreshold("Default dark");
    run("Convert to Mask");
    run("Watershed");            // split touching bodies
    run("Set Measurements...", "area centroid redirect=None decimal=4");
    run("Analyze Particles...", "size=" + DENSE_MIN + "-" + DENSE_MAX +
        " circularity=0.20-1.00 show=Nothing clear");

    selectWindow(orig); roiManager("select", bestRoi); Roi.getBounds(nx, ny, nw, nh);

    for (r = 0; r < nResults; r++) {
        cx = getResult("X", r); cy = getResult("Y", r); a = getResult("Area", r);
        cxp = cx; cyp = cy; if (calibrated) { cxp = cx / px; cyp = cy / px; }
        // skip bodies inside the nucleus bounding box (those are chromatin, not organelles)
        insideNuc = (cxp >= nx && cxp <= nx+nw && cyp >= ny && cyp <= ny+nh);
        if (insideNuc) {
            selectWindow(orig); roiManager("select", bestRoi);
            if (Roi.contains(round(cxp), round(cyp))) continue;
        }
        dg_count = dg_count + 1; dg_area = dg_area + a;
        d = distEdge(cxp, cyp, px, calibrated);
        if (!isNaN(d) && d <= NEAR_FAR_UM) dg_near = dg_near + 1; else dg_far = dg_far + 1;
        dX = Array.concat(dX, cxp); dY = Array.concat(dY, cyp);
    }
    if (isOpen("dwork")) { selectWindow("dwork"); close(); }
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    qc_dX = dX; qc_dY = dY;
    run("Set Measurements...", "area mean standard min centroid perimeter shape redirect=None decimal=4");
}
function distEdge(cxp, cyp, px, calibrated) {
    if (!isOpen(DMAP)) return NaN;
    cur = getTitle(); selectWindow(DMAP);
    ix = round(cxp); iy = round(cyp);
    if (ix < 0) ix = 0; if (iy < 0) iy = 0;
    if (ix >= getWidth()) ix = getWidth() - 1;
    if (iy >= getHeight()) iy = getHeight() - 1;
    v = getPixel(ix, iy); selectWindow(cur);
    if (calibrated) return v * px;
    return v;
}

// ---- QC overlay ----
function saveQC(orig, bestRoi, name) {
    selectWindow(orig);
    run("Select None"); run("Remove Overlay");
    roiManager("select", bestRoi);
    Roi.setStrokeColor("yellow"); Roi.setStrokeWidth(4);
    run("Add Selection...");
    if (qc_micX.length > 0) { makeSelection("point", qc_micX, qc_micY);
        Roi.setStrokeColor("cyan"); run("Add Selection..."); }
    if (qc_dX.length > 0) { makeSelection("point", qc_dX, qc_dY);
        Roi.setStrokeColor("magenta"); run("Add Selection..."); }
    run("Flatten");
    saveAs("PNG", qcDir + noExt(name) + "_QC.png");
    close();
}

// ---- helpers ----
function areaUnit(calibrated, px) { return 1; }   // NUC_MIN already in the right units
function isTiff(name) { n = toLowerCase(name); return endsWith(n, ".tif") || endsWith(n, ".tiff"); }
function countTiffs(dir) {
    c = 0; list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) c = c + countTiffs(dir + list[i]);
        else if (isTiff(list[i])) c = c + 1;
    }
    return c;
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
    parts = split(dir, "/"); n = parts.length;
    if (n >= 1) return parts[n-1];
    return dir;
}
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
function closeAll() {
    roiManager("reset");
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    while (nImages > 0) { selectImage(nImages); close(); }
    run("Clear Results");
}
