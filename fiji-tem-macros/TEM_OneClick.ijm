// ============================================================================
//  TEM_OneClick.ijm  -- everything automatic in ONE macro.
//  Run it, pick your image folder, pick an output folder. That's all.
//
//  It goes through every .tif/.tiff (including subfolders) and measures, per
//  nucleus: circularity, aspect ratio, roundness, solidity, area, perimeter,
//  and heterochromatin vs euchromatin %.
//
//  SCALE: circularity, shape and chromatin % do NOT need a scale and are always
//  correct. Area/perimeter need a pixel size. If your images have the scale in
//  their metadata it is used automatically. Otherwise fill MAG/PX below (one
//  line per magnification) and it converts area to microns; if you leave them
//  empty, area is reported in PIXELS (still fine for comparing between groups).
// ============================================================================

// ---- optional: magnification -> microns-per-pixel (leave as-is to skip) ----
// Example: MAG = newArray(2000, 20000);  PX = newArray(0.0055, 0.00055);
var MAG = newArray();     // magnifications, read from the filename (e.g. 20000)
var PX  = newArray();     // matching pixel size in microns per pixel
// ---------------------------------------------------------------------------

var HET_FACTOR = 0.5;     // heterochromatin = pixels darker than mean-0.5*SD
var MIN_NUC    = 2.0;     // ignore blobs smaller than this (microns^2, or px if uncalibrated)
var SAVE_QC    = true;    // save an outline image per nucleus so you can check

var nTotal = 0;
var nDone  = 0;

setBatchMode(true);
inDir  = getDirectory("Choose the folder with your TEM images");
outDir = getDirectory("Choose a folder to save the results");
qcDir  = outDir + "QC_overlays/";
if (SAVE_QC) File.makeDirectory(qcDir);

if (isOpen("Results_TEM")) { selectWindow("Results_TEM"); run("Close"); }
Table.create("Results_TEM");
run("Set Measurements...", "area mean standard min centroid perimeter shape redirect=None decimal=4");

nTotal = countTiffs(inDir);
print("\\Clear");
print("Found " + nTotal + " images. Starting...");
walk(inDir);

selectWindow("Results_TEM");
Table.save(outDir + "TEM_Results.csv");
setBatchMode(false);
showMessage("Done!",
    "Measured " + nDone + " of " + nTotal + " images.\n" +
    "Saved: " + outDir + "TEM_Results.csv\n" +
    (SAVE_QC ? "Check outlines in: " + qcDir : ""));

// ============================================================================
function walk(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/"))
            walk(dir + list[i]);
        else if (isTiff(list[i]))
            processOne(dir + list[i], dir, list[i]);
    }
}

function processOne(path, dir, name) {
    print("[" + (nDone + 1) + "/" + nTotal + "] " + name);
    showProgress(nDone, nTotal);
    open(path);
    if (bitDepth() == 24) run("8-bit");

    // --- figure out the scale ---
    scaleSource = "pixels";
    px = 1;
    if (hasMetadataScale()) { getPixelSize(u, px, ph2); scaleSource = "metadata"; }
    else {
        mag = readMag(name);
        p = pxForMag(mag);
        if (p > 0) { px = p; setVoxelSize(px, px, 1, "micron"); scaleSource = "mag=" + mag; }
    }
    calibrated = (scaleSource != "pixels");

    orig = getTitle();
    w = getWidth();
    h = getHeight();

    // --- find the nucleus (largest smooth dark-ish region) ---
    run("Duplicate...", "title=work");
    run("Gaussian Blur...", "sigma=3");
    setAutoThreshold("Otsu dark");
    run("Convert to Mask");
    run("Options...", "iterations=3 count=1 black");
    run("Fill Holes");
    roiManager("reset");
    minSize = MIN_NUC;
    run("Analyze Particles...", "size=" + minSize + "-Infinity circularity=0.1-1.00 show=Nothing add");

    best = -1;
    bestA = 0;
    frameA = w * h;
    if (calibrated) frameA = w * h * px * px;
    for (r = 0; r < roiManager("count"); r++) {
        roiManager("select", r);
        getStatistics(a);
        if (a > bestA && a < 0.95 * frameA) { bestA = a; best = r; }
    }

    row = Table.size("Results_TEM");
    fld = folderOf(dir);
    Table.set("Image", row, name, "Results_TEM");
    Table.set("Folder", row, fld, "Results_TEM");
    Table.set("ScaleSource", row, scaleSource, "Results_TEM");

    if (best < 0) {
        Table.set("Status", row, "no_nucleus_found", "Results_TEM");
        Table.update("Results_TEM");
        closeAll();
        return;
    }

    selectWindow(orig);
    roiManager("select", best);
    run("Measure");
    m = nResults - 1;
    area = getResult("Area", m);
    getStatistics(dummyA, nucMean, nucMin, nucMax, nucSD);

    areaPx = area;
    areaUm = NaN;
    if (calibrated) { areaUm = area; areaPx = area / (px * px); }

    // --- chromatin: dark = heterochromatin, light = euchromatin ---
    thr = round(nucMean - HET_FACTOR * nucSD);
    Roi.getBounds(bx, by, bw, bh);
    dark = 0;
    tot = 0;
    for (y = by; y < by + bh; y++) {
        for (x = bx; x < bx + bw; x++) {
            if (Roi.contains(x, y)) {
                tot = tot + 1;
                if (getPixel(x, y) <= thr) dark = dark + 1;
            }
        }
    }
    hetPct = NaN;
    if (tot > 0) hetPct = 100.0 * dark / tot;

    Table.set("Status", row, "ok", "Results_TEM");
    Table.set("Circularity", row, getResult("Circ.", m), "Results_TEM");
    Table.set("AspectRatio", row, getResult("AR", m), "Results_TEM");
    Table.set("Roundness", row, getResult("Round", m), "Results_TEM");
    Table.set("Solidity", row, getResult("Solidity", m), "Results_TEM");
    Table.set("Perimeter", row, getResult("Perim.", m), "Results_TEM");
    Table.set("Area_px", row, areaPx, "Results_TEM");
    Table.set("Area_um2", row, areaUm, "Results_TEM");
    Table.set("Heterochromatin_pct", row, hetPct, "Results_TEM");
    Table.set("Euchromatin_pct", row, 100 - hetPct, "Results_TEM");
    Table.update("Results_TEM");
    nDone = nDone + 1;

    if (SAVE_QC) {
        selectWindow(orig);
        roiManager("select", best);
        run("Flatten");
        saveAs("PNG", qcDir + noExt(name) + "_QC.png");
        close();
    }
    closeAll();
}

// ---- helpers ----
function isTiff(name) {
    n = toLowerCase(name);
    return endsWith(n, ".tif") || endsWith(n, ".tiff");
}
function countTiffs(dir) {
    c = 0;
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/"))
            c = c + countTiffs(dir + list[i]);
        else if (isTiff(list[i]))
            c = c + 1;
    }
    return c;
}
function hasMetadataScale() {
    getPixelSize(unit, pw, ph);
    u = toLowerCase(unit);
    if (u == "micron" || u == "microns" || u == "um") return (pw != 1);
    if (u == "nm" || u == "nanometer" || u == "nanometers") { setVoxelSize(pw/1000.0, ph/1000.0, 1, "micron"); return true; }
    if (u == "a" || u == "angstrom") { setVoxelSize(pw/10000.0, ph/10000.0, 1, "micron"); return true; }
    return false;
}
function readMag(name) {
    s = name;
    if (matches(s, ".*[0-9]+[kK]?[xX].*")) {
        num = replace(s, ".*?([0-9]+)([kK]?)[xX].*", "$1");
        kf = replace(s, ".*?([0-9]+)([kK]?)[xX].*", "$2");
        v = parseFloat(num);
        if (kf == "k" || kf == "K") v = v * 1000;
        return v;
    }
    if (matches(s, ".*[xX][0-9]+.*"))
        return parseFloat(replace(s, ".*[xX]([0-9]+).*", "$1"));
    return -1;
}
function pxForMag(mag) {
    for (i = 0; i < MAG.length; i++)
        if (MAG[i] == mag) return PX[i];
    return -1;
}
function folderOf(dir) {
    parts = split(dir, "/");
    n = parts.length;
    if (n >= 2) return parts[n-2] + "/" + parts[n-1];
    if (n == 1) return parts[0];
    return dir;
}
function noExt(name) {
    d = lastIndexOf(name, ".");
    if (d > 0) return substring(name, 0, d);
    return name;
}
function closeAll() {
    roiManager("reset");
    while (nImages > 0) { selectImage(nImages); close(); }
    run("Clear Results");
}
