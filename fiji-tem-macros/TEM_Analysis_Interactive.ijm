// =============================================================================
//  TEM_Analysis_Interactive.ijm
//  Semi-automated TEM cell/organelle analysis for Fiji (ImageJ) macro language
// -----------------------------------------------------------------------------
//  WHAT THIS DOES
//  A menu-driven assistant you run ON THE CURRENTLY OPEN TEM IMAGE. Each action
//  appends a row to a persistent results table (one per metric) that you save to
//  CSV at the end. It combines:
//    * automatable metrics  -> nuclear shape, chromatin (hetero/euchromatin)
//    * human-in-the-loop     -> organelle identification, perinuclear space,
//                               micronuclei (you draw/click, macro measures)
//
//  WHY SEMI-AUTOMATED
//  TEM is single-channel grayscale. ER / Golgi / mitochondria / lipid bodies /
//  vacuoles have no reliable universal intensity signature, so automatic
//  segmentation misclassifies them. This macro lets you identify them by eye and
//  does all the measuring, counting, distance and statistics for you. For true
//  auto-segmentation of organelles use Trainable Weka Segmentation / Labkit /
//  ilastik (see README) and feed the resulting masks into the batch macro.
//
//  CALIBRATION IS MANDATORY. Every quantitative result (µm, circularity is
//  scale-free but area/width/distance are not) depends on correct pixel size.
//  Action [1] sets it. Do it before anything else on each new magnification.
//
//  OUTPUT TABLES (saved as CSV via action [9]):
//    Nuclei.csv        one row per nucleus  (shape, area, multinucleation flag)
//    Chromatin.csv     one row per nucleus  (%heterochromatin, hetero/eu ratio)
//    Perinuclear.csv   one row per line     (nuclear-envelope gap width)
//    Organelles.csv    one row per organelle(type, area, distance to nucleus)
//    Micronuclei.csv   one row per micronuc.(area, circularity, distance)
//    ImageSummary.csv  one row per image    (roll-up counts for stats)
// =============================================================================

// ---------------------------------------------------------------------------
// GLOBAL CONFIG  (edit defaults here if you like)
// ---------------------------------------------------------------------------
var OUT_DIR          = "";     // where CSVs are written; asked on first save
var NEAR_FAR_UM      = 1.0;    // organelle centroid <= this (µm) from nuclear
                               // boundary is classified "near", else "far"
var MICRONUC_MAX_UM2 = 5.0;    // objects smaller than this treated as micronuclei
var TASK_LABELS = newArray(
    "1  Set / check scale (calibration)",
    "2  Add nucleus  (measure shape + multinucleation)",
    "3  Chromatin: heterochromatin vs euchromatin (needs a nucleus ROI)",
    "4  Perinuclear space width (draw lines)",
    "5  Organelles: ER / mito / Golgi / vacuole / lipid (click or trace)",
    "6  Micronuclei (trace small chromatin bodies)",
    "7  Finalise this image (write ImageSummary row)",
    "8  Show tables",
    "9  Save all tables to CSV",
    "0  Quit");

// tables live for the whole Fiji session so results accumulate across images
initTables();

// ---------------------------------------------------------------------------
// MAIN LOOP
// ---------------------------------------------------------------------------
keepGoing = true;
while (keepGoing) {
    if (nImages == 0) {
        showMessage("TEM Analysis",
            "Open a TEM image first, then run the macro again.");
        keepGoing = false;
    } else {
        Dialog.create("TEM Analysis - choose an action");
        Dialog.addMessage("Image: " + getTitle());
        Dialog.addChoice("Action:", TASK_LABELS, TASK_LABELS[1]);
        Dialog.show();
        choice = Dialog.getChoice();
        code = substring(choice, 0, 1);

        if      (code == "1") actionCalibrate();
        else if (code == "2") actionNucleus();
        else if (code == "3") actionChromatin();
        else if (code == "4") actionPerinuclear();
        else if (code == "5") actionOrganelles();
        else if (code == "6") actionMicronuclei();
        else if (code == "7") actionImageSummary();
        else if (code == "8") showTables();
        else if (code == "9") saveTables();
        else if (code == "0") keepGoing = false;
    }
}

// =============================================================================
//  ACTIONS
// =============================================================================

// --- [1] Calibration -------------------------------------------------------
function actionCalibrate() {
    getPixelSize(unit, pw, ph);
    Dialog.create("Calibration");
    Dialog.addMessage("Current: " + pw + " x " + ph + " " + unit + "/pixel");
    Dialog.addMessage("TEM pixel size is usually in the image metadata or a\n" +
                      "scale bar. Enter the size of ONE pixel, OR use a scale\n" +
                      "bar: measure its length in pixels with the line tool,\n" +
                      "then set known length below.");
    Dialog.addNumber("Pixel size (leave 0 to use scale-bar method):", 0);
    Dialog.addString("Unit:", "micron");
    Dialog.addMessage("--- scale-bar method (optional) ---");
    Dialog.addNumber("Scale bar length drawn (pixels):", 0);
    Dialog.addNumber("Scale bar real length:", 0);
    Dialog.show();
    px    = Dialog.getNumber();
    unit2 = Dialog.getString();
    barPx = Dialog.getNumber();
    barReal = Dialog.getNumber();

    if (px > 0) {
        setVoxelSize(px, px, 1, unit2);
    } else if (barPx > 0 && barReal > 0) {
        run("Set Scale...", "distance=" + barPx + " known=" + barReal +
            " pixel=1 unit=" + unit2);
    }
    getPixelSize(u2, pw2, ph2);
    showMessage("Calibration set", "Now: " + pw2 + " " + u2 + "/pixel");
}

// --- [2] Nucleus shape + multinucleation -----------------------------------
function actionNucleus() {
    requireCalibration();
    setStdMeasurements();
    if (selectionType() < 0) {
        setTool("freehand");
        waitForUser("Trace the NUCLEUS outline",
            "Use the freehand (or wand) tool to outline the nuclear envelope,\n" +
            "then click OK. Add ONE nucleus per click of this action; run the\n" +
            "action again for each additional nucleus in the same cell.");
    }
    if (selectionType() < 0) { showMessage("No selection made."); return; }

    roiManager("add");
    run("Measure");
    row = nResults - 1;
    img = getTitle();

    // multinucleation: how many nucleus ROIs already added for THIS image?
    nHere = countNucleiForImage(img) + 1;

    Table.set("Image",        Table.size("Nuclei"), img,                    "Nuclei");
    r = Table.size("Nuclei") - 1;
    Table.set("Folder",       r, imgFolder(),                               "Nuclei");
    Table.set("NucleusIndex", r, nHere,                                     "Nuclei");
    Table.set("Area_um2",     r, getResult("Area", row),                    "Nuclei");
    Table.set("Perimeter_um", r, getResult("Perim.", row),                  "Nuclei");
    Table.set("Circularity",  r, getResult("Circ.", row),                   "Nuclei");
    Table.set("AspectRatio",  r, getResult("AR", row),                      "Nuclei");
    Table.set("Roundness",    r, getResult("Round", row),                   "Nuclei");
    Table.set("Solidity",     r, getResult("Solidity", row),               "Nuclei");
    Table.set("MeanDensity",  r, getResult("Mean", row),                    "Nuclei");
    Table.update("Nuclei");
    run("Select None");
    showMessage("Nucleus #" + nHere + " recorded",
        "Circularity = " + d2s(getResult("Circ.", row), 3) +
        "\nArea = " + d2s(getResult("Area", row), 3) + " um^2" +
        "\nThis image now has " + nHere + " nucleus/nuclei." +
        (nHere > 1 ? "\n>> Multinucleated." : ""));
}

// --- [3] Chromatin: heterochromatin (dense/dark) vs euchromatin ------------
function actionChromatin() {
    requireCalibration();
    if (selectionType() < 0) {
        // let the user pick a previously stored nucleus ROI
        if (roiManager("count") == 0) {
            waitForUser("Select the nucleus",
                "Trace the nucleus (freehand) OR add it via action [2] first, " +
                "then run this action with that selection active.");
        } else {
            waitForUser("Select the nucleus ROI",
                "In the ROI Manager click the nucleus ROI you want to analyse, " +
                "then click OK.");
        }
    }
    if (selectionType() < 0) { showMessage("No nucleus selection."); return; }

    img = getTitle();
    // restrict analysis to inside the nucleus
    getStatistics(nucArea, nucMean, nucMin, nucMax, nucStd);

    // Heterochromatin is electron-dense = DARK in a standard TEM image.
    // Auto-threshold the dark pixels within the nuclear ROI.
    Dialog.create("Chromatin threshold");
    Dialog.addMessage("Nucleus mean gray = " + d2s(nucMean,1) +
                      "  (min " + nucMin + ", max " + nucMax + ")");
    Dialog.addMessage("Heterochromatin = darker (denser) pixels.\n" +
                      "Pixels with gray value <= threshold are counted as\n" +
                      "heterochromatin. Auto suggests a value; adjust if needed.");
    autoT = round(nucMean - 0.5 * nucStd);
    Dialog.addNumber("Heterochromatin threshold (<=):", autoT);
    Dialog.show();
    thr = Dialog.getNumber();

    // measure area fraction below threshold within the ROI, robust to bit depth
    getHistogram(values, counts, 256);
    // build a mask-free count by scanning ROI bounds
    Roi.getBounds(bx, by, bw, bh);
    heteroPix = 0; totalPix = 0;
    for (y = by; y < by + bh; y++) {
        for (x = bx; x < bx + bw; x++) {
            if (Roi.contains(x, y)) {
                v = getPixel(x, y);
                totalPix++;
                if (v <= thr) heteroPix++;
            }
        }
    }
    euPix = totalPix - heteroPix;
    heteroPct = 100.0 * heteroPix / totalPix;
    ratio = (euPix > 0) ? (heteroPix / euPix) : NaN;

    r = Table.size("Chromatin");
    Table.set("Image",            r, img,        "Chromatin");
    Table.set("Folder",           r, imgFolder(),"Chromatin");
    Table.set("NucleusMeanGray",  r, nucMean,    "Chromatin");
    Table.set("Threshold",        r, thr,        "Chromatin");
    Table.set("Heterochrom_pct",  r, heteroPct,  "Chromatin");
    Table.set("Euchrom_pct",      r, 100-heteroPct,"Chromatin");
    Table.set("Hetero_Eu_ratio",  r, ratio,      "Chromatin");
    Table.set("NucleusArea_um2",  r, nucArea,    "Chromatin");
    Table.update("Chromatin");
    showMessage("Chromatin recorded",
        "Heterochromatin = " + d2s(heteroPct,1) + " %\n" +
        "Euchromatin     = " + d2s(100-heteroPct,1) + " %\n" +
        "Hetero/Eu ratio = " + d2s(ratio,3) + "\n\n" +
        (heteroPct > 50 ? "More heterochromatin (condensed)." :
                          "More euchromatin (open)."));
    run("Select None");
}

// --- [4] Perinuclear space width -------------------------------------------
function actionPerinuclear() {
    requireCalibration();
    img = getTitle();
    setTool("line");
    n = getNumber("How many width measurements on this nucleus?", 5);
    for (i = 1; i <= n; i++) {
        waitForUser("Perinuclear space " + i + "/" + n,
            "Draw a straight line ACROSS the perinuclear space: from the inner\n" +
            "nuclear membrane to the outer nuclear membrane (the ER-continuous\n" +
            "gap). Click OK to record. High magnification needed.");
        if (selectionType() == 5) {   // straight line
            getLine(x1, y1, x2, y2, lw);
            len = lineLengthCalibrated(x1, y1, x2, y2);
            r = Table.size("Perinuclear");
            Table.set("Image",        r, img,        "Perinuclear");
            Table.set("Folder",       r, imgFolder(),"Perinuclear");
            Table.set("Measurement",  r, i,    "Perinuclear");
            Table.set("Width_um",     r, len,  "Perinuclear");
            Table.update("Perinuclear");
        } else {
            i--; // no line drawn, repeat this measurement
        }
    }
    run("Select None");
    showMessage("Perinuclear widths recorded for this nucleus.");
}

// --- [5] Organelles ---------------------------------------------------------
function actionOrganelles() {
    requireCalibration();
    img = getTitle();
    types = newArray("ER","Mitochondria","Golgi","Vacuole","LipidBody","Other");
    Dialog.create("Organelle type");
    Dialog.addChoice("Which organelle are you marking now?", types, types[0]);
    Dialog.addMessage("You will trace/click each one; the macro records count,\n" +
                      "area and distance to the nuclear boundary.\n" +
                      "'Near' = centroid within " + NEAR_FAR_UM + " um of the nucleus edge.");
    Dialog.addNumber("Near/Far cutoff (um):", NEAR_FAR_UM);
    Dialog.show();
    otype  = Dialog.getChoice();
    NEAR_FAR_UM = Dialog.getNumber();

    // need a nucleus ROI to measure distance from. Build/pick one.
    nucRoi = ensureNucleusForDistance();  // returns true if we have a distance map
    setStdMeasurements();

    setTool("freehand");
    idx = 0;
    more = true;
    while (more) {
        idx++;
        waitForUser(otype + " #" + idx,
            "Trace this " + otype + " (freehand), or use point tool for a rough\n" +
            "mark. Click OK to record. You'll be asked whether to continue.");
        if (selectionType() < 0) { idx--; }
        else {
            run("Measure");
            row = nResults - 1;
            cx = getResult("X", row);   // centroid, calibrated
            cy = getResult("Y", row);
            dist = NaN; nearFar = "NA";
            if (nucRoi) {
                dpx = distanceFromNucleusEdge(cx, cy);
                dist = dpx;
                nearFar = (dist <= NEAR_FAR_UM) ? "near" : "far";
            }
            r = Table.size("Organelles");
            Table.set("Image",        r, img,                       "Organelles");
            Table.set("Folder",       r, imgFolder(),               "Organelles");
            Table.set("Type",         r, otype,                     "Organelles");
            Table.set("Index",        r, idx,                       "Organelles");
            Table.set("Area_um2",     r, getResult("Area", row),    "Organelles");
            Table.set("MeanDensity",  r, getResult("Mean", row),    "Organelles");
            Table.set("Dist_to_nucleus_um", r, dist,                "Organelles");
            Table.set("NearFar",      r, nearFar,                   "Organelles");
            Table.update("Organelles");
        }
        m = getBoolean("Mark another " + otype + " in this image?");
        more = m;
    }
    run("Select None");
    showMessage(idx + " " + otype + "(s) recorded for this image.");
}

// --- [6] Micronuclei --------------------------------------------------------
function actionMicronuclei() {
    requireCalibration();
    img = getTitle();
    setStdMeasurements();
    ensureNucleusForDistance();
    setTool("freehand");
    idx = 0; more = true;
    while (more) {
        idx++;
        waitForUser("Micronucleus #" + idx,
            "Trace a micronucleus: a small, membrane-bound chromatin body\n" +
            "SEPARATE from the main nucleus (< " + MICRONUC_MAX_UM2 + " um^2 typ.).\n" +
            "Click OK to record.");
        if (selectionType() < 0) { idx--; }
        else {
            run("Measure");
            row = nResults - 1;
            cx = getResult("X", row); cy = getResult("Y", row);
            dist = distanceFromNucleusEdge(cx, cy);
            r = Table.size("Micronuclei");
            Table.set("Image",       r, img,                      "Micronuclei");
            Table.set("Folder",      r, imgFolder(),              "Micronuclei");
            Table.set("Index",       r, idx,                      "Micronuclei");
            Table.set("Area_um2",    r, getResult("Area", row),   "Micronuclei");
            Table.set("Circularity", r, getResult("Circ.", row),  "Micronuclei");
            Table.set("Dist_to_nucleus_um", r, dist,              "Micronuclei");
            Table.update("Micronuclei");
        }
        more = getBoolean("Mark another micronucleus?");
    }
    run("Select None");
    showMessage(idx + " micronucleus/nuclei recorded.");
}

// --- [7] Per-image roll-up --------------------------------------------------
function actionImageSummary() {
    img = getTitle();
    nNuc  = countNucleiForImage(img);
    nMic  = countRows("Micronuclei", img);
    nER   = countOrganelle(img, "ER");
    nMito = countOrganelle(img, "Mitochondria");
    nGol  = countOrganelle(img, "Golgi");
    nVac  = countOrganelle(img, "Vacuole");
    nLip  = countOrganelle(img, "LipidBody");

    r = Table.size("ImageSummary");
    Table.set("Image",            r, img,   "ImageSummary");
    Table.set("Sample",           r, guessSample(img), "ImageSummary");
    Table.set("NucleusCount",     r, nNuc,  "ImageSummary");
    Table.set("Multinucleated",   r, (nNuc > 1 ? 1 : 0), "ImageSummary");
    Table.set("MicronucleusCount",r, nMic,  "ImageSummary");
    Table.set("ER_count",         r, nER,   "ImageSummary");
    Table.set("Mito_count",       r, nMito, "ImageSummary");
    Table.set("Golgi_count",      r, nGol,  "ImageSummary");
    Table.set("Vacuole_count",    r, nVac,  "ImageSummary");
    Table.set("Lipid_count",      r, nLip,  "ImageSummary");
    Table.update("ImageSummary");
    showMessage("Image summary written for " + img);
}

// =============================================================================
//  HELPERS
// =============================================================================
function initTables() {
    ensureTable("Nuclei");
    ensureTable("Chromatin");
    ensureTable("Perinuclear");
    ensureTable("Organelles");
    ensureTable("Micronuclei");
    ensureTable("ImageSummary");
}
function ensureTable(name) {
    if (!isOpen(name)) { Table.create(name); }
}
function setStdMeasurements() {
    run("Set Measurements...",
        "area mean min centroid perimeter shape redirect=None decimal=3");
}
function requireCalibration() {
    normalizeToMicron();                 // use embedded nm/Angstrom scale as um
    getPixelSize(unit, pw, ph);
    if (pw == 1 && (unit == "pixel" || unit == "pixels" || unit=="")) {
        showMessage("Not calibrated",
            "Image is in pixels (the scale is not in the metadata - maybe it's\n" +
            "only a drawn scale bar). Run action [1] to set the scale first, or\n" +
            "areas/distances will be in pixels, not microns.");
    }
}
// convert an embedded length-unit scale to microns; leave pixels/unknown as-is
function normalizeToMicron() {
    getPixelSize(unit, pw, ph);
    u = toLowerCase(unit);
    if (u == "micron" || u == "microns" || u == "um" || u == "µm") return;
    if (u == "nm" || u == "nanometer" || u == "nanometre" || u == "nanometers")
        setVoxelSize(pw/1000.0, ph/1000.0, 1, "micron");
    else if (u == "a" || u == "angstrom" || u == "ang" || u == "å")
        setVoxelSize(pw/10000.0, ph/10000.0, 1, "micron");
    else if (u == "mm" || u == "millimeter")
        setVoxelSize(pw*1000.0, ph*1000.0, 1, "micron");
}
function lineLengthCalibrated(x1, y1, x2, y2) {
    getPixelSize(unit, pw, ph);
    dx = (x2 - x1) * pw;
    dy = (y2 - y1) * ph;
    return sqrt(dx*dx + dy*dy);
}
function countNucleiForImage(img) {
    return countRows("Nuclei", img);
}
function countRows(tbl, img) {
    n = Table.size(tbl); c = 0;
    for (i = 0; i < n; i++) if (Table.getString("Image", i, tbl) == img) c++;
    return c;
}
function countOrganelle(img, otype) {
    n = Table.size("Organelles"); c = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, "Organelles") == img &&
            Table.getString("Type",  i, "Organelles") == otype) c++;
    return c;
}
// derive a sample/group label from the file path or name; edit to taste
function guessSample(img) {
    return imgFolder();   // condition/cell path of the open image
}
// last two path components of the OPEN image's directory (e.g. "Control/cell07")
// so the stats macro can match your condition keyword to the folder path.
function imgFolder() {
    d = getDirectory("image");
    if (d == "") return "unknown";
    parts = split(d, "/");
    n = parts.length;
    if (n >= 2) return parts[n-2] + "/" + parts[n-1];
    if (n == 1) return parts[0];
    return d;
}

// ---- distance-to-nucleus machinery ----
// We build a binary mask of the nucleus, compute a Euclidean distance map on
// the OUTSIDE, so reading the map at an organelle centroid gives the calibrated
// distance from the nuclear boundary.
var DMAP_TITLE = "__nucDistMap";
var HAVE_DMAP  = false;

function ensureNucleusForDistance() {
    if (isOpen(DMAP_TITLE) && HAVE_DMAP) return true;
    if (selectionType() < 0) {
        ok = getBoolean("A nucleus outline is needed to measure organelle\n" +
                        "distances. Do you want to trace it now?\n" +
                        "(No = distances will be blank.)");
        if (!ok) return false;
        setTool("freehand");
        waitForUser("Trace the nucleus", "Outline the nuclear envelope, click OK.");
        if (selectionType() < 0) return false;
    }
    buildNucleusDistanceMap();
    return HAVE_DMAP;
}
function buildNucleusDistanceMap() {
    orig = getTitle();
    getPixelSize(unit, pw, ph);
    w = getWidth(); h = getHeight();
    // create an 8-bit mask: nucleus = 255
    newImage(DMAP_TITLE, "8-bit black", w, h, 1);
    selectWindow(orig);
    if (selectionType() >= 0) {
        roiManager("add"); ri = roiManager("count") - 1;
        selectWindow(DMAP_TITLE);
        roiManager("select", ri);
        setColor(255); fill();
        roiManager("select", ri); roiManager("delete");
        run("Select None");
    }
    // distance map OUTSIDE the nucleus: invert so background(outside)=object
    setOption("BlackBackground", true);
    run("Options...", "iterations=1 count=1 black");
    // EDM measures distance to nearest background pixel; we want distance from
    // nucleus edge for outside pixels -> compute EDM of the inverted mask.
    selectWindow(DMAP_TITLE);
    run("Invert");            // now outside=255, nucleus=0
    run("Distance Map");      // each outside pixel = pixels to nearest nucleus px
    HAVE_DMAP = true;
    selectWindow(orig);
}
// returns calibrated distance (um) from nucleus edge at calibrated point (cx,cy)
function distanceFromNucleusEdge(cx, cy) {
    if (!HAVE_DMAP || !isOpen(DMAP_TITLE)) return NaN;
    getPixelSize(unit, pw, ph);
    px = round(cx / pw); py = round(cy / ph);
    cur = getTitle();
    selectWindow(DMAP_TITLE);
    if (px < 0) px = 0; if (py < 0) py = 0;
    if (px >= getWidth())  px = getWidth()-1;
    if (py >= getHeight()) py = getHeight()-1;
    dpx = getPixel(px, py);          // distance in pixels (0 if inside nucleus)
    selectWindow(cur);
    return dpx * pw;                 // approximate calibrated distance
}

function showTables() {
    tabs = newArray("Nuclei","Chromatin","Perinuclear","Organelles",
                    "Micronuclei","ImageSummary");
    for (i = 0; i < tabs.length; i++)
        if (isOpen(tabs[i])) { selectWindow(tabs[i]); }
}
function saveTables() {
    if (OUT_DIR == "") OUT_DIR = getDirectory("Choose an output folder for CSVs");
    tabs = newArray("Nuclei","Chromatin","Perinuclear","Organelles",
                    "Micronuclei","ImageSummary");
    for (i = 0; i < tabs.length; i++) {
        if (isOpen(tabs[i]) && Table.size(tabs[i]) > 0) {
            selectWindow(tabs[i]);
            Table.save(OUT_DIR + tabs[i] + ".csv");
        }
    }
    showMessage("Saved", "CSV tables written to:\n" + OUT_DIR);
}
