// =============================================================================
//  TEM_Manual_Code2.ijm  = CODE 2, fully manual (no ilastik, no memory limits)
//  Run it on an OPEN high-magnification image. A menu lets you, per image:
//    1  Set scale (draw a line on the scale bar + type its real length)
//    2  Trace the nucleus (shape + counts multinucleation; sets the distance ref)
//    3  Mark organelles (ER/mito/Golgi/vacuole/lipid): size + distance + near/far
//    4  Mark micronuclei (small bodies outside the nucleus): size
//    5  Measure nuclear membrane / perinuclear thickness (draw short lines)
//    6  Finish this image (writes a per-image summary row)
//    7  Save all CSVs
//    0  Quit
//  Tables persist across images in one Fiji session. Save with [7] at the end.
//  Output: Nuclei.csv, Organelles.csv, Micronuclei.csv, Membrane.csv,
//          ImageSummary.csv  (works with the stats macro).
// =============================================================================

var NEAR_FAR_UM = 1.0;          // organelle centre <= this (um) from nucleus = "near"
var DMAP = "__nucDistMap";
var HAVE_DMAP = false;

initTables();

var go = true;
while (go) {
    if (nImages == 0) {
        showMessage("Open a high-mag image first, then run the macro again.");
        go = false;
    } else {
        Dialog.create("TEM manual - choose action");
        Dialog.addMessage("Image: " + getTitle());
        Dialog.addChoice("Action:", newArray(
            "1  Set scale (scale bar)",
            "2  Trace nucleus (shape + chromatin + multinucleation)",
            "3  Mark organelles (type, size, near/far)",
            "4  Mark micronuclei",
            "5  Membrane / perinuclear thickness",
            "6  Nuclear pores (count + density)",
            "7  Finish this image (summary row)",
            "8  Save all CSVs",
            "0  Quit"), "1  Set scale (scale bar)");
        Dialog.show();
        c = substring(Dialog.getChoice(), 0, 1);
        if      (c == "1") actSetScale();
        else if (c == "2") actNucleus();
        else if (c == "3") actOrganelles();
        else if (c == "4") actMicronuclei();
        else if (c == "5") actMembrane();
        else if (c == "6") actPores();
        else if (c == "7") actSummary();
        else if (c == "8") saveAll();
        else if (c == "0") go = false;
    }
}

// =============================================================================
function actSetScale() {
    setTool("line");
    waitForUser("Scale bar",
        "Draw a straight line EXACTLY along the scale bar of this image, then OK.");
    if (selectionType() != 5) { showMessage("No line drawn."); return; }
    getLine(x1, y1, x2, y2, lw);
    barPx = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    Dialog.create("Scale bar length");
    Dialog.addMessage("Bar = " + d2s(barPx,1) + " pixels.");
    Dialog.addNumber("Scale bar REAL length:", 500);
    Dialog.addChoice("Unit:", newArray("nm", "micron"), "nm");
    Dialog.show();
    known = Dialog.getNumber();
    unit = Dialog.getChoice();
    knownUm = known;
    if (unit == "nm") knownUm = known / 1000.0;
    run("Set Scale...", "distance=" + barPx + " known=" + knownUm + " pixel=1 unit=micron");
    run("Select None");
    showMessage("Scale set: " + d2s(knownUm/barPx, 6) + " um/pixel.");
}

function actNucleus() {
    if (!checkScale()) return;
    setMeas();
    setTool("freehand");
    waitForUser("Trace the NUCLEUS",
        "Outline the nuclear envelope with the freehand tool, then OK.\n" +
        "(Run this action again for each extra nucleus in the same cell.)");
    if (selectionType() < 0) { showMessage("No selection."); return; }
    run("Measure");
    row = nResults - 1;
    img = getTitle();
    n = countImg("Nuclei", img) + 1;
    r = Table.size("Nuclei");
    Table.set("Image", r, img, "Nuclei");
    Table.set("NucleusIndex", r, n, "Nuclei");
    Table.set("Area_um2", r, getResult("Area", row), "Nuclei");
    Table.set("Perimeter_um", r, getResult("Perim.", row), "Nuclei");
    Table.set("Circularity", r, getResult("Circ.", row), "Nuclei");
    Table.set("AspectRatio", r, getResult("AR", row), "Nuclei");
    Table.set("Roundness", r, getResult("Round", row), "Nuclei");
    Table.set("Solidity", r, getResult("Solidity", row), "Nuclei");
    // --- AUTO: heterochromatin (dark) vs euchromatin (light) inside the nucleus ---
    getStatistics(nArea, nMean, nMin, nMax, nStd);
    thr = round(nMean - 0.5 * nStd);
    Roi.getBounds(bx, by, bw, bh);
    dark = 0; tot = 0;
    for (yy = by; yy < by + bh; yy++) {
        for (xx = bx; xx < bx + bw; xx++) {
            if (Roi.contains(xx, yy)) {
                tot = tot + 1;
                if (getPixel(xx, yy) <= thr) dark = dark + 1;
            }
        }
    }
    het = NaN; if (tot > 0) het = 100.0 * dark / tot;
    Table.set("Heterochrom_pct", r, het, "Nuclei");
    Table.set("Euchrom_pct", r, 100 - het, "Nuclei");
    Table.update("Nuclei");
    buildDistMap();     // use this nucleus as the distance reference
    run("Select None");
    showMessage("Nucleus " + n + " recorded (circularity " +
        d2s(getResult("Circ.", row),3) + ")." +
        (n > 1 ? "\n>> This cell is multinucleated." : ""));
}

function actOrganelles() {
    if (!checkScale()) return;
    types = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody");
    Dialog.create("Organelle type");
    Dialog.addChoice("Marking which organelle?", types, types[0]);
    Dialog.addNumber("Near/far cutoff (um):", NEAR_FAR_UM);
    Dialog.show();
    otype = Dialog.getChoice();
    NEAR_FAR_UM = Dialog.getNumber();
    if (!HAVE_DMAP)
        showMessage("Tip: trace the nucleus first (action 2) so distances are recorded.");
    setMeas();
    setTool("freehand");
    img = getTitle();
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser(otype + " #" + i,
            "Trace this " + otype + " (freehand), then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure");
            row = nResults - 1;
            cx = getResult("X", row); cy = getResult("Y", row);
            dist = NaN; nf = "NA";
            if (HAVE_DMAP) { dist = distEdge(cx, cy); if (dist <= NEAR_FAR_UM) nf = "near"; else nf = "far"; }
            r = Table.size("Organelles");
            Table.set("Image", r, img, "Organelles");
            Table.set("Type", r, otype, "Organelles");
            Table.set("Index", r, i, "Organelles");
            Table.set("Area_um2", r, getResult("Area", row), "Organelles");
            Table.set("Dist_to_nucleus_um", r, dist, "Organelles");
            Table.set("NearFar", r, nf, "Organelles");
            Table.update("Organelles");
        }
        more = getBoolean("Mark another " + otype + "?");
    }
    run("Select None");
    showMessage(i + " " + otype + "(s) recorded.");
}

function actMicronuclei() {
    if (!checkScale()) return;
    setMeas();
    setTool("freehand");
    img = getTitle();
    i = 0; more = true;
    while (more) {
        i = i + 1;
        waitForUser("Micronucleus #" + i,
            "Trace a small chromatin body SEPARATE from the main nucleus, then OK.");
        if (selectionType() < 0) { i = i - 1; }
        else {
            run("Measure");
            row = nResults - 1;
            r = Table.size("Micronuclei");
            Table.set("Image", r, img, "Micronuclei");
            Table.set("Index", r, i, "Micronuclei");
            Table.set("Area_um2", r, getResult("Area", row), "Micronuclei");
            Table.set("Circularity", r, getResult("Circ.", row), "Micronuclei");
            Table.update("Micronuclei");
        }
        more = getBoolean("Mark another micronucleus?");
    }
    run("Select None");
    showMessage(i + " micronucleus/nuclei recorded.");
}

function actMembrane() {
    if (!checkScale()) return;
    img = getTitle();
    setTool("line");
    n = getNumber("How many thickness measurements on this image?", 5);
    i = 0;
    while (i < n) {
        waitForUser("Thickness " + (i+1) + " of " + n,
            "Draw a short line ACROSS the nuclear membrane / perinuclear space\n" +
            "(inner to outer), then OK.");
        if (selectionType() == 5) {
            lenUm = getValue("Length");
            r = Table.size("Membrane");
            Table.set("Image", r, img, "Membrane");
            Table.set("Measurement", r, i+1, "Membrane");
            Table.set("Thickness_nm", r, lenUm*1000.0, "Membrane");
            Table.set("Thickness_um", r, lenUm, "Membrane");
            Table.update("Membrane");
            i = i + 1;
        }
    }
    run("Select None");
    showMessage(n + " thickness measurements recorded.");
}

function actPores() {
    if (!checkScale()) return;
    img = getTitle();
    setTool("multipoint");
    waitForUser("Nuclear pores",
        "With the multi-point tool, CLICK on each nuclear pore along the\n" +
        "envelope, then click OK. (Every click = one pore.)");
    if (selectionType() != 10) { showMessage("Use the multi-point tool and click the pores."); return; }
    getSelectionCoordinates(xs, ys);
    count = xs.length;
    perim = nucPerim(img);          // total traced-nucleus envelope length (um)
    dens = NaN; if (perim > 0) dens = count / perim;
    r = Table.size("Pores");
    Table.set("Image", r, img, "Pores");
    Table.set("PoreCount", r, count, "Pores");
    Table.set("NucPerimeter_um", r, perim, "Pores");
    Table.set("Pores_per_um", r, dens, "Pores");
    Table.update("Pores");
    run("Select None");
    showMessage(count + " pores recorded." +
        (perim > 0 ? "\nDensity = " + d2s(dens,3) + " pores/um of envelope." :
        "\n(Trace the nucleus first for pore density.)"));
}

function actSummary() {
    img = getTitle();
    nNuc = countImg("Nuclei", img);
    r = Table.size("ImageSummary");
    Table.set("Image", r, img, "ImageSummary");
    Table.set("Sample", r, sampleOf(), "ImageSummary");
    Table.set("NucleusCount", r, nNuc, "ImageSummary");
    Table.set("Multinucleated", r, (nNuc > 1 ? 1 : 0), "ImageSummary");
    Table.set("MicronucleusCount", r, countImg("Micronuclei", img), "ImageSummary");
    Table.set("ER_count", r, countType(img, "ER"), "ImageSummary");
    Table.set("Mito_count", r, countType(img, "Mitochondria"), "ImageSummary");
    Table.set("Golgi_count", r, countType(img, "Golgi"), "ImageSummary");
    Table.set("Vacuole_count", r, countType(img, "Vacuole"), "ImageSummary");
    Table.set("Lipid_count", r, countType(img, "LipidBody"), "ImageSummary");
    Table.set("Pore_count", r, poreCount(img), "ImageSummary");
    Table.update("ImageSummary");
    showMessage("Summary written for " + img + " (" + nNuc + " nucleus/nuclei).");
}

function saveAll() {
    dir = getDirectory("Choose an output folder for the CSVs");
    saveTab("Nuclei", dir);
    saveTab("Organelles", dir);
    saveTab("Micronuclei", dir);
    saveTab("Membrane", dir);
    saveTab("Pores", dir);
    saveTab("ImageSummary", dir);
    showMessage("Saved all CSVs to:\n" + dir);
}
// total perimeter (um) of nuclei traced for this image, for pore density
function nucPerim(img) {
    n = Table.size("Nuclei"); s = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, "Nuclei") == img) s = s + Table.get("Perimeter_um", i, "Nuclei");
    return s;
}
function poreCount(img) {
    n = Table.size("Pores"); c = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, "Pores") == img) c = c + Table.get("PoreCount", i, "Pores");
    return c;
}

// ---- distance map from the nucleus ROI ----
function buildDistMap() {
    if (selectionType() < 0) return;
    orig = getTitle();
    getPixelSize(u, pw, ph);
    w = getWidth(); h = getHeight();
    if (isOpen(DMAP)) { selectWindow(DMAP); close(); }
    newImage(DMAP, "8-bit black", w, h, 1);
    selectWindow(orig);
    roiManager("reset");
    roiManager("add");
    selectWindow(DMAP);
    roiManager("select", 0);
    setColor(255); fill();
    run("Select None");
    setOption("BlackBackground", true);
    run("Invert");
    run("Distance Map");
    HAVE_DMAP = true;
    selectWindow(orig);
    roiManager("reset");
}
function distEdge(cx, cy) {
    if (!HAVE_DMAP || !isOpen(DMAP)) return NaN;
    cur = getTitle();
    selectWindow(DMAP);
    getPixelSize(u, pw, ph);
    px = round(cx / pw); py = round(cy / ph);
    if (px < 0) px = 0; if (py < 0) py = 0;
    if (px >= getWidth()) px = getWidth() - 1;
    if (py >= getHeight()) py = getHeight() - 1;
    d = getPixel(px, py) * pw;
    selectWindow(cur);
    return d;
}

// ---- helpers ----
function initTables() {
    mk("Nuclei"); mk("Organelles"); mk("Micronuclei"); mk("Membrane"); mk("Pores"); mk("ImageSummary");
}
function mk(name) { if (!isOpen(name)) Table.create(name); }
function setMeas() {
    run("Set Measurements...", "area mean centroid perimeter shape redirect=None decimal=3");
}
function checkScale() {
    getPixelSize(u, pw, ph);
    if (pw == 1 && (u == "pixel" || u == "pixels" || u == "")) {
        showMessage("Set the scale first (action 1) - draw the scale bar so sizes are in um.");
        return false;
    }
    return true;
}
function countImg(tbl, img) {
    n = Table.size(tbl); c = 0;
    for (i = 0; i < n; i++) if (Table.getString("Image", i, tbl) == img) c = c + 1;
    return c;
}
function countType(img, otype) {
    n = Table.size("Organelles"); c = 0;
    for (i = 0; i < n; i++)
        if (Table.getString("Image", i, "Organelles") == img &&
            Table.getString("Type", i, "Organelles") == otype) c = c + 1;
    return c;
}
function sampleOf() {
    d = getDirectory("image");
    if (d == "") return "unknown";
    p = split(d, "/"); n = p.length;
    if (n >= 2) return p[n-2] + "/" + p[n-1];
    if (n >= 1) return p[n-1];
    return d;
}
function saveTab(name, dir) {
    if (isOpen(name) && Table.size(name) > 0) { selectWindow(name); Table.save(dir + name + ".csv"); }
}
