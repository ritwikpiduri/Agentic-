// ============================================================================
//  TEM_ConvertArea.ijm  --  add a micron^2 area column to a Code 1 CSV.
//
//  Code 1 (TEM_OneClick) reported nucleus area in PIXELS because its
//  calibration table was empty. This reads that CSV, figures out the
//  magnification from each image's filename, and fills an Area_um2 column
//  using your pixel sizes.  No images, no re-clicking.
//
//  area (um^2) = area (px) * pixelsize_um * pixelsize_um
//
//  It writes the result next to your file as <name>_um2.csv
// ============================================================================

// ---- your pixel sizes (um/pixel) ----
var CAL_MAG = newArray(2000,   2600,     3400,    11000,    17500,    22000);
var CAL_PX  = newArray(0.0072, 0.005454, 0.00425, 0.001355, 0.000825, 0.0006479);
var CAM_CONST = 14.4;      // any other mag: um/px = 14.4 / magnification

path = File.openDialog("Select your Code 1 CSV (the one with Area_px)");
open(path);
T = File.getName(path);
N = Table.size(T);

// which column holds the pixel area? try common names.
areaCol = "Area_px";
if (!hasColumn(T, areaCol)) areaCol = "Area";       // fallback
if (!hasColumn(T, areaCol)) { showMessage("Could not find an 'Area_px' or 'Area' column."); exit; }

for (i = 0; i < N; i++) {
    img = Table.getString("Image", i, T);
    mag = readMag(img);
    px  = pxForMag(mag);
    apx = Table.get(areaCol, i, T);
    aum = NaN;
    if (px > 0 && !isNaN(apx)) aum = apx * px * px;
    Table.set("Magnification", i, mag, T);
    Table.set("PixelSize_um", i, px, T);
    Table.set("Area_um2", i, aum, T);
}
Table.update(T);

// save next to the original
dir = File.getParent(path);
base = noExt(File.getName(path));
outPath = dir + File.separator + base + "_um2.csv";
Table.save(outPath);
showMessage("Done", "Added Area_um2 (and Magnification, PixelSize_um).\nSaved:\n" + outPath);

// ---- helpers ----
function hasColumn(tbl, col) {
    headings = Table.headings(tbl);
    return indexOf(headings, col) >= 0;
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
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
