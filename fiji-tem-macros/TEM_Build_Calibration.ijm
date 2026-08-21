// =============================================================================
//  TEM_Build_Calibration.ijm
//  Build a magnification -> pixel-size table ONCE, by measuring the scale bar
//  on one representative image per magnification. The other macros then read
//  the magnification from each filename and apply the matching pixel size
//  automatically -- so your mixed-magnification dataset calibrates itself.
//
//  You only do this a few times (once per distinct magnification), NOT per
//  image. Output: magnification_calibration.csv  (columns: magnification,
//  pixel_size_um). Put that CSV in your dataset root or output folder.
//
//  HOW TO USE
//   1. Run this macro. Choose where to save the CSV.
//   2. For each magnification: open a representative image (File > Open), draw a
//      straight line exactly along its scale bar, and enter the bar's real
//      length (e.g. 500 nm). The macro reads the magnification from the
//      filename (you can correct it) and records the pixel size.
//   3. Repeat for every magnification, then finish.
// =============================================================================

outDir = getDirectory("Choose where to save magnification_calibration.csv");
csvPath = outDir + "magnification_calibration.csv";

if (isOpen("Calibration")) { selectWindow("Calibration"); run("Close"); }
Table.create("Calibration");
setOption("ExpandableArrays", true);

unitChoices = newArray("nm", "micron", "A (angstrom)");

more = true;
while (more) {
    waitForUser("Open an image",
        "Open ONE representative image for a magnification you haven't done yet\n" +
        "(File > Open), then click OK. Cancel-free: just OK when it's open.");
    if (nImages == 0) { if (!getBoolean("No image open. Try again?")) more = false; continue; }

    name = getTitle();
    mag = parseMag(name);
    Dialog.create("Magnification");
    Dialog.addMessage("File: " + name);
    Dialog.addNumber("Magnification (auto-read; correct if wrong):", (mag > 0 ? mag : 0));
    Dialog.show();
    mag = Dialog.getNumber();

    setTool("line");
    waitForUser("Measure the scale bar",
        "Draw a straight line EXACTLY along the scale bar of this image,\n" +
        "then click OK.");
    if (selectionType() != 5) {   // not a straight line
        if (!getBoolean("No line drawn. Skip this image?")) continue;
        else continue;
    }
    getLine(x1, y1, x2, y2, lw);
    lenPx = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

    Dialog.create("Scale bar length");
    Dialog.addMessage("Line length = " + d2s(lenPx,1) + " pixels");
    Dialog.addNumber("Scale bar REAL length:", 500);
    Dialog.addChoice("Unit:", unitChoices, "nm");
    Dialog.show();
    known = Dialog.getNumber();
    unit  = Dialog.getChoice();

    knownUm = known;
    if (unit == "nm") knownUm = known / 1000.0;
    else if (unit == "A (angstrom)") knownUm = known / 10000.0;

    pxUm = knownUm / lenPx;    // microns per pixel

    r = Table.size("Calibration");
    Table.set("magnification", r, mag,   "Calibration");
    Table.set("pixel_size_um", r, pxUm,  "Calibration");
    Table.set("source_image",  r, name,  "Calibration");
    Table.update("Calibration");

    print("Mag " + mag + "x  ->  " + d2s(pxUm,6) + " um/pixel  (bar " +
          known + " " + unit + " = " + d2s(lenPx,1) + " px)");

    more = getBoolean("Recorded mag " + mag + "x. Add another magnification?");
}

// write a clean 2-column CSV the other macros read
txt = "magnification,pixel_size_um\n";
for (i = 0; i < Table.size("Calibration"); i++) {
    txt += Table.get("magnification", i, "Calibration") + "," +
           Table.get("pixel_size_um", i, "Calibration") + "\n";
}
File.saveString(txt, csvPath);
showMessage("Saved", "Calibration table written to:\n" + csvPath +
    "\n\nThe analysis macros will read the magnification from each filename\n" +
    "and apply the matching pixel size automatically.");

// -----------------------------------------------------------------------------
// Extract the magnification number from a filename. Handles tokens like
// 5000x, 5000X, x5000, 20kx, 5kx, Mag20000 ('k' means x1000). -1 if not found.
function parseMag(name) {
    s = name;
    if (matches(s, ".*[0-9]+[.]?[0-9]*[kK]?[xX].*")) {
        num = replace(s, ".*?([0-9]+[.]?[0-9]*)([kK]?)[xX].*", "$1");
        kfl = replace(s, ".*?([0-9]+[.]?[0-9]*)([kK]?)[xX].*", "$2");
        v = parseFloat(num);
        if (kfl == "k" || kfl == "K") v = v * 1000;
        return v;
    }
    if (matches(s, ".*[xX][0-9]+.*"))
        return parseFloat(replace(s, ".*[xX]([0-9]+).*", "$1"));
    return -1;
}
