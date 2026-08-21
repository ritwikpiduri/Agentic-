// ============================================================================
//  TEM_MembraneThickness.ijm  -- the by-hand part of Code 2.
//  Measures nuclear membrane / perinuclear-space thickness in nanometres.
//
//  This one is manual by necessity: the gap between the inner and outer nuclear
//  membrane is only tens of nm wide, so you draw the line across it -- no
//  software finds it reliably. The macro handles the scale and the maths and
//  logs every measurement to MembraneThickness.csv.
//
//  For EACH high-magnification image (where the envelope is clearly visible):
//    1. draw a line along that image's own SCALE BAR + type its real length,
//    2. draw short lines ACROSS the membrane (inner -> outer); each is recorded.
//  Because you set the scale from each image's own bar, mixed magnifications are
//  fine -- every image calibrates itself.
// ============================================================================

var rows = 0;
outDir = getDirectory("Choose a folder to save the results");
if (isOpen("Membrane")) { selectWindow("Membrane"); run("Close"); }
Table.create("Membrane");
setOption("ExpandableArrays", true);
run("Set Measurements...", "redirect=None decimal=4");

more = true;
while (more) {
    if (nImages == 0) {
        waitForUser("Open an image",
            "Open a high-magnification image (File > Open) where you can see the\n" +
            "nuclear membrane, then click OK.");
        if (nImages == 0) { more = getBoolean("No image open. Keep going?"); continue; }
    }
    name = getTitle();

    // ---- 1. set the scale from this image's scale bar ----
    setTool("line");
    waitForUser("Scale bar",
        "Draw a straight line EXACTLY along the SCALE BAR of this image,\n" +
        "then click OK.  (If this image is already calibrated, just click OK.)");
    if (selectionType() == 5) {
        getLine(x1, y1, x2, y2, lw);
        barPx = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        Dialog.create("Scale bar length");
        Dialog.addMessage("Bar = " + d2s(barPx, 1) + " pixels long.");
        Dialog.addNumber("Scale bar REAL length:", 500);
        Dialog.addChoice("Unit:", newArray("nm", "micron"), "nm");
        Dialog.show();
        realLen = Dialog.getNumber();
        unit = Dialog.getChoice();
        knownUm = realLen;
        if (unit == "nm") knownUm = realLen / 1000.0;
        run("Set Scale...", "distance=" + barPx + " known=" + knownUm + " pixel=1 unit=micron");
    }

    // ---- 2. measure membrane thickness lines ----
    nLines = getNumber("How many thickness measurements on this image?", 5);
    i = 0;
    while (i < nLines) {
        waitForUser("Thickness " + (i + 1) + " of " + nLines,
            "Draw a short line ACROSS the nuclear membrane / perinuclear space\n" +
            "(from the inner membrane to the outer membrane), then click OK.");
        if (selectionType() == 5) {
            lenUm = getValue("Length");     // calibrated length in microns
            lenNm = lenUm * 1000.0;
            fld = folderName();
            r = Table.size("Membrane");
            Table.set("Image", r, name, "Membrane");
            Table.set("Folder", r, fld, "Membrane");
            Table.set("Measurement", r, i + 1, "Membrane");
            Table.set("Thickness_nm", r, lenNm, "Membrane");
            Table.set("Thickness_um", r, lenUm, "Membrane");
            Table.update("Membrane");
            i = i + 1;
        }
    }
    run("Select None");
    more = getBoolean("Done with this image.\nOpen another image and continue?  (No = finish & save)");
}

selectWindow("Membrane");
Table.save(outDir + "MembraneThickness.csv");
showMessage("Saved", "Saved: " + outDir + "MembraneThickness.csv");

function folderName() {
    d = getDirectory("image");
    if (d == "") return "unknown";
    parts = split(d, "/");
    n = parts.length;
    if (n >= 2) return parts[n - 2] + "/" + parts[n - 1];
    if (n >= 1) return parts[n - 1];
    return d;
}
