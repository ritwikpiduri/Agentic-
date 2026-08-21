// ============================================================================
//  TEM_MembraneThickness.ijm  -- the by-hand part of Code 2.
//  Measures nuclear membrane / perinuclear-space thickness in nanometres.
//
//  It OPENS each image in a folder for you and asks: measure / skip / finish.
//  For an image you measure: draw a line along its SCALE BAR (+ its real
//  length), then draw short lines ACROSS the membrane. Each is logged to
//  MembraneThickness.csv. Because the scale comes from each image's own bar,
//  mixed magnifications are fine.
//
//  TIP: the membrane is only clear in some high-mag tiles, so the software
//  cannot know which to use. Easiest: copy the tiles that show the membrane
//  into one folder and point this macro there. (Or point it at everything and
//  click "Skip" on the ones without a visible membrane.)
// ============================================================================

inDir  = getDirectory("Choose the folder of images to review");
outDir = getDirectory("Choose a folder to save the results");
if (isOpen("Membrane")) { selectWindow("Membrane"); run("Close"); }
Table.create("Membrane");
setOption("ExpandableArrays", true);
run("Set Measurements...", "redirect=None decimal=4");

list = getFileList(inDir);
stop = false;
for (f = 0; f < list.length; f++) {
    if (stop) f = list.length;
    else if (isTiff(list[f])) {
        open(inDir + list[f]);
        name = getTitle();
        run("Enhance Contrast", "saturated=0.35");   // easier to see, view only

        Dialog.create("This image");
        Dialog.addMessage("Image: " + name);
        Dialog.addChoice("What do you want to do?",
            newArray("Measure it", "Skip it", "Finish now"), "Measure it");
        Dialog.show();
        choice = Dialog.getChoice();

        if (choice == "Finish now") { close(); stop = true; }
        else if (choice == "Skip it") { close(); }
        else {
            measureImage(name);
            close();
        }
    }
}

selectWindow("Membrane");
Table.save(outDir + "MembraneThickness.csv");
showMessage("Saved", "Saved: " + outDir + "MembraneThickness.csv\n" +
    Table.size("Membrane") + " measurements recorded.");

// ----------------------------------------------------------------------------
function measureImage(name) {
    // 1. set scale from this image's own scale bar
    setTool("line");
    waitForUser("Scale bar",
        "Draw a straight line EXACTLY along the SCALE BAR of this image,\n" +
        "then click OK.");
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

    // 2. draw lines across the membrane
    nLines = getNumber("How many thickness measurements on this image?", 5);
    i = 0;
    while (i < nLines) {
        waitForUser("Thickness " + (i + 1) + " of " + nLines,
            "Draw a short line ACROSS the nuclear membrane / perinuclear space\n" +
            "(inner membrane to outer membrane), then click OK.");
        if (selectionType() == 5) {
            lenUm = getValue("Length");
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
}
function isTiff(n) {
    n = toLowerCase(n);
    return endsWith(n, ".tif") || endsWith(n, ".tiff");
}
function folderName() {
    d = getDirectory("image");
    if (d == "") return "unknown";
    parts = split(d, "/");
    n = parts.length;
    if (n >= 2) return parts[n - 2] + "/" + parts[n - 1];
    if (n >= 1) return parts[n - 1];
    return d;
}
