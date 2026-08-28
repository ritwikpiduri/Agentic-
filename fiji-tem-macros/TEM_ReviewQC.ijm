// ============================================================================
//  TEM_ReviewQC.ijm  --  keep/exclude Code 1 nuclei by looking at the outline.
//
//  It goes row by row through a Code 1 CSV, opens that image's QC picture
//  (the one with the yellow nucleus outline), and asks you ONE thing:
//     "Is the yellow outline on the NUCLEUS?"   Yes = keep,  No = exclude.
//
//  You decide from the OUTLINE, never from the circularity number -- so real
//  dysmorphic nuclei (low circularity, correct outline) are KEPT.
//
//  High-magnification images (no whole nucleus in frame) are auto-excluded
//  from nuclear-shape analysis -- you set the cutoff below.
//
//  It adds two columns and saves <name>_reviewed.csv (original untouched):
//     Keep        1 = use it, 0 = exclude
//     ReviewNote  why (kept / excluded_bad_outline / high_mag_skip / no_QC_image)
//
//  Afterwards: filter Keep = 1 for your nuclear-shape analysis. Press Cancel
//  any time to stop -- progress is saved after every row.
// ============================================================================

csvPath = File.openDialog("Select your Code 1 CSV");
open(csvPath);
T = File.getName(csvPath);
N = Table.size(T);

qcDir = getDirectory("Choose the QC folder (Code 1's yellow-outline pictures)");

Dialog.create("Review settings");
Dialog.addNumber("Only review images at/below this magnification:", 4000);
Dialog.addMessage("Whole nucleus is only seen at low mag. Images above this\nare auto-excluded from nuclear-shape analysis.");
Dialog.show();
magCut = Dialog.getNumber();

reviewed = 0; kept = 0;
for (i = 0; i < N; i++) {
    img = Table.getString("Image", i, T);
    mag = readMag(img);
    keep = 0; note = "";

    if (mag > magCut && mag > 0) { note = "high_mag_skip"; }
    else {
        qcPath = qcDir + noExt(img) + "_QC.png";
        if (!File.exists(qcPath)) { note = "no_QC_image"; }
        else {
            open(qcPath);
            setLocation(20, 20);
            k = getBoolean("Row " + (i+1) + " / " + N + "\n \n" +
                "Is the YELLOW outline on the NUCLEUS?\n \n" +
                "Yes = keep (even if it looks irregular / low circularity)\n" +
                "No  = outline is on the cell / a fragment / junk -> exclude");
            if (k) { keep = 1; note = "kept"; } else { keep = 0; note = "excluded_bad_outline"; }
            if (isOpen(getTitle())) close();
            reviewed = reviewed + 1;
        }
    }
    if (keep == 1) kept = kept + 1;
    Table.set("Keep", i, keep, T);
    Table.set("ReviewNote", i, note, T);
    Table.update(T);
    // save progress after each row
    outPath = File.getParent(csvPath) + File.separator + noExt(File.getName(csvPath)) + "_reviewed.csv";
    selectWindow(T); Table.save(outPath);
}

outPath = File.getParent(csvPath) + File.separator + noExt(File.getName(csvPath)) + "_reviewed.csv";
selectWindow(T); Table.save(outPath);
showMessage("Done",
    "Reviewed " + reviewed + " nuclei.\nKept: " + kept + "\n \nSaved: " + outPath +
    "\n \nNow filter Keep = 1 for your nuclear-shape analysis.");

// ---- helpers ----
function readMag(name) {
    s = name;
    if (matches(s, ".*[0-9]+[ ]?[kK]?[ ]?[xX].*")) {
        num = replace(s, ".*?([0-9]+)[ ]?([kK]?)[ ]?[xX].*", "$1");
        kf  = replace(s, ".*?([0-9]+)[ ]?([kK]?)[ ]?[xX].*", "$2");
        v = parseFloat(num); if (kf == "k" || kf == "K") v = v * 1000; return v;
    }
    return -1;
}
function noExt(name) { d = lastIndexOf(name, "."); if (d > 0) return substring(name, 0, d); return name; }
