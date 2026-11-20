// ============================================================================
//  LabPack_Summary.ijm  --  compares conditions from LabPack_Nuclear.csv files.
//
//  Point it at ONE folder that holds your LabPack_Nuclear.csv file(s) (it also
//  searches sub-folders, so one CSV per condition-folder is fine). It merges
//  them by the "Sample" column (your folder name) and prints + saves, per
//  condition:
//    - Nuclei measured
//    - mean Area, Circularity, Aspect Ratio, Solidity
//    - mean Heterochromatin % / Euchromatin % / MeanGrey
//    - % of imaged cells that are multinucleated (images with >= 2 nuclei)
//    - Micronuclei per nucleus
//    - mean Perinuclear space (nm)
//
//  Output: LabPack_Summary.csv (+ a readable list in the Log window).
//  Raw per-nucleus comparison and stats are still best done in Excel/Prism;
//  this just gives you the per-condition means quickly.
// ============================================================================

srcDir = getDirectory("Choose the FOLDER holding your LabPack_Nuclear.csv file(s)");
csvs = findCsvs(srcDir);
if (csvs.length == 0) { showMessage("No LabPack_Nuclear.csv found under:\n" + srcDir); exit; }

// global row store
var rSample = newArray(0);
var rImage  = newArray(0);
var rRec    = newArray(0);
var rArea   = newArray(0);
var rCirc   = newArray(0);
var rAR     = newArray(0);
var rSol    = newArray(0);
var rHet    = newArray(0);
var rEu     = newArray(0);
var rMG     = newArray(0);
var rPeri   = newArray(0);

setBatchMode(true);
loaded = 0;
for (i = 0; i < csvs.length; i++) if (loadCsv(csvs[i])) loaded = loaded + 1;
setBatchMode(false);
if (rSample.length == 0) { showMessage("CSV found but no usable rows."); exit; }

samples = newArray(0);
for (i = 0; i < rSample.length; i++)
    if (idxOf(samples, rSample[i]) < 0) samples = Array.concat(samples, rSample[i]);

OUT = "LabPack_Summary";
if (isOpen(OUT)) { selectWindow(OUT); run("Close"); }
Table.create(OUT);

print("\\Clear");
print("========== LabPack SUMMARY  (" + loaded + " file(s)) ==========");

for (si = 0; si < samples.length; si++) {
    smp = samples[si];
    nNuc = 0; sA = 0; sC = 0; sAR = 0; sS = 0; sHet = 0; sEu = 0; sMG = 0;
    nMem = 0; sMem = 0; nMic = 0;
    imgN = newArray(0); imgC = newArray(0);

    for (i = 0; i < rSample.length; i++) {
        if (rSample[i] != smp) continue;
        if (rRec[i] == "Nucleus") {
            nNuc = nNuc + 1;
            sA = sA + rArea[i]; sC = sC + rCirc[i]; sAR = sAR + rAR[i]; sS = sS + rSol[i];
            if (!isNaN(rHet[i])) sHet = sHet + rHet[i];
            if (!isNaN(rEu[i]))  sEu = sEu + rEu[i];
            if (!isNaN(rMG[i]))  sMG = sMG + rMG[i];
            k = idxOf(imgN, rImage[i]);
            if (k < 0) { imgN = Array.concat(imgN, rImage[i]); imgC = Array.concat(imgC, 1); }
            else imgC[k] = imgC[k] + 1;
        }
        else if (rRec[i] == "Perinuclear") { if (!isNaN(rPeri[i])) { nMem = nMem + 1; sMem = sMem + rPeri[i]; } }
        else if (rRec[i] == "Micronucleus") { nMic = nMic + 1; }
    }

    nMulti = 0; for (k = 0; k < imgC.length; k++) if (imgC[k] >= 2) nMulti = nMulti + 1;
    pctMulti = NaN; if (imgN.length > 0) pctMulti = 100.0 * nMulti / imgN.length;
    micPer = NaN; if (nNuc > 0) micPer = nMic / nNuc;

    r = Table.size(OUT);
    Table.set("Sample", r, smp, OUT);
    Table.set("Nuclei_measured", r, nNuc, OUT);
    Table.set("Area_um2_mean", r, avg(sA, nNuc), OUT);
    Table.set("Circularity_mean", r, avg(sC, nNuc), OUT);
    Table.set("AspectRatio_mean", r, avg(sAR, nNuc), OUT);
    Table.set("Solidity_mean", r, avg(sS, nNuc), OUT);
    Table.set("Heterochrom_pct_mean", r, avg(sHet, nNuc), OUT);
    Table.set("Euchrom_pct_mean", r, avg(sEu, nNuc), OUT);
    Table.set("MeanGrey_mean", r, avg(sMG, nNuc), OUT);
    Table.set("Pct_multinucleated", r, pctMulti, OUT);
    Table.set("Micronuclei_per_nucleus", r, micPer, OUT);
    Table.set("Perinuclear_nm_mean", r, avg(sMem, nMem), OUT);
    Table.update(OUT);

    print("");
    print("--- " + smp + " ---");
    print("  Nuclei: " + nNuc + "   Area " + d2s(avg(sA,nNuc),3) + " um2   Circ " + d2s(avg(sC,nNuc),3) +
          "   AR " + d2s(avg(sAR,nNuc),2) + "   Solidity " + d2s(avg(sS,nNuc),3));
    print("  Heterochromatin " + d2s(avg(sHet,nNuc),1) + "%   Euchromatin " + d2s(avg(sEu,nNuc),1) +
          "%   MeanGrey " + d2s(avg(sMG,nNuc),1));
    print("  Multinucleated " + d2s(pctMulti,1) + "%   Micronuclei/nucleus " + d2s(micPer,2) +
          "   Perinuclear " + d2s(avg(sMem,nMem),1) + " nm");
}

selectWindow(OUT);
Table.save(srcDir + "LabPack_Summary.csv");
print("");
print("Saved: " + srcDir + "LabPack_Summary.csv");
showMessage("Summary done", "Compared " + samples.length + " condition(s).\nSaved LabPack_Summary.csv\n(See the Log window.)");

// ============================================================================
function loadCsv(path) {
    open(path); t = File.getName(path);
    if (!isOpen(t)) return false;
    n = Table.size(t); heads = Table.headings; added = false;
    for (i = 0; i < n; i++) {
        rSample = Array.concat(rSample, getStr(t, "Sample", i, heads));
        rImage  = Array.concat(rImage, getStr(t, "Image", i, heads));
        rRec    = Array.concat(rRec, getStr(t, "RecordType", i, heads));
        rArea   = Array.concat(rArea, getNum(t, "Area_um2", i, heads));
        rCirc   = Array.concat(rCirc, getNum(t, "Circularity", i, heads));
        rAR     = Array.concat(rAR, getNum(t, "AspectRatio", i, heads));
        rSol    = Array.concat(rSol, getNum(t, "Solidity", i, heads));
        rHet    = Array.concat(rHet, getNum(t, "Heterochrom_pct", i, heads));
        rEu     = Array.concat(rEu, getNum(t, "Euchrom_pct", i, heads));
        rMG     = Array.concat(rMG, getNum(t, "MeanGrey", i, heads));
        rPeri   = Array.concat(rPeri, getNum(t, "PerinuclearSpace_nm", i, heads));
        added = true;
    }
    selectWindow(t); run("Close");
    return added;
}
function hasCol(heads, name) {
    parts = split(heads, "\t");
    for (i = 0; i < parts.length; i++) if (parts[i] == name) return true;
    return false;
}
function getStr(tbl, col, row, heads) { if (!hasCol(heads, col)) return ""; return Table.getString(col, row, tbl); }
function getNum(tbl, col, row, heads) { if (!hasCol(heads, col)) return NaN; return Table.get(col, row, tbl); }
function findCsvs(dir) {
    out = newArray(0); l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        if (endsWith(l[i], "/")) out = Array.concat(out, findCsvs(dir + l[i]));
        else {
            low = toLowerCase(l[i]);
            if (!endsWith(low, ".csv")) continue;
            if (indexOf(low, "summary") >= 0) continue;
            if (indexOf(low, "_prev") >= 0) continue;
            out = Array.concat(out, dir + l[i]);
        }
    }
    return out;
}
function avg(sum, n) { if (n > 0) return sum / n; return NaN; }
function idxOf(arr, val) { for (i = 0; i < arr.length; i++) if (arr[i] == val) return i; return -1; }
