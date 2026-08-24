// ============================================================================
//  TEM_Summary.ijm  --  reads TEM_Results.csv and COMPARES the 4 samples.
//  Run this AFTER you have collected data with TEM_Complete.ijm.
//
//  It prints, per sample (folder), and saves TEM_Summary.csv:
//    - Nuclei measured
//    - mean Circularity, Roundness, Solidity, Aspect Ratio
//    - mean Heterochromatin % / Euchromatin %
//    - % of imaged cells that are multinucleated
//    - Micronuclei per nucleus
//    - mean Membrane thickness (nm)
//    - Organelle counts BY TYPE (ER / Mito / Golgi / Vacuole / LipidBody)
//      with how many are NEAR vs FAR from the nucleus
//    - mean Nuclear pores per micron
//
//  This is what tells you "membrane is thinner in shLMNA", "more mitochondria
//  far from the nucleus in K32R", etc.
// ============================================================================

path = File.openDialog("Select your TEM_Results.csv");
open(path);                       // opens as a table
T = File.getName(path);           // table title = file name

// ---- collect unique samples ----
N = Table.size(T);
samples = newArray(0);
for (i = 0; i < N; i++) {
    s = Table.getString("Sample", i, T);
    if (indexOfArray(samples, s) < 0) samples = Array.concat(samples, s);
}

orgTypes = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody");

// ---- output table ----
OUT = "TEM_Summary";
if (isOpen(OUT)) { selectWindow(OUT); run("Close"); }
Table.create(OUT);

print("\\Clear");
print("========== TEM SUMMARY (by sample) ==========");

for (si = 0; si < samples.length; si++) {
    smp = samples[si];

    // nucleus shape + chromatin
    nNuc = 0; sCirc = 0; sRound = 0; sSol = 0; sAR = 0; sHet = 0; sEu = 0;
    // membrane
    nMem = 0; sMem = 0;
    // pores
    nPor = 0; sPor = 0;
    // micronuclei
    nMic = 0;

    // per-image nucleus counts for multinucleation
    imgNames = newArray(0); imgNuc = newArray(0);

    for (i = 0; i < N; i++) {
        if (Table.getString("Sample", i, T) != smp) continue;
        rec = Table.getString("RecordType", i, T);
        if (rec == "Nucleus") {
            nNuc = nNuc + 1;
            sCirc = sCirc + getNum("Circularity", i, T);
            sRound = sRound + getNum("Roundness", i, T);
            sSol = sSol + getNum("Solidity", i, T);
            sAR = sAR + getNum("AspectRatio", i, T);
            sHet = sHet + getNum("Heterochrom_pct", i, T);
            sEu = sEu + getNum("Euchrom_pct", i, T);
            img = Table.getString("Image", i, T);
            k = indexOfArray(imgNames, img);
            if (k < 0) { imgNames = Array.concat(imgNames, img); imgNuc = Array.concat(imgNuc, 1); }
            else { imgNuc[k] = imgNuc[k] + 1; }
        }
        else if (rec == "Membrane") { v = getNum("Thickness_nm", i, T); if (!isNaN(v)) { nMem = nMem + 1; sMem = sMem + v; } }
        else if (rec == "Pore")     { v = getNum("Pores_per_um", i, T); if (!isNaN(v)) { nPor = nPor + 1; sPor = sPor + v; } }
        else if (rec == "Micronucleus") { nMic = nMic + 1; }
    }

    // multinucleation: images with >= 2 nuclei
    nImg = imgNames.length; nMulti = 0;
    for (k = 0; k < imgNuc.length; k++) if (imgNuc[k] >= 2) nMulti = nMulti + 1;

    // organelle counts by type + near/far
    orgN = newArray(orgTypes.length);
    orgNear = newArray(orgTypes.length);
    orgFar = newArray(orgTypes.length);
    for (i = 0; i < N; i++) {
        if (Table.getString("Sample", i, T) != smp) continue;
        if (Table.getString("RecordType", i, T) != "Organelle") continue;
        ty = Table.getString("Type", i, T);
        ti = indexOfArray(orgTypes, ty);
        if (ti < 0) continue;
        orgN[ti] = orgN[ti] + 1;
        nf = Table.getString("NearFar", i, T);
        if (nf == "near") orgNear[ti] = orgNear[ti] + 1;
        else if (nf == "far") orgFar[ti] = orgFar[ti] + 1;
    }

    // averages (guard divide-by-zero)
    mCirc = avg(sCirc, nNuc); mRound = avg(sRound, nNuc); mSol = avg(sSol, nNuc);
    mAR = avg(sAR, nNuc); mHet = avg(sHet, nNuc); mEu = avg(sEu, nNuc);
    mMem = avg(sMem, nMem); mPor = avg(sPor, nPor);
    pctMulti = NaN; if (nImg > 0) pctMulti = 100.0 * nMulti / nImg;
    micPerNuc = NaN; if (nNuc > 0) micPerNuc = nMic / nNuc;

    // write summary row
    r = Table.size(OUT);
    Table.set("Sample", r, smp, OUT);
    Table.set("Nuclei_measured", r, nNuc, OUT);
    Table.set("Circularity_mean", r, mCirc, OUT);
    Table.set("Roundness_mean", r, mRound, OUT);
    Table.set("Solidity_mean", r, mSol, OUT);
    Table.set("AspectRatio_mean", r, mAR, OUT);
    Table.set("Heterochrom_pct_mean", r, mHet, OUT);
    Table.set("Euchrom_pct_mean", r, mEu, OUT);
    Table.set("Pct_multinucleated", r, pctMulti, OUT);
    Table.set("Micronuclei_per_nucleus", r, micPerNuc, OUT);
    Table.set("Membrane_nm_mean", r, mMem, OUT);
    Table.set("Pores_per_um_mean", r, mPor, OUT);
    for (t = 0; t < orgTypes.length; t++) {
        Table.set(orgTypes[t] + "_count", r, orgN[t], OUT);
        Table.set(orgTypes[t] + "_near", r, orgNear[t], OUT);
        Table.set(orgTypes[t] + "_far", r, orgFar[t], OUT);
    }
    Table.update(OUT);

    // plain-language print
    print("");
    print("--- " + smp + " ---");
    print("  Nuclei measured: " + nNuc);
    print("  Circularity: " + d2s(mCirc,3) + "   Roundness: " + d2s(mRound,3) + "   Solidity: " + d2s(mSol,3));
    print("  Heterochromatin: " + d2s(mHet,1) + "%   Euchromatin: " + d2s(mEu,1) + "%");
    print("  Multinucleated cells: " + d2s(pctMulti,1) + "%   Micronuclei/nucleus: " + d2s(micPerNuc,2));
    print("  Membrane thickness: " + d2s(mMem,1) + " nm   Pores/um: " + d2s(mPor,3));
    line = "  Organelles: ";
    for (t = 0; t < orgTypes.length; t++)
        line = line + orgTypes[t] + "=" + orgN[t] + "(near " + orgNear[t] + "/far " + orgFar[t] + ")  ";
    print(line);
}

// save summary next to the CSV
dir = File.getParent(path);
selectWindow(OUT);
Table.save(dir + File.separator + "TEM_Summary.csv");
print("");
print("Saved: " + dir + File.separator + "TEM_Summary.csv");
showMessage("Summary done", "Compared " + samples.length + " samples.\nSaved TEM_Summary.csv\n(See the Log window for the readout.)");

// ---- helpers ----
function getNum(col, row, tbl) { return Table.get(col, row, tbl); }
function avg(sum, n) { if (n > 0) return sum / n; return NaN; }
function indexOfArray(arr, val) {
    for (i = 0; i < arr.length; i++) if (arr[i] == val) return i;
    return -1;
}
