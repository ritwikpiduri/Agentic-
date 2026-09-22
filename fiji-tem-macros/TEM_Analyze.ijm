// ============================================================================
//  TEM_Analyze.ijm  --  reads your saved CSVs and DOES THE MATH for you.
//
//  Point it at ONE folder that contains your result CSVs (it searches inside
//  sub-folders too), and it produces a per-sample comparison you can paste
//  straight into your figures / paper.
//
//  It understands BOTH files you may have made:
//    - TEM_Organelle.csv       (near/far ovals, per-field counts, pores, and
//                               any nucleus rows if present)
//    - TEM_OrganelleCount.csv  (the standalone per-field counter)
//  It merges them by the "Sample" column (your folder name), so having one
//  CSV per sample folder is fine.
//
//  PER SAMPLE it reports (and saves to TEM_Analyze_Summary.csv):
//    For each organelle type (ER / Mito / Golgi / Vacuole / LipidBody / Lysosome):
//      _count      = how many you outlined with an oval (near/far session)
//      _near/_far  = of those, how many were near vs far from the nucleus
//      _dens/um2   = DENSITY from the per-field counter = total counted
//                    / total field area  (THIS is the fair "more or less?"
//                    number, because it corrects for how much area you looked at)
//      _fields     = how many fields contributed to that density
//    Nuclear pores: mean pores/um, and a pooled density (total pores / total
//                   envelope length measured).
//    Plus, if nucleus rows are present: circularity, solidity, aspect ratio,
//    heterochromatin %, membrane nm, % multinucleated, micronuclei/nucleus.
//
//  Raw counts are NOT comparable between samples (you image different amounts).
//  DENSITY (per um2 for organelles, per um for pores) is the comparable number.
// ============================================================================

var TYPES = newArray("ER", "Mitochondria", "Golgi", "Vacuole", "LipidBody", "Lysosome");

// ---- global row store (every row from every CSV, unified) ----
var rSample   = newArray(0);
var rRec      = newArray(0);   // Nucleus / Micronucleus / Membrane / Organelle / OrganelleCount / Pore
var rType     = newArray(0);   // organelle type ("" if n/a)
var rImage    = newArray(0);
var rNF       = newArray(0);   // near / far / ""
var rCount    = newArray(0);   // per-field count (OrganelleCount) OR pore count
var rArea     = newArray(0);   // organelle area, OR field area (OrganelleCount)
var rPerim    = newArray(0);   // envelope length (Pore)
var rCirc     = newArray(0);
var rAR       = newArray(0);
var rSol      = newArray(0);
var rHet      = newArray(0);
var rEu       = newArray(0);
var rThick    = newArray(0);
var rPoreDens = newArray(0);

srcDir = getDirectory("Choose the FOLDER that holds your result CSVs (subfolders are searched)");
csvs = findCsvs(srcDir);
if (csvs.length == 0) { showMessage("No result CSVs found under:\n" + srcDir); exit; }

setBatchMode(true);
loaded = 0;
for (i = 0; i < csvs.length; i++) if (loadCsv(csvs[i])) loaded = loaded + 1;
setBatchMode(false);

if (rSample.length == 0) { showMessage("Found CSVs but no usable rows."); exit; }

// unique samples, in first-seen order
samples = newArray(0);
for (i = 0; i < rSample.length; i++)
    if (indexOfArray(samples, rSample[i]) < 0) samples = Array.concat(samples, rSample[i]);

OUT = "TEM_Analyze_Summary";
if (isOpen(OUT)) { selectWindow(OUT); run("Close"); }
Table.create(OUT);

print("\\Clear");
print("========== TEM ANALYSIS  (" + loaded + " CSV file(s), " + rSample.length + " rows) ==========");

for (si = 0; si < samples.length; si++) {
    smp = samples[si];
    r = Table.size(OUT);
    Table.set("Sample", r, smp, OUT);

    print("");
    print("--- " + smp + " ---");

    // ---------- organelles: ovals (count + near/far) ----------
    for (t = 0; t < TYPES.length; t++) {
        ty = TYPES[t];
        cOval = 0; cNear = 0; cFar = 0; sArea = 0; nArea = 0;
        for (i = 0; i < rSample.length; i++) {
            if (rSample[i] != smp) continue;
            if (rRec[i] != "Organelle") continue;
            if (rType[i] != ty) continue;
            cOval = cOval + 1;
            if (rNF[i] == "near") cNear = cNear + 1;
            else if (rNF[i] == "far") cFar = cFar + 1;
            if (!isNaN(rArea[i])) { sArea = sArea + rArea[i]; nArea = nArea + 1; }
        }
        // ---------- organelles: per-field counter -> density ----------
        totCnt = 0; totField = 0; nFields = 0;
        for (i = 0; i < rSample.length; i++) {
            if (rSample[i] != smp) continue;
            if (rRec[i] != "OrganelleCount") continue;
            if (rType[i] != ty) continue;
            if (!isNaN(rCount[i])) totCnt = totCnt + rCount[i];
            if (!isNaN(rArea[i]) && rArea[i] > 0) { totField = totField + rArea[i]; nFields = nFields + 1; }
        }
        dens = NaN; if (totField > 0) dens = totCnt / totField;
        meanA = NaN; if (nArea > 0) meanA = sArea / nArea;

        Table.set(ty + "_count", r, cOval, OUT);
        Table.set(ty + "_near",  r, cNear, OUT);
        Table.set(ty + "_far",   r, cFar,  OUT);
        Table.set(ty + "_meanArea_um2", r, meanA, OUT);
        Table.set(ty + "_dens_per_um2", r, dens, OUT);
        Table.set(ty + "_fields", r, nFields, OUT);

        line = "  " + ty + ": ovals=" + cOval + " (near " + cNear + "/far " + cFar + ")";
        if (nFields > 0) line = line + "   density=" + d2s(dens,4) + " /um2 (" + totCnt + " in " + nFields + " fields)";
        print(line);
    }

    // ---------- pores ----------
    nPor = 0; sPorDens = 0; totPore = 0; totEnv = 0;
    for (i = 0; i < rSample.length; i++) {
        if (rSample[i] != smp) continue;
        if (rRec[i] != "Pore") continue;
        if (!isNaN(rPoreDens[i])) { nPor = nPor + 1; sPorDens = sPorDens + rPoreDens[i]; }
        if (!isNaN(rCount[i])) totPore = totPore + rCount[i];
        if (!isNaN(rPerim[i]) && rPerim[i] > 0) totEnv = totEnv + rPerim[i];
    }
    mPorDens = NaN; if (nPor > 0) mPorDens = sPorDens / nPor;
    pooledPor = NaN; if (totEnv > 0) pooledPor = totPore / totEnv;
    Table.set("Pores_per_um_mean", r, mPorDens, OUT);
    Table.set("Pores_per_um_pooled", r, pooledPor, OUT);
    Table.set("Pore_measurements", r, nPor, OUT);
    if (nPor > 0) print("  Pores: mean " + d2s(mPorDens,3) + " /um   pooled " + d2s(pooledPor,3) + " /um  (n=" + nPor + ")");

    // ---------- nucleus shape / chromatin (only if present) ----------
    nNuc = 0; sC = 0; sS = 0; sAR = 0; sHet = 0; sEu = 0;
    imgN = newArray(0); imgC = newArray(0);
    for (i = 0; i < rSample.length; i++) {
        if (rSample[i] != smp) continue;
        if (rRec[i] != "Nucleus") continue;
        nNuc = nNuc + 1;
        if (!isNaN(rCirc[i])) sC = sC + rCirc[i];
        if (!isNaN(rSol[i]))  sS = sS + rSol[i];
        if (!isNaN(rAR[i]))   sAR = sAR + rAR[i];
        if (!isNaN(rHet[i]))  sHet = sHet + rHet[i];
        if (!isNaN(rEu[i]))   sEu = sEu + rEu[i];
        k = indexOfArray(imgN, rImage[i]);
        if (k < 0) { imgN = Array.concat(imgN, rImage[i]); imgC = Array.concat(imgC, 1); }
        else imgC[k] = imgC[k] + 1;
    }
    if (nNuc > 0) {
        nMulti = 0; for (k = 0; k < imgC.length; k++) if (imgC[k] >= 2) nMulti = nMulti + 1;
        pctMulti = NaN; if (imgN.length > 0) pctMulti = 100.0 * nMulti / imgN.length;
        Table.set("Nuclei", r, nNuc, OUT);
        Table.set("Circularity_mean", r, sC/nNuc, OUT);
        Table.set("Solidity_mean", r, sS/nNuc, OUT);
        Table.set("AspectRatio_mean", r, sAR/nNuc, OUT);
        Table.set("Heterochrom_pct_mean", r, sHet/nNuc, OUT);
        Table.set("Euchrom_pct_mean", r, sEu/nNuc, OUT);
        Table.set("Pct_multinucleated", r, pctMulti, OUT);
        print("  Nuclei: " + nNuc + "  circ " + d2s(sC/nNuc,3) + "  solidity " + d2s(sS/nNuc,3) +
              "  AR " + d2s(sAR/nNuc,2) + "  het " + d2s(sHet/nNuc,1) + "%  multinuc " + d2s(pctMulti,1) + "%");
    }

    // ---------- membrane ----------
    nMem = 0; sMem = 0;
    for (i = 0; i < rSample.length; i++) {
        if (rSample[i] != smp) continue;
        if (rRec[i] != "Membrane") continue;
        if (!isNaN(rThick[i])) { nMem = nMem + 1; sMem = sMem + rThick[i]; }
    }
    if (nMem > 0) { Table.set("Membrane_nm_mean", r, sMem/nMem, OUT); print("  Membrane: " + d2s(sMem/nMem,1) + " nm (n=" + nMem + ")"); }

    // ---------- micronuclei ----------
    nMic = 0;
    for (i = 0; i < rSample.length; i++) { if (rSample[i] == smp && rRec[i] == "Micronucleus") nMic = nMic + 1; }
    if (nNuc > 0) Table.set("Micronuclei_per_nucleus", r, nMic/nNuc, OUT);
    if (nMic > 0) print("  Micronuclei: " + nMic);

    Table.update(OUT);
}

selectWindow(OUT);
Table.save(srcDir + "TEM_Analyze_Summary.csv");
print("");
print("Saved: " + srcDir + "TEM_Analyze_Summary.csv");
showMessage("Analysis done",
    "Compared " + samples.length + " sample(s) from " + loaded + " CSV file(s).\n" +
    "Saved: TEM_Analyze_Summary.csv\n(See the Log window for the readout.)");

// ============================================================================
// load one CSV into the global row store; returns true if it added rows
function loadCsv(path) {
    open(path);
    t = File.getName(path);
    if (!isOpen(t)) return false;
    n = Table.size(t);
    heads = Table.headings;
    added = false;

    hasRec  = hasCol(heads, "RecordType");
    hasCnt  = hasCol(heads, "Count");        // schema B (TEM_OrganelleCount.csv)
    hasFld  = hasCol(heads, "FieldArea_um2");

    for (i = 0; i < n; i++) {
        smp = getStr(t, "Sample", i, heads);
        img = getStr(t, "Image", i, heads);
        ty  = getStr(t, "Type", i, heads);

        if (hasRec) {
            // schema A: TEM_Organelle.csv (or a merged file with nucleus rows)
            rec = getStr(t, "RecordType", i, heads);
            cnt = NaN;
            if (rec == "OrganelleCount") cnt = getNum(t, "Index", i, heads);   // count stored in Index
            else if (rec == "Pore")      cnt = getNum(t, "PoreCount", i, heads);
            pushRow(smp, rec, ty, img,
                getStr(t, "NearFar", i, heads),
                cnt,
                getNum(t, "Area_um2", i, heads),
                getNum(t, "Perimeter_um", i, heads),
                getNum(t, "Circularity", i, heads),
                getNum(t, "AspectRatio", i, heads),
                getNum(t, "Solidity", i, heads),
                getNum(t, "Heterochrom_pct", i, heads),
                getNum(t, "Euchrom_pct", i, heads),
                getNum(t, "Thickness_nm", i, heads),
                getNum(t, "Pores_per_um", i, heads));
            added = true;
        }
        else if (hasCnt && hasFld) {
            // schema B: dedicated per-field counter
            pushRow(smp, "OrganelleCount", ty, img, "",
                getNum(t, "Count", i, heads),
                getNum(t, "FieldArea_um2", i, heads),
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN);
            added = true;
        }
    }
    selectWindow(t); run("Close");
    return added;
}

function pushRow(smp, rec, ty, img, nf, cnt, area, perim, circ, ar, sol, het, eu, thick, pdens) {
    rSample   = Array.concat(rSample, smp);
    rRec      = Array.concat(rRec, rec);
    rType     = Array.concat(rType, ty);
    rImage    = Array.concat(rImage, img);
    rNF       = Array.concat(rNF, nf);
    rCount    = Array.concat(rCount, cnt);
    rArea     = Array.concat(rArea, area);
    rPerim    = Array.concat(rPerim, perim);
    rCirc     = Array.concat(rCirc, circ);
    rAR       = Array.concat(rAR, ar);
    rSol      = Array.concat(rSol, sol);
    rHet      = Array.concat(rHet, het);
    rEu       = Array.concat(rEu, eu);
    rThick    = Array.concat(rThick, thick);
    rPoreDens = Array.concat(rPoreDens, pdens);
}

// ---- safe column access (returns "" / NaN if the column is absent) ----
function hasCol(heads, name) {
    parts = split(heads, "\t");
    for (i = 0; i < parts.length; i++) if (parts[i] == name) return true;
    return false;
}
function getStr(tbl, col, row, heads) {
    if (!hasCol(heads, col)) return "";
    return Table.getString(col, row, tbl);
}
function getNum(tbl, col, row, heads) {
    if (!hasCol(heads, col)) return NaN;
    return Table.get(col, row, tbl);
}

// ---- find result CSVs recursively, skipping outputs / backups ----
function findCsvs(dir) {
    out = newArray(0); l = getFileList(dir);
    for (i = 0; i < l.length; i++) {
        if (endsWith(l[i], "/")) out = Array.concat(out, findCsvs(dir + l[i]));
        else {
            nm = l[i]; low = toLowerCase(nm);
            if (!endsWith(low, ".csv")) continue;
            if (indexOf(low, "summary") >= 0) continue;
            if (indexOf(low, "analyze") >= 0) continue;
            if (indexOf(low, "_prev") >= 0) continue;
            out = Array.concat(out, dir + nm);
        }
    }
    return out;
}
function indexOfArray(arr, val) {
    for (i = 0; i < arr.length; i++) if (arr[i] == val) return i;
    return -1;
}
