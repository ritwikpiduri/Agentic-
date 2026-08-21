// =============================================================================
//  TEM_Stats.ijm   -- statistics for the TEM CSVs, PURE IMAGEJ (no Python)
// -----------------------------------------------------------------------------
//  Reads the CSV tables written by the other macros and tests whether each
//  metric differs between your samples/conditions. Everything runs inside Fiji.
//
//  TESTS (implemented from scratch in macro code):
//    * numeric metrics  -> Mann-Whitney U (non-parametric; good for EM counts &
//                          non-normal data), normal approximation with tie
//                          correction -> two-sided p-value.
//    * proportions      -> 2x2 chi-square (df=1) -> p-value.
//    * >2 groups        -> all pairwise Mann-Whitney comparisons.
//  p-values use an error-function approximation (Abramowitz & Stegun 7.1.26).
//
//  GROUPING: you type your condition keywords (e.g. "Control,Treated"). Each
//  row is assigned to the first keyword that appears anywhere in its
//  Image/Folder/Sample text. So it works whether the condition is in the folder
//  path or the filename. Rows matching no keyword are ignored.
//
//  OUTPUT: a readable report in the Log window + two tables you can save:
//    TEM_Descriptives (n/mean/median/SD per group) and TEM_Tests (p-values).
//
//  USAGE: Plugins > Macros > Run... > this file. Point it at the CSV folder.
// =============================================================================

var GROUPS = newArray(0);          // condition keywords (filled from dialog)
var ALPHA  = 0.05;

// dialog
csvDir = getDirectory("Choose the folder containing the TEM CSV files");
Dialog.create("TEM statistics - define your groups");
Dialog.addString("Condition keywords (comma-separated):", "Control,Treated", 40);
Dialog.addMessage("Each data row is assigned to the FIRST keyword found in its\n" +
                  "Image / Folder / Sample text. Case-insensitive. Rows that\n" +
                  "match no keyword are skipped.");
Dialog.addNumber("Significance threshold (alpha):", 0.05);
Dialog.show();
gstr = Dialog.getString();
ALPHA = Dialog.getNumber();
GROUPS = splitTrim(gstr, ",");

freshTable("TEM_Descriptives");
freshTable("TEM_Tests");
print("\\Clear");
print("=== TEM statistics ===");
print("Groups: " + gstr);
print("CSV folder: " + csvDir);
print("");

// ---- run each metric if its CSV is present ----
runNumeric(csvDir + "Nuclei.csv",      "Circularity",      "Nuclear circularity");
runNumeric(csvDir + "Nuclei.csv",      "Area_um2",         "Nuclear area (um^2)");
runNumeric(csvDir + "Nuclei.csv",      "Solidity",         "Nuclear solidity");

runNumeric(csvDir + "TEM_Batch_NucleusChromatin.csv", "Circularity",     "Nuclear circularity (batch)");
runNumeric(csvDir + "TEM_Batch_NucleusChromatin.csv", "Nucleus_Area_um2","Nuclear area (batch)");
runNumeric(csvDir + "TEM_Batch_NucleusChromatin.csv", "Heterochrom_pct", "Heterochromatin % (batch)");

runNumeric(csvDir + "Chromatin.csv",   "Heterochrom_pct",  "Heterochromatin %");
runNumeric(csvDir + "Chromatin.csv",   "Hetero_Eu_ratio",  "Hetero/Eu ratio");

runNumeric(csvDir + "Perinuclear.csv", "Width_um",         "Perinuclear space width (um)");

runProportion(csvDir + "ImageSummary.csv", "Multinucleated",   "Multinucleated cells (fraction)");
runNumeric(csvDir + "ImageSummary.csv",    "MicronucleusCount","Micronuclei per image");
runNumeric(csvDir + "ImageSummary.csv",    "ER_count",         "ER per image");
runNumeric(csvDir + "ImageSummary.csv",    "Mito_count",       "Mitochondria per image");
runNumeric(csvDir + "ImageSummary.csv",    "Golgi_count",      "Golgi per image");
runNumeric(csvDir + "ImageSummary.csv",    "Vacuole_count",    "Vacuoles per image");
runNumeric(csvDir + "ImageSummary.csv",    "Lipid_count",      "Lipid bodies per image");

runNumeric(csvDir + "Organelles.csv",  "Dist_to_nucleus_um","Organelle distance to nucleus (um)");
runNearFar(csvDir + "Organelles.csv",  "Organelles near nucleus (fraction)");

// save
selectWindow("TEM_Descriptives"); Table.save(csvDir + "TEM_Descriptives.csv");
selectWindow("TEM_Tests");        Table.save(csvDir + "TEM_Tests.csv");
print("");
print("Saved: TEM_Descriptives.csv and TEM_Tests.csv in the CSV folder.");
print("NOTE: many metrics tested -> apply a multiple-comparison correction");
print("(e.g. Benjamini-Hochberg) before reporting individual p-values.");
showMessage("Statistics done",
    "Report is in the Log window.\nTEM_Descriptives.csv and TEM_Tests.csv saved.");

// =============================================================================
//  METRIC RUNNERS
// =============================================================================
function runNumeric(path, col, label) {
    if (!File.exists(path)) return;
    if (!loadValues(path, col)) return;          // fills VALS[], GRP[]
    if (VALS.length == 0) { print("[skip] " + label + ": no matching rows"); return; }
    print("-- " + label + " --");
    present = groupsPresent();
    // descriptives
    for (g = 0; g < present.length; g++) {
        vals = valuesForGroup(present[g]);
        if (vals.length == 0) continue;
        Array.getStatistics(vals, mn, mx, mean, sd);
        med = median(vals);
        print("   " + pad(present[g],14) + " n=" + vals.length +
              "  mean=" + d2s(mean,4) + "  median=" + d2s(med,4) + "  sd=" + d2s(sd,4));
        addDesc(label, present[g], vals.length, mean, med, sd);
    }
    // tests: all pairwise
    for (a = 0; a < present.length; a++) {
        for (b = a + 1; b < present.length; b++) {
            va = valuesForGroup(present[a]);
            vb = valuesForGroup(present[b]);
            if (va.length < 3 || vb.length < 3) continue;
            p = mannWhitney(va, vb);
            reportTest(label, present[a] + " vs " + present[b], "Mann-Whitney", p);
        }
    }
    print("");
}

function runProportion(path, col, label) {
    if (!File.exists(path)) return;
    if (!loadValues(path, col)) return;
    if (VALS.length == 0) { print("[skip] " + label + ": no matching rows"); return; }
    print("-- " + label + " --");
    present = groupsPresent();
    succ = newArray(present.length);
    tot  = newArray(present.length);
    for (g = 0; g < present.length; g++) {
        vals = valuesForGroup(present[g]);
        s = 0;
        for (i = 0; i < vals.length; i++) if (vals[i] != 0) s++;
        succ[g] = s; tot[g] = vals.length;
        frac = (vals.length > 0) ? s / vals.length : 0;
        print("   " + pad(present[g],14) + " " + s + "/" + vals.length +
              "  fraction=" + d2s(frac,4));
        addDesc(label, present[g], vals.length, frac, frac, NaN);
    }
    for (a = 0; a < present.length; a++)
        for (b = a + 1; b < present.length; b++) {
            p = chi2x2(succ[a], tot[a], succ[b], tot[b]);
            reportTest(label, present[a] + " vs " + present[b], "chi-square(1)", p);
        }
    print("");
}

// special: Organelles NearFar column -> 1 if "near"
function runNearFar(path, label) {
    if (!File.exists(path)) return;
    if (!loadNearFar(path)) return;
    if (VALS.length == 0) { print("[skip] " + label + ": no matching rows"); return; }
    print("-- " + label + " --");
    present = groupsPresent();
    succ = newArray(present.length);
    tot  = newArray(present.length);
    for (g = 0; g < present.length; g++) {
        vals = valuesForGroup(present[g]);
        s = 0; for (i = 0; i < vals.length; i++) if (vals[i] != 0) s++;
        succ[g] = s; tot[g] = vals.length;
        frac = (vals.length > 0) ? s / vals.length : 0;
        print("   " + pad(present[g],14) + " " + s + "/" + vals.length +
              " near  fraction=" + d2s(frac,4));
        addDesc(label, present[g], vals.length, frac, frac, NaN);
    }
    for (a = 0; a < present.length; a++)
        for (b = a + 1; b < present.length; b++) {
            p = chi2x2(succ[a], tot[a], succ[b], tot[b]);
            reportTest(label, present[a] + " vs " + present[b], "chi-square(1)", p);
        }
    print("");
}

// =============================================================================
//  DATA LOADING  (fills globals VALS[] numeric and GRP[] strings, aligned)
// =============================================================================
var VALS = newArray(0);
var GRP  = newArray(0);

function loadValues(path, col) {
    title = File.getName(path);
    Table.open(path);
    heads = Table.headings;
    if (indexOf("\t" + heads + "\t", "\t" + col + "\t") < 0) {
        print("[skip] " + col + " not found in " + title);
        selectWindow(title); run("Close"); return false;
    }
    n = Table.size;
    VALS = newArray(0); GRP = newArray(0);
    for (r = 0; r < n; r++) {
        v = Table.get(col, r);
        if (isNaN(v)) continue;
        g = assignGroup(rowText(r));
        if (g == "") continue;
        VALS = Array.concat(VALS, v);
        GRP  = Array.concat(GRP, g);
    }
    selectWindow(title); run("Close");
    return true;
}

function loadNearFar(path) {
    title = File.getName(path);
    Table.open(path);
    heads = Table.headings;
    if (indexOf("\t" + heads + "\t", "\tNearFar\t") < 0) {
        selectWindow(title); run("Close"); return false;
    }
    n = Table.size;
    VALS = newArray(0); GRP = newArray(0);
    for (r = 0; r < n; r++) {
        nf = Table.getString("NearFar", r);
        if (nf == "NA" || nf == "") continue;
        g = assignGroup(rowText(r));
        if (g == "") continue;
        VALS = Array.concat(VALS, (nf == "near") ? 1 : 0);
        GRP  = Array.concat(GRP, g);
    }
    selectWindow(title); run("Close");
    return true;
}

// concatenate the identifier columns of the active table's row r into one
// searchable string so group keywords can match a folder OR a filename
function rowText(r) {
    s = "";
    s = s + safeStr("Image", r);
    s = s + " " + safeStr("Folder", r);
    s = s + " " + safeStr("Sample", r);
    return s;
}
function safeStr(col, r) {
    heads = Table.headings;
    if (indexOf("\t" + heads + "\t", "\t" + col + "\t") < 0) return "";
    return Table.getString(col, r);
}

// =============================================================================
//  GROUPING
// =============================================================================
function assignGroup(text) {
    lt = toLowerCase(text);
    for (i = 0; i < GROUPS.length; i++)
        if (indexOf(lt, toLowerCase(GROUPS[i])) >= 0) return GROUPS[i];
    return "";
}
function groupsPresent() {
    out = newArray(0);
    for (i = 0; i < GROUPS.length; i++) {
        cnt = 0;
        for (j = 0; j < GRP.length; j++) if (GRP[j] == GROUPS[i]) cnt++;
        if (cnt > 0) out = Array.concat(out, GROUPS[i]);
    }
    return out;
}
function valuesForGroup(g) {
    out = newArray(0);
    for (j = 0; j < GRP.length; j++) if (GRP[j] == g) out = Array.concat(out, VALS[j]);
    return out;
}

// =============================================================================
//  STATISTICS
// =============================================================================
// Mann-Whitney U, normal approx with tie correction -> two-sided p
function mannWhitney(a, b) {
    n1 = a.length; n2 = b.length; N = n1 + n2;
    all = Array.concat(a, b);
    ranks = rankAvg(all);            // average ranks, ties handled
    R1 = 0;
    for (i = 0; i < n1; i++) R1 += ranks[i];
    U1 = R1 - n1 * (n1 + 1) / 2.0;
    mu = n1 * n2 / 2.0;
    // tie correction term
    tieSum = tieCorrection(all);
    sigma = sqrt( (n1 * n2 / 12.0) * ((N + 1) - tieSum / (N * (N - 1.0))) );
    if (sigma == 0) return NaN;
    z = (abs(U1 - mu) - 0.5) / sigma;    // continuity correction
    if (z < 0) z = 0;
    return 2.0 * (1.0 - normCDF(z));
}
// average ranks of values (1..N), ties get mean rank
function rankAvg(values) {
    n = values.length;
    ranks = newArray(n);
    order = Array.rankPositions(values);   // ascending sort indices
    i = 0;
    while (i < n) {
        j = i;
        while (j + 1 < n && values[order[j+1]] == values[order[i]]) j++;
        avg = ((i + 1) + (j + 1)) / 2.0;
        for (k = i; k <= j; k++) ranks[order[k]] = avg;
        i = j + 1;
    }
    return ranks;
}
function tieCorrection(values) {
    n = values.length;
    order = Array.rankPositions(values);
    sum = 0; i = 0;
    while (i < n) {
        j = i;
        while (j + 1 < n && values[order[j+1]] == values[order[i]]) j++;
        t = j - i + 1;
        sum += t*t*t - t;
        i = j + 1;
    }
    return sum;
}
// 2x2 chi-square (df=1) -> p; groups (s1/n1) vs (s2/n2)
function chi2x2(s1, n1, s2, n2) {
    f1 = n1 - s1; f2 = n2 - s2;         // failures
    N = n1 + n2;
    rowSucc = s1 + s2; rowFail = f1 + f2;
    if (N == 0 || rowSucc == 0 || rowFail == 0) return NaN;
    // expected
    e11 = n1 * rowSucc / N; e12 = n1 * rowFail / N;
    e21 = n2 * rowSucc / N; e22 = n2 * rowFail / N;
    chi = sq(s1 - e11)/e11 + sq(f1 - e12)/e12 + sq(s2 - e21)/e21 + sq(f2 - e22)/e22;
    // df=1: P(X>chi) = 2*(1 - Phi(sqrt(chi)))
    return 2.0 * (1.0 - normCDF(sqrt(chi)));
}
function sq(x) { return x * x; }
// standard normal CDF via erf approximation (A&S 7.1.26)
function normCDF(z) { return 0.5 * (1.0 + erf(z / sqrt(2.0))); }
function erf(x) {
    s = 1; if (x < 0) { s = -1; x = -x; }
    t = 1.0 / (1.0 + 0.3275911 * x);
    y = 1.0 - (((((1.061405429*t - 1.453152027)*t) + 1.421413741)*t
              - 0.284496736)*t + 0.254829592) * t * exp(-x*x);
    return s * y;
}
function median(a) {
    s = Array.copy(a); Array.sort(s); n = s.length;
    if (n == 0) return NaN;
    if (n % 2 == 1) return s[(n-1)/2];
    return (s[n/2 - 1] + s[n/2]) / 2.0;
}

// =============================================================================
//  REPORTING / UTILITIES
// =============================================================================
function reportTest(metric, comparison, test, p) {
    sig = (!isNaN(p) && p < ALPHA) ? "  *SIGNIFICANT*" : "";
    print("   " + test + "  " + comparison + "  p=" + d2s(p,5) + sig);
    r = Table.size("TEM_Tests");
    Table.set("Metric",     r, metric,     "TEM_Tests");
    Table.set("Comparison", r, comparison, "TEM_Tests");
    Table.set("Test",       r, test,       "TEM_Tests");
    Table.set("p_value",    r, p,          "TEM_Tests");
    Table.set("Significant",r, (!isNaN(p) && p < ALPHA) ? 1 : 0, "TEM_Tests");
    Table.update("TEM_Tests");
}
function addDesc(metric, group, n, mean, med, sd) {
    r = Table.size("TEM_Descriptives");
    Table.set("Metric", r, metric, "TEM_Descriptives");
    Table.set("Group",  r, group,  "TEM_Descriptives");
    Table.set("N",      r, n,      "TEM_Descriptives");
    Table.set("Mean",   r, mean,   "TEM_Descriptives");
    Table.set("Median", r, med,    "TEM_Descriptives");
    Table.set("SD",     r, sd,     "TEM_Descriptives");
    Table.update("TEM_Descriptives");
}
function splitTrim(s, delim) {
    parts = split(s, delim);
    out = newArray(0);
    for (i = 0; i < parts.length; i++) {
        t = String.trim(parts[i]);
        if (t != "") out = Array.concat(out, t);
    }
    return out;
}
function pad(s, w) { while (lengthOf(s) < w) s = s + " "; return s; }
function freshTable(name) {
    if (isOpen(name)) { selectWindow(name); run("Close"); }
    Table.create(name);
}
