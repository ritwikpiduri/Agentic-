#!/usr/bin/env python3
"""
tem_stats.py - aggregate the CSVs produced by the Fiji TEM macros and test
whether metrics differ between samples/groups.

Handles the questions:
  * nuclear circularity / area change between samples?
  * more heterochromatin or euchromatin, and does it shift between samples?
  * multinucleated-cell frequency, micronucleus frequency change?
  * organelle counts (ER, mito, Golgi, vacuole, lipid) up or down?
  * organelle proximity to nucleus (near vs far) shift?
  * perinuclear space width change?

USAGE
  python tem_stats.py --dir /path/to/csv_folder --group-col Sample
  # or map each folder/file to a group explicitly, see --group-map

Each input CSV has an "Image" column. A "group" (sample/condition) is derived
from that image name or its Folder column. Edit `assign_group()` to match how
YOUR files/folders encode the condition (e.g. "Control", "Treated_24h").

Comparisons: 2 groups -> Mann-Whitney U (non-parametric, safe for EM counts);
>2 groups -> Kruskal-Wallis. Proportions (multinucleation, near/far) -> chi-sq.
Requires: pandas, scipy.  pip install pandas scipy
"""
import argparse
import os
import sys
import pandas as pd
from scipy import stats


def assign_group(name: str) -> str:
    """Map an image name / folder path to a sample/condition label.

    EDIT THIS to your naming. Example assumes the top-level folder name is the
    condition, e.g. 'Control/cell03/img012.tif' -> 'Control'.
    """
    name = str(name).replace("\\", "/")
    parts = [p for p in name.split("/") if p]
    return parts[0] if parts else "all"


def load(csv_dir: str, fname: str) -> pd.DataFrame:
    p = os.path.join(csv_dir, fname)
    if not os.path.exists(p):
        return pd.DataFrame()
    df = pd.read_csv(p)
    key = "Folder" if "Folder" in df.columns else "Image"
    df["group"] = df[key].map(assign_group)
    return df


def compare_numeric(df: pd.DataFrame, value_col: str, label: str):
    if df.empty or value_col not in df.columns:
        return
    groups = [g[value_col].dropna().values for _, g in df.groupby("group")]
    names = list(df.groupby("group").groups.keys())
    if len([g for g in groups if len(g)]) < 2:
        return
    print(f"\n== {label} ({value_col}) ==")
    for n, g in zip(names, groups):
        if len(g):
            print(f"  {n:20s} n={len(g):4d}  "
                  f"mean={g.mean():.4f}  median={pd.Series(g).median():.4f}  "
                  f"sd={g.std(ddof=1) if len(g)>1 else 0:.4f}")
    try:
        if len(groups) == 2:
            u, p = stats.mannwhitneyu(groups[0], groups[1], alternative="two-sided")
            print(f"  Mann-Whitney U={u:.1f}  p={p:.4g}  "
                  f"{'*SIGNIFICANT*' if p < 0.05 else '(n.s.)'}")
        else:
            h, p = stats.kruskal(*[g for g in groups if len(g)])
            print(f"  Kruskal-Wallis H={h:.2f}  p={p:.4g}  "
                  f"{'*SIGNIFICANT*' if p < 0.05 else '(n.s.)'}")
    except ValueError as e:
        print(f"  (test skipped: {e})")


def compare_proportion(df: pd.DataFrame, flag_col: str, label: str):
    if df.empty or flag_col not in df.columns:
        return
    tab = df.groupby("group")[flag_col].agg(["sum", "count"])
    tab["fraction"] = tab["sum"] / tab["count"]
    print(f"\n== {label} (proportion) ==")
    print(tab.to_string())
    if len(tab) >= 2:
        contingency = [[r["sum"], r["count"] - r["sum"]] for _, r in tab.iterrows()]
        try:
            chi2, p, _, _ = stats.chi2_contingency(contingency)
            print(f"  chi-square={chi2:.2f}  p={p:.4g}  "
                  f"{'*SIGNIFICANT*' if p < 0.05 else '(n.s.)'}")
        except ValueError as e:
            print(f"  (test skipped: {e})")


def organelle_counts_per_image(summary: pd.DataFrame):
    cols = ["ER_count", "Mito_count", "Golgi_count", "Vacuole_count", "Lipid_count"]
    for c in cols:
        compare_numeric(summary, c, f"Organelle count per image")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="folder containing the CSVs")
    args = ap.parse_args()
    d = args.dir

    nuclei      = load(d, "Nuclei.csv")
    chromatin   = load(d, "Chromatin.csv")
    perinuclear = load(d, "Perinuclear.csv")
    organelles  = load(d, "Organelles.csv")
    micronuclei = load(d, "Micronuclei.csv")
    summary     = load(d, "ImageSummary.csv")
    batch       = load(d, "TEM_Batch_NucleusChromatin.csv")

    # Nuclear shape
    compare_numeric(nuclei, "Circularity", "Nuclear circularity")
    compare_numeric(nuclei, "Area_um2",    "Nuclear area")
    compare_numeric(nuclei, "Solidity",    "Nuclear solidity (envelope irregularity)")

    # Batch shape/chromatin (if you ran the auto pass)
    compare_numeric(batch, "Circularity",     "Nuclear circularity (batch)")
    compare_numeric(batch, "Heterochrom_pct", "Heterochromatin % (batch)")

    # Chromatin (interactive)
    compare_numeric(chromatin, "Heterochrom_pct", "Heterochromatin %")
    compare_numeric(chromatin, "Hetero_Eu_ratio", "Hetero/Eu ratio")

    # Perinuclear space
    compare_numeric(perinuclear, "Width_um", "Perinuclear space width")

    # Multinucleation & micronuclei frequency (per image)
    compare_proportion(summary, "Multinucleated", "Multinucleated cells")
    compare_numeric(summary, "MicronucleusCount", "Micronuclei per image")

    # Organelle counts and proximity
    organelle_counts_per_image(summary)
    if not organelles.empty and "NearFar" in organelles.columns:
        organelles["is_near"] = (organelles["NearFar"] == "near").astype(int)
        compare_proportion(organelles, "is_near", "Organelles near nucleus (<=cutoff)")
        compare_numeric(organelles, "Dist_to_nucleus_um", "Organelle distance to nucleus")

    print("\nDone. '*SIGNIFICANT*' means p<0.05 (uncorrected). With many metrics,"
          "\napply a multiple-comparison correction (e.g. Benjamini-Hochberg).")


if __name__ == "__main__":
    main()
