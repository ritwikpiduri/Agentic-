# LabPack — Nuclear TEM measurements (Fiji macros)

Automated nuclear morphometry for TEM images, for **any microscope / cell
model**. You trace and click through a menu; the macros do the measuring,
calibration math and bookkeeping, and save a CSV after every entry.

This is the **nuclear** workflow only (shape, chromatin, nuclear envelope,
micronuclei, multinucleation). Organelles are not included here.

## Files

| File | What it does |
|---|---|
| `LabPack_Nuclear.ijm` | Measures each image via a menu → `LabPack_Nuclear.csv` + QC pictures |
| `LabPack_Summary.ijm` | Reads the CSV(s) and prints per-condition means → `LabPack_Summary.csv` |

## Before you start

1. Install **Fiji**: https://imagej.net/software/fiji/ (unzip and run).
2. **One folder per condition**, e.g. `Control/`, `Treatment_A/`, … The folder
   name becomes the `Sample` label in the output, so keep it clean.
3. Know how your images are scaled — you have three options (chosen in the first
   dialog):
   - **Scale bar** (recommended): you draw a line on the burned-in scale bar and
     type its length. Works even if the image has no embedded calibration.
   - **Fixed pixel size**: if all images are the same magnification and you know
     the µm/pixel, type it once.
   - **Embedded**: if your microscope already saved the scale into the TIFF.

## Running `LabPack_Nuclear.ijm`

1. Fiji ▸ **Plugins ▸ Macros ▸ Run…** ▸ pick `LabPack_Nuclear.ijm`
   (or drag it into the Fiji window and press **Run**).
2. Choose the **image folder** for one condition, then an **output folder**.
3. In **Settings**, pick your **calibration mode** and the **heterochromatin
   grey cutoff** (see the chromatin note below).
4. For each image a menu appears:
   - **Set scale** — only in scale-bar mode: draw a straight line along the
     scale bar, type its length + unit. The menu header shows
     `calibrated: … um/px` once it's set.
   - **Nucleus** — trace around the nuclear envelope (freehand). You get area,
     perimeter, circularity, aspect ratio, roundness, solidity, and
     heterochromatin %. Trace a **2nd/3rd nucleus in the same image** → that
     image is scored as **multinucleated**.
   - **Micronuclei** — trace each small body separate from the nucleus
     (area + circularity).
   - **Perinuclear space** — draw short lines across the gap between the inner
     and outer membranes → nm (take several per nucleus; they're averaged later).
   - **NEXT / SKIP / QUIT**.
5. Output is saved continuously to `LabPack_Nuclear.csv`. A QC picture of every
   traced nucleus (yellow outline) is saved in `QC_nucleus/` so anyone can check
   what was measured. If you stop, rerun and choose **RESUME**, then start at the
   image number shown in the title bar.

> Run **one condition folder at a time** (or point it at a parent folder — each
> image's `Sample` comes from its own folder name, so conditions stay separate).

## Getting per-condition results

Run **`LabPack_Summary.ijm`**, point it at the folder holding your CSV(s). It
prints and saves `LabPack_Summary.csv` with, per condition: nuclei measured;
mean area, circularity, aspect ratio, solidity; mean heterochromatin / euchromatin
% and mean grey; % multinucleated; micronuclei per nucleus; mean perinuclear
space (nm). Do your statistics (ANOVA / Kruskal–Wallis for >2 groups) in
Excel/Prism on the per-nucleus values; **n = number of nuclei**, not images.

## What each measurement means (quick reference)

| Column | Meaning | Unit |
|---|---|---|
| `Area_um2` | nuclear size | µm² |
| `Perimeter_um` | boundary length | µm |
| `Circularity` | closeness to a circle; drops with ruffling/folds (4π·Area/Perim²) | 0–1 |
| `AspectRatio` | elongation (long axis ÷ short axis of fitted ellipse) | ≥1 |
| `Roundness` | 1/AspectRatio (redundant — report one) | 0–1 |
| `Solidity` | Area ÷ convex-hull area; drops with lobes/indentations | 0–1 |
| `Heterochrom_pct` | % of nuclear area darker than the fixed cutoff | % |
| `Euchrom_pct` | 100 − heterochromatin | % |
| `MeanGrey` | mean grey of the nucleus (cutoff-free backup; lower = darker) | 0–255 |
| `PerinuclearSpace_nm` | gap between inner and outer nuclear membranes | nm |

## ⚠️ Chromatin: use ONE fixed cutoff, and keep imaging constant

Heterochromatin (dense chromatin) is **dark** on TEM. The macro counts the
fraction of the nucleus **darker than a fixed grey cutoff** (e.g. 100), set once
and used for every image.

**Do not** use a per-nucleus / automatic threshold: it re-centres on each nucleus
and makes every condition look roughly equal — an artifact. For the comparison to
be valid, **staining and exposure must be identical across all conditions**, and
the **same cutoff** used throughout. `MeanGrey` is saved as a cutoff-free backup;
if two methods disagree, trust the mean grey.

To choose a cutoff: open a few representative images, trace a nucleus, and use
**Analyze ▸ Histogram** to see where the dark chromatin sits; pick a value that
separates dense from loose chromatin, then keep it fixed.

## Calibration reminders

- **Calibrate per magnification.** If your images are at different mags, set the
  scale separately for each (the macro lets you re-set the scale per image).
- **Don't reuse another lab's pixel-size table** — pixel size is specific to your
  detector + magnification. The scale bar is ground truth.
- **Sanity-check**: after calibrating, measure something of known size
  (a mitochondrion ≈ 0.5–1 µm, or the scale bar itself) and confirm it's sensible.

## Common mistakes

- Tracing the **whole cell** instead of the **nucleus** (they can look similar in
  low-contrast TEM) — follow the nuclear envelope.
- Calling the perinuclear gap "membrane thickness" — it's the **space between two
  membranes** (nm), not a single bilayer (~7–8 nm).
- Counting **images or cell lines as n** instead of nuclei (pseudoreplication).
- Changing staining/exposure between conditions, then comparing chromatin.
