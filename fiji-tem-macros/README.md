# TEM image analysis toolkit (Fiji/ImageJ)

Quantifies nuclear morphology, chromatin, perinuclear space, organelles and
micronuclei from TEM images, then tests differences between samples.

Dataset in mind: **4 folders × ~20 cells × 243 TIFFs**. Organise as
`ROOT/<Condition>/<cellNN>/<image>.tif` so grouping and stats work
automatically.

## What each metric needs — read this first

TEM is single-channel grayscale, so some things automate cleanly and some
genuinely need your eye. Be realistic about which is which:

| Metric | How | Automatable? |
|---|---|---|
| Nuclear circularity, area, solidity | threshold the nucleus, shape descriptors | **Yes** (`TEM_Batch_Auto.ijm`), spot-check |
| Heterochromatin vs euchromatin | dark(dense) vs light pixel fraction inside nucleus | **Yes**, spot-check the threshold |
| Multinucleation | count nuclei per cell/image | Semi — confirm by eye |
| Micronuclei | small chromatin body separate from nucleus | Manual trace (validation needed) |
| Perinuclear space width | line across inner→outer nuclear membrane | Manual, high-mag images only |
| ER / Golgi / mito / vacuole / lipid | identify by ultrastructure | **Manual** — no reliable auto signature |
| Organelle distance to nucleus (near/far) | centroid vs nuclear-edge distance map | Auto **once you mark the organelle** |

> Organelles have no universal intensity signature in TEM, so nothing here
> auto-detects them. You mark them; the macro measures, counts, and computes
> distances. If you want true auto-segmentation, train **Trainable Weka
> Segmentation**, **Labkit**, or **ilastik** on a handful of images and feed the
> masks into the "measure from mask" section of `TEM_Batch_Auto.ijm`.

## Files

- `TEM_Analysis_Interactive.ijm` — menu-driven, run on the open image. Handles
  everything above; appends to per-metric tables you save as CSV.
- `TEM_Batch_Auto.ijm` — unattended pass over the whole tree for nuclear shape +
  chromatin. Writes `TEM_Batch_NucleusChromatin.csv` and QC overlay PNGs.
- `TEM_Stats.ijm` — **statistics inside Fiji, no Python needed.** Reads the
  CSVs and runs Mann-Whitney U (numeric) and chi-square (proportions) between
  your conditions, writing a Log report plus `TEM_Descriptives.csv` /
  `TEM_Tests.csv`.
- `tem_stats.py` — *optional* same analysis in Python (only if you have pandas +
  scipy). Ignore it if you only have ImageJ.
- `ilastik/` — **automatic organelle detection so you don't trace by hand.**
  Train a classifier once, batch-segment all images, then
  `TEM_Measure_From_Ilastik.ijm` measures counts/area/distance into the same
  CSVs. See `ilastik/README_ilastik_pipeline.md`. This is the recommended route
  for the organelle work at your scale (972 images).

## Step 0 — Calibration (mostly automatic here)

Every area/width/distance depends on the correct **pixel size** (µm/pixel).
**Check which case you're in:** open one TIFF → `Image ▸ Properties`.

- **Pixel width shows a real value + a unit (micron/nm)** → the scale is in the
  metadata. The macros **read it automatically, per image** (`USE_EMBEDDED_SCALE`
  is on), and normalise nm/Å to microns. Mixed magnifications just work, since
  each image uses its own scale. Nothing to set. `TEM_Batch_Auto.ijm` records a
  `ScaleSource` column and warns if any image lacked a scale.
- **Pixel width shows `1 pixel`** → the scale is only a *drawn scale bar*, which
  ImageJ can't read. Measure the bar once (line tool → its pixel length), set
  `PIXEL_SIZE_UM` in the macros as the fallback, or use interactive action
  **[1]**. If magnifications are mixed here, process each separately.

> The ilastik label images lose the original scale, so
> `TEM_Measure_From_Ilastik.ijm` asks for your **original TIFF root** and reads
> the real scale from each matching original automatically.

## Workflow

**A. Automated first pass (all images):**
1. Fiji → `Plugins ▸ Macros ▸ Run…` → `TEM_Batch_Auto.ijm`.
2. Pick the ROOT folder and an output folder.
3. Open the `QC_overlays/` PNGs, delete rows in the CSV where the nucleus
   outline is obviously wrong. This gives you circularity + heterochromatin %
   for the whole dataset fast.

**B. Interactive pass (the things needing your eye):**
1. Open an image → `Run…` → `TEM_Analysis_Interactive.ijm`.
2. `[1]` calibrate → `[2]` trace each nucleus (run once per nucleus; ≥2 ⇒
   multinucleated) → `[3]` chromatin → `[4]` perinuclear lines → `[5]`
   organelles (pick type, trace each, near/far auto) → `[6]` micronuclei →
   `[7]` write the image summary. Re-open the macro for the next image; tables
   persist across images in the same Fiji session.
3. `[9]` save all CSVs when done (or periodically — it re-writes the same files).

**C. Statistics (all in Fiji — no Python):**
1. `Plugins ▸ Macros ▸ Run…` → `TEM_Stats.ijm`.
2. Point it at your CSV folder and type your condition keywords, e.g.
   `Control,Treated` (each data row is assigned to the first keyword found in
   its Image/Folder/Sample text, so it works whether the condition is in the
   folder path or the filename).
3. Read the report in the Log window; `TEM_Descriptives.csv` and
   `TEM_Tests.csv` are saved next to your data. Numeric metrics use
   Mann-Whitney U, proportions (multinucleation, near/far) use chi-square.
   `*SIGNIFICANT*` = p<0.05; with many metrics apply a Benjamini-Hochberg
   correction before reporting.

> `tem_stats.py` does the same thing in Python but is optional — skip it if you
> only have ImageJ.

## Tips & caveats

- **Consistency drives validity**: same threshold logic, same magnification per
  comparison, blinded tracing where possible.
- **Heterochromatin = electron-dense = dark**. If your images are inverted,
  invert the comparison.
- **Distance map** gives distance from the nuclear *edge* in µm; 0 means the
  organelle centroid sits inside the traced nucleus.
- **Sample size**: 243 images/cell is a lot of manual organelle work — automate
  A across everything, then do B on a randomised, adequately-powered subset
  rather than all 243.
- These macros measure; they don't diagnose. Have ultrastructure calls
  (what *is* a vacuole vs autophagosome) confirmed by an EM microscopist.
