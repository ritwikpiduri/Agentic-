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
- `tem_stats.py` — aggregates the CSVs and runs group comparisons
  (Mann-Whitney / Kruskal-Wallis / chi-square).

## Step 0 — Calibrate (do NOT skip)

Every area/width/distance is meaningless without the correct **pixel size**
(µm/pixel). Find it from:
- the TIFF metadata (FEI/Gatan TEMs embed it), or
- a scale bar: draw a line along it, note its pixel length, divide the known
  length by that.

Set `PIXEL_SIZE_UM` in `TEM_Batch_Auto.ijm`. In the interactive macro use
action **[1]**. **If magnifications are mixed, process each magnification
separately** — one pixel size cannot be right for all of them.

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

**C. Statistics:**
```bash
pip install pandas scipy
python tem_stats.py --dir /path/to/csv_folder
```
Edit `assign_group()` in `tem_stats.py` so it reads your condition label from
the path (default: top-level folder name). Output flags `*SIGNIFICANT*` at
p<0.05 — with this many metrics, apply a Benjamini-Hochberg correction before
reporting.

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
