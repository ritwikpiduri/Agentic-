# Automatic organelle segmentation with ilastik → Fiji

Goal: stop tracing organelles by hand. You label a *handful* of images once;
ilastik learns the ultrastructure and segments **all 972 images**
automatically. Fiji then measures the result (counts, area, distance to
nucleus, near/far) into the same CSVs `tem_stats.py` reads.

Pipeline: **train (once) → batch-segment (unattended) → measure in Fiji →
stats**.

---

## 0. Install

Download ilastik (free, Win/Mac/Linux): https://www.ilastik.org/download.html
Open it once so you know where it installed — the batch script needs that path.

## 1. Create a Pixel Classification project

1. ilastik → **Pixel Classification** → save project as `organelles.ilp`
   (keep it next to this README).
2. **Add your images** for training only — pick **8–15 representative images**
   spanning your conditions and magnifications (more variety > more images).

## 2. Define classes — ORDER MATTERS

Add labels in **exactly this order** (ilastik numbers them 1,2,3… in the order
you add them, and the Fiji macro assumes this order):

| # | Label | Paint on… |
|---|-------|-----------|
| 1 | Background | cytoplasm/empty space, resin |
| 2 | Nucleus | nucleoplasm inside the nuclear envelope |
| 3 | ER | ribbon/cisternae of endoplasmic reticulum |
| 4 | Mitochondria | mitochondria (cristae texture) |
| 5 | Golgi | Golgi stacks |
| 6 | Vacuole | vacuoles / large empty vesicles |
| 7 | LipidBody | lipid droplets (uniform dense/round) |

> If you don't have a class (e.g. no Golgi), just omit it — but then also
> remove it from `CLASS_VALUES`/`CLASS_NAMES` in `TEM_Measure_From_Ilastik.ijm`
> so the numbering still lines up. The macro needs the class called **Nucleus**
> for distance measurements.

## 3. Choose features

**Features** tab → select a broad set across scales (σ ≈ 0.7–10 px):
Gaussian smoothing, Laplacian, Gaussian gradient magnitude, Difference of
Gaussians, **Structure Tensor Eigenvalues**, and **Hessian of Gaussian
Eigenvalues**. The Hessian/structure-tensor features capture membrane texture
that separates mito/ER/Golgi — include the larger scales.

## 4. Train interactively

1. **Training** tab. Turn on **Live Update**.
2. Brush a few strokes per class on clear examples. Watch the prediction
   overlay. **Correct mistakes** by painting the right class where it's wrong —
   this is the fastest way to improve it.
3. Iterate until borders look right on several images. You do **not** need
   pixel-perfect; the Fiji step filters specks and measures whole objects.
4. Save the project.

Tips for TEM specifically:
- Label plenty of **Background** and **membrane edges** — most errors are
  organelle-vs-background confusion at boundaries.
- If two organelles get confused, add corrective strokes on the confused
  regions rather than more strokes on easy areas.
- Keep magnification consistent within a comparison; if mixed, train on both.

## 5. Batch-segment everything (unattended)

Two options:

**A. GUI:** *Batch Processing* applet → add all 972 images → set
**Export Source = Simple Segmentation** → **Export**. In *Export Settings*
choose **Format = TIFF** and a filename pattern that preserves subfolders.

**B. Headless (recommended for 972 files):** use `run_ilastik_batch.sh`
(edit the two paths at the top), which runs ilastik from the command line over
your whole tree and writes one `*_seg.tif` per image. See that script's header.

"Simple Segmentation" output = a label image where each pixel value is the
class index (1=Background, 2=Nucleus, …). That's exactly what the Fiji macro
expects.

## 6. Measure in Fiji

1. Fiji → `Plugins ▸ Macros ▸ Run…` → `../TEM_Measure_From_Ilastik.ijm`.
2. Set `PIXEL_SIZE_UM` and confirm `CLASS_VALUES`/`CLASS_NAMES` match your
   training order.
3. Point it at the folder of segmentation TIFFs and an output folder.
4. It writes `Organelles.csv`, `Nuclei.csv`, `ImageSummary.csv` — automatically,
   no tracing — with per-object area, circularity, distance to nucleus, and
   near/far, plus nuclear shape and multinucleation.

## 7. Statistics (in Fiji, no Python)

`Plugins ▸ Macros ▸ Run…` → `../TEM_Stats.ijm`, point it at the CSV folder, and
type your condition keywords (e.g. `Control,Treated`). It writes a Log report
and `TEM_Descriptives.csv` / `TEM_Tests.csv`.

## Validate before you trust it

Auto-segmentation is powerful but not infallible. Open 10–20 QC pairs (raw vs
segmentation) across conditions. If a class is systematically wrong, add
training strokes on those cases and re-run steps 4–6 — that loop is much faster
than tracing 972 images by hand.

## What still isn't automatic

- **Perinuclear space width** and **micronuclei** stay in the interactive macro
  (`../TEM_Analysis_Interactive.ijm`) — they need judgement at high mag. The
  ImageSummary micronucleus count from this pipeline is 0 by design; fill it
  from the interactive pass if you need it.
- Ultrastructure *identity* calls should still be confirmed by an EM
  microscopist.
