# TEM Nuclear Morphometry Protocol (Fiji / ImageJ)

A step-by-step protocol for measuring **nuclear shape, chromatin, nuclear
envelope, micronuclei and multinucleation** from transmission electron
microscopy (TEM) images, using **Fiji (ImageJ)** only — no coding or Python
required. It is written so it can be applied to **any cell model**: you supply
your own images and your own scale calibration.

> Organelle analysis (ER / mitochondria / Golgi / vacuoles / lipid bodies /
> lysosomes, near-vs-far positioning, counts) is intentionally **not** covered
> here — this document is the nuclear workflow only.

---

## 0. What you will measure and why

| Measurement | What it tells you | Unit |
|---|---|---|
| **Area** | Nuclear size | µm² |
| **Perimeter** | Length of the nuclear boundary | µm |
| **Circularity** | How close the outline is to a perfect circle (sensitive to boundary ruffling/folds) | none (0–1) |
| **Aspect ratio (AR)** | How elongated the nucleus is (long axis ÷ short axis) | none (≥1) |
| **Roundness** | Elongation, expressed as 1/AR | none (0–1) |
| **Solidity** | How "dented"/lobed the outline is (area ÷ its convex hull) | none (0–1) |
| **Heterochromatin %** | Fraction of the nucleus occupied by dense (dark) chromatin | % of area |
| **Euchromatin %** | Fraction occupied by light chromatin (= 100 − heterochromatin) | % of area |
| **Perinuclear space** | Gap between the inner and outer nuclear membranes | nm |
| **Micronuclei** | Small chromatin bodies separate from the main nucleus (genomic instability) | count / nucleus |
| **Multinucleation** | Cells with ≥2 nuclei | % of cells |

**Key principle for all of these:** a measurement is only comparable between
samples if **everything upstream is held constant** — fixation, sectioning,
staining, microscope, magnification, and (critically) **calibration**. Decide
your imaging settings once and keep them identical across all conditions.

---

## 1. Software and one-time setup

1. Download **Fiji** (ImageJ bundled with plugins): https://imagej.net/software/fiji/
   Unzip and run; nothing to install.
2. Set the measurements you want **once per session**:
   **Analyze ▸ Set Measurements…** and tick:
   - **Area**
   - **Perimeter**
   - **Shape descriptors** (gives Circularity, Aspect Ratio, Roundness, Solidity)
   - **Fit ellipse** (gives Major/Minor axis — the basis of AR)
   - **Mean gray value** (used for chromatin)
   - Set **Decimal places = 3–4**.
   - Leave *Redirect to* = None.
3. Useful tools on the toolbar:
   - **Freehand selection** (for tracing the nucleus)
   - **Straight line** (for perinuclear space / membrane gaps)
   - **Multi-point** (for counting pores, micronuclei)
   - **Wand** is *not* reliable on TEM (the cell is often darker than
     background and the wand floods the whole cell) — trace by hand.

---

## 2. Organise your images (do this before measuring)

```
TEM_data/
├── Control/        <- all images for condition 1
├── Treatment_A/    <- condition 2
├── Treatment_B/    <- condition 3
└── Treatment_C/    <- condition 4
```

- **One folder per condition.** The folder name becomes the sample label in
  your spreadsheet, so keep them clean (no spaces if you can avoid it).
- Record the **magnification** of every image (it's usually in the filename or
  the microscope metadata). You need it for calibration.
- Low-magnification images (whole nucleus visible) are used for **shape,
  chromatin, micronuclei, multinucleation**. High-magnification images are used
  for **perinuclear space / membrane** detail.
- Aim for **n ≈ 20 nuclei per condition** (more is better). Count **nuclei as
  your n**, not cell lines or images — otherwise you pseudoreplicate.

---

## 3. Calibration — the single most important step

Every µm and nm value depends on telling Fiji **how big one pixel is**. If this
is wrong, *all* size measurements are wrong (shape descriptors like circularity
are unitless and survive, but Area, Perimeter, perinuclear space do not).

You have two ways to calibrate. Use whichever your images support.

### Option A — from the scale bar burned into the image (most robust)

1. Open an image. Select the **Straight line** tool.
2. Draw a line **exactly along the scale bar**, matching its length.
3. **Analyze ▸ Set Scale…**
   - *Distance in pixels* is filled in automatically from your line.
   - *Known distance* = the number printed on the scale bar (e.g. 2).
   - *Unit of length* = the scale bar's unit (e.g. `um` for µm, or `nm`).
   - Tick **Global** so it applies to every image you open this session
     **only if all your images are at the same magnification**. If
     magnifications differ, do **not** use Global — calibrate each
     magnification separately.
4. Fiji now reports measurements in your chosen unit.

### Option B — from a known pixel size (if your microscope software gives it)

Many TEM cameras report a fixed **pixel size at each magnification** (µm/pixel).
If you have that table:

1. **Image ▸ Properties…**
2. Set **Pixel width** and **Pixel height** to the µm/pixel value for that
   image's magnification, **Unit = micron**. Leave Voxel depth = 1.

> ⚠️ Do **not** reuse another lab's pixel-size table or camera constant unless
> it is **your** microscope and camera. Pixel size is specific to the
> detector + magnification. When in doubt, use **Option A** (the scale bar is
> ground truth).

**Sanity check:** measure a structure of known size (e.g. a mitochondrion
≈ 0.5–1 µm, or the scale bar itself) and confirm the number is sensible before
trusting a whole dataset.

---

## 4. Measure nuclear shape (Area, Perimeter, Circularity, AR, Roundness, Solidity)

These all come from **one traced outline** of the nucleus.

### Steps
1. Open a **low-mag** image where the **whole nucleus** is visible.
2. Confirm it's calibrated (Section 3). The scale bar should read in µm.
3. Select the **Freehand selection** tool.
4. **Trace around the nuclear envelope** — follow the outer edge of the nucleus
   all the way around, back to the start. Take your time; the outline is the
   measurement. Zoom in (+) for ruffled edges.
5. Press **Ctrl+M** (Analyze ▸ Measure). A Results row appears.
6. Repeat for each nucleus. **Save the Results table** as CSV when done
   (Results window ▸ File ▸ Save As).

### What each number means (so you can explain it)

- **Area (µm²)** — the amount of 2-D space the nucleus covers in this section.
  Bigger = larger nucleus (or a section cut through its widest part).

- **Perimeter (µm)** — the length of the boundary you traced. A folded/ruffled
  envelope has a longer perimeter for the same area.

- **Circularity** = 4π × Area ÷ Perimeter².
  - Range **0 to 1**. **1 = a perfect circle.**
  - Drops when the outline is **irregular or folded** (perimeter grows faster
    than area). Very sensitive to small ruffles.
  - *Manual intuition:* a smooth round nucleus ≈ 0.9+; a crinkled/indented one
    is lower.
  - **Unitless** (it's a ratio of an area to a length², so the units cancel).

- **Aspect Ratio (AR)** = Major axis ÷ Minor axis of the **best-fit ellipse**.
  - **≥ 1.** AR = 1 means the fitted ellipse is a circle (not elongated);
    AR = 2 means twice as long as it is wide.
  - Captures **elongation/stretching**, independent of ruffling.
  - *Manual intuition:* fit an oval to the nucleus; divide its long diameter by
    its short diameter.
  - **Unitless** (a length ÷ a length).

- **Roundness** = 4 × Area ÷ (π × Major axis²) = **1 ÷ AR**.
  - Range 0–1; 1 = not elongated. It is just AR expressed the other way round —
    **redundant with AR**, so report one of them, not both. (Roundness ignores
    ruffling, unlike circularity.)
  - **Unitless.**

- **Solidity** = Area ÷ **Convex-hull area**.
  - The **convex hull** is the smallest convex (no-dents) shape that wraps the
    outline — imagine stretching a rubber band around the nucleus. Solidity
    compares the real area to that rubber-band area.
  - Range **0 to 1**. **1 = perfectly convex** (no concavities). Lower values
    mean **deep indentations / lobes / invaginations** — the nucleus has
    "bites" taken out of its outline.
  - *Difference from circularity:* circularity reacts to **fine edge
    roughness**; solidity reacts to **large-scale lobes/indentations**. A
    nucleus can be smooth-edged but kidney-shaped (good circularity-ish, lower
    solidity), or round overall but finely ruffled (good solidity, lower
    circularity). Reporting both is informative.
  - **Unitless.**
  - *Literature for context:* nuclear irregularity / solidity as a morphometric
    marker of nuclear envelope abnormality is used in cancer cytopathology and
    laminopathy work — cite a solidity/nuclear-contour paper from your field so
    readers know what "low solidity" is pointing at.

**What to plot:** Area, Circularity, and **either** AR **or** Solidity (they
capture different things — elongation vs. lobulation). Don't plot both AR and
Roundness (same information).

---

## 5. Measure chromatin (heterochromatin vs euchromatin %)

**Biology:** heterochromatin is **dense and dark** on TEM; euchromatin is
**loose and light**. So "% heterochromatin" = the fraction of the nuclear area
whose pixels are **darker** than a chosen cutoff.

> ⚠️ The big pitfall: if you let Fiji pick the threshold **per nucleus**
> (auto / relative methods), it re-centres on each nucleus and you get roughly
> the **same %** everywhere — an artifact. Use a **single fixed grey cutoff**
> for the whole dataset, and only compare images taken with **identical
> staining and exposure**.

### Steps
1. Trace the nucleus outline (Section 4) — or reuse the saved outline.
2. **Edit ▸ Selection ▸ Restore Selection** if you need to bring the outline
   back, then **Image ▸ Crop** is *not* needed — work inside the selection.
3. Make sure the image is **8-bit** (Image ▸ Type ▸ 8-bit) so grey values run
   0 (black) – 255 (white).
4. **Analyze ▸ Histogram** on the selection to see the chromatin's grey
   distribution. Pick a **fixed cutoff** (e.g. pixels ≤ 100 = heterochromatin)
   and **use that same number for every image**.
5. To get the percentage: with the nucleus selected,
   **Image ▸ Adjust ▸ Threshold…**, set the lower/upper to highlight the dark
   chromatin (e.g. 0–100), **don't Apply** — instead read the thresholded
   fraction, or use **Analyze ▸ Measure** with *Limit to threshold* ticked in
   Set Measurements to get the dark **Area**, then:

   **Heterochromatin % = (dark area ÷ total nuclear area) × 100**
   **Euchromatin % = 100 − Heterochromatin %**

6. Also record the **Mean gray value** of the whole nucleus as a
   **cutoff-free backup** (lower mean = darker = more heterochromatin overall).
   If two methods disagree, the mean grey is the honest fallback.

**Interpretation tip (lamin / chromatin studies):** total heterochromatin %
can look flat across samples even when molecular assays (ATAC, MNase) say
otherwise — because TEM sees **where** chromatin sits, not just how much.
If your biology predicts redistribution, also measure **peripheral vs interior**
heterochromatin (split the nucleus into an outer rim and a core and compare the
dark fraction in each) rather than only the global %.

---

## 6. Measure the nuclear envelope (perinuclear space)

**Definition:** the **perinuclear space** is the **gap between the inner and
outer nuclear membranes** (the two parallel lines of the nuclear envelope).
Normal is ~**30–50 nm**. This is **not** "membrane thickness" — a single lipid
bilayer is only ~**7–8 nm**; you are measuring the **space between the two
membranes**, which is why nm is the right unit.

### Steps
1. Open a **high-mag** image where both membranes of the envelope are clearly
   resolved as two lines.
2. Calibrate in **nm** (Section 3) — at high mag your scale bar is usually in nm.
3. Select the **Straight line** tool.
4. Draw a short line **perpendicular to the envelope**, from the **outer edge of
   the inner membrane to the inner edge of the outer membrane** (i.e. spanning
   the gap). Press **Ctrl+M**; read **Length**.
5. Take **several measurements per nucleus** (e.g. 3–5 positions around the
   envelope) and **average** them — the gap varies around the nucleus.
6. Report the mean perinuclear space (nm) per nucleus, then per sample.

**On your figure axis**, label it **"Perinuclear space (nm)"** — clearer to
readers than "membrane thickness," and it's what you actually measured.

---

## 7. Count micronuclei and score multinucleation

### Micronuclei
- A **micronucleus** is a **small, round, membrane-bound chromatin body**,
  **separate** from the main nucleus, usually nearby, with similar density to
  the nucleus. It reflects chromosome mis-segregation / genomic instability.
- **Not** the same as **multinucleation** (which is ≥2 **full-size** nuclei in
  one cell).
- **Count** them with the **Multi-point** tool (one click each) or just tally by
  eye, and record **micronuclei ÷ number of nuclei** per sample.

### Multinucleation
- Score **per cell**: does this cell have **≥2 nuclei**? Record the fraction of
  imaged cells that are multinucleated:
  **% multinucleated = (cells with ≥2 nuclei ÷ total cells scored) × 100**.
- This needs **whole-cell** (low-mag) fields so you can see the whole cell.

---

## 8. Record your data consistently

Make one spreadsheet with **one row per nucleus** and columns:

```
Sample | Image | Area_um2 | Perimeter_um | Circularity | AspectRatio |
Solidity | Heterochrom_pct | Euchrom_pct | PerinuclearSpace_nm |
Micronuclei | Multinucleated(Y/N)
```

- Save Fiji's Results table as **CSV** frequently (you can't undo a lost table).
- Keep the **units in the column name** so nobody has to guess.
- Blinding: if possible, have the folders renamed by someone else so you measure
  without knowing the condition — removes unconscious bias.

---

## 9. Analysis and statistics (brief)

- Compare **per-nucleus** values across conditions. Your **n = number of
  nuclei**.
- With **4 groups**, don't run a pile of t-tests. Use **one-way ANOVA**
  (if roughly normal) or **Kruskal–Wallis** (if not), then the appropriate
  post-hoc test. For a single planned 2-group comparison, a **two-tailed
  t-test** (Welch's, unequal variance) is fine.
- Plot **mean ± SEM** (or better, show the individual points) with the sample
  size stated.
- **Sanity-check surprising p-values:** if two groups have nearly identical
  means but a tiny p-value, something is off (often an outlier or a
  units/calibration mistake) — go back and look at the raw ranges.

---

## 10. Common mistakes to avoid

- **Mismatched magnifications** sharing one Global calibration → all sizes
  wrong. Calibrate per magnification.
- **Auto/relative chromatin threshold** → everything looks equal (artifact).
  Use one fixed cutoff + consistent imaging; keep mean grey as backup.
- **Tracing the whole cell instead of the nucleus** (the two can look similar in
  low-contrast TEM). Trace the **nuclear envelope**.
- **Calling the perinuclear gap "membrane thickness."** It's the space between
  two membranes (nm), not a bilayer.
- **Counting lines/images as n** instead of nuclei (pseudoreplication).
- **Changing staining or exposure between conditions** — then chromatin
  intensity comparisons are meaningless.

---

*Prepared for lab use. Adapt the thresholds, magnifications and n to your own
cell model and microscope; the calibration step (Section 3) is what makes it
portable.*
