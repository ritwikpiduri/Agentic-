#!/usr/bin/env bash
# =============================================================================
#  run_ilastik_batch.sh
#  Headless batch segmentation of a whole TEM dataset tree with a trained
#  ilastik Pixel Classification project. Writes one "*_seg.tif" Simple
#  Segmentation image per input TIFF, mirroring the folder structure, so
#  TEM_Measure_From_Ilastik.ijm can measure them.
#
#  PREREQUISITES
#   - You trained and saved organelles.ilp (see README_ilastik_pipeline.md).
#   - ilastik is installed.
#
#  EDIT THE TWO PATHS BELOW, then:  bash run_ilastik_batch.sh
# =============================================================================
set -euo pipefail

# ---- EDIT THESE ----
# Path to the ilastik headless launcher:
#   Linux : /opt/ilastik-<ver>-Linux/run_ilastik.sh
#   macOS : /Applications/ilastik-<ver>-OSX.app/Contents/ilastik-release/run_ilastik.sh
#   Windows (Git Bash): "/c/Program Files/ilastik-<ver>/ilastik.exe" --headless ...
ILASTIK="/opt/ilastik-1.4.0-Linux/run_ilastik.sh"
PROJECT="$(dirname "$0")/organelles.ilp"

# Root of your raw TEM dataset (the 4 folders live under here):
INPUT_ROOT="/path/to/TEM_dataset"
# Where segmentation images go (structure is mirrored under here):
OUTPUT_ROOT="/path/to/TEM_segmentation"
# --------------------

if [[ ! -x "$ILASTIK" ]]; then
  echo "ERROR: ilastik launcher not found/executable at: $ILASTIK" >&2
  echo "Edit ILASTIK at the top of this script." >&2
  exit 1
fi
if [[ ! -f "$PROJECT" ]]; then
  echo "ERROR: project not found: $PROJECT (train it first, see README)" >&2
  exit 1
fi

# Collect all TIFFs (handles spaces via -print0)
mapfile -d '' FILES < <(find "$INPUT_ROOT" -type f \( -iname '*.tif' -o -iname '*.tiff' \) -print0)
echo "Found ${#FILES[@]} TIFF images under $INPUT_ROOT"

count=0
for f in "${FILES[@]}"; do
  count=$((count+1))
  rel="${f#"$INPUT_ROOT"/}"                 # path relative to input root
  outdir="$OUTPUT_ROOT/$(dirname "$rel")"
  mkdir -p "$outdir"
  base="$(basename "${f%.*}")"
  outfile="$outdir/${base}_seg.tif"

  if [[ -f "$outfile" ]]; then
    echo "[$count/${#FILES[@]}] skip (exists): $rel"
    continue
  fi
  echo "[$count/${#FILES[@]}] segmenting: $rel"

  "$ILASTIK" --headless \
    --project="$PROJECT" \
    --export_source="Simple Segmentation" \
    --output_format=tif \
    --output_filename_format="$outfile" \
    --raw_data="$f"
done

echo "Done. Segmentation images in: $OUTPUT_ROOT"
echo "Now run TEM_Measure_From_Ilastik.ijm in Fiji on that folder."
