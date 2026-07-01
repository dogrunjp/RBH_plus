#!/usr/bin/env bash
set -euo pipefail

# Example wrapper for generating a human -> Heterocephalus glaber mapping table
# using RBH+in-paralog from DIAMOND all-vs-all hits.
#
# Usage:
#   ./src/run_human_to_hglaber.sh HUMAN_FASTA TARGET_FASTA [OUTDIR] [TARGET_ID_MAP]
#
# The script will:
#   1. Create ID lists from FASTA headers
#   2. Build a combined DIAMOND protein DB
#   3. Run all-vs-all DIAMOND blastp
#   4. Run rbh_plus_orthologs.py to emit a mapping table
#   5. Convert the mapping to PathLift's provided-table format
#
# Output:
#   <outdir>/map_H_to_T.tsv
#   <outdir>/pathlift_H_to_T.tsv

if [[ $# -lt 2 || $# -gt 4 ]]; then
  echo "Usage: $0 HUMAN_FASTA TARGET_FASTA [OUTDIR] [TARGET_ID_MAP]" >&2
  echo "Example: $0 /path/to/human.faa /path/to/hglaber.faa ./output/human_hglaber" >&2
  exit 1
fi

HUMAN_FASTA=$1
TARGET_FASTA=$2
OUTDIR=${3:-./output/human_hglaber}
TARGET_ID_MAP=${4:-}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  PYTHON_BIN="$CONDA_PREFIX/bin/python"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN=python
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN=python3
else
  echo "Error: python interpreter not found" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

HUMAN_IDS="$OUTDIR/human.ids"
TARGET_IDS="$OUTDIR/hglaber.ids"
COMBINED_FASTA="$OUTDIR/human_hglaber.faa"
DB_PATH="$OUTDIR/human_hglaber"
HITS_TSV="$OUTDIR/human_hglaber_hits.tsv"
RBH_MAPPING="$OUTDIR/map_H_to_T.tsv"
PATHLIFT_MAPPING="$OUTDIR/pathlift_H_to_T.tsv"

extract_ids_from_fasta() {
  local fasta=$1
  local out=$2
  awk '
    /^>/ {
      id=$0
      sub(/^>/, "", id)
      sub(/[[:space:]].*$/, "", id)
      if (id != "") print id
    }
  ' "$fasta" > "$out"
}

extract_ids_from_fasta "$HUMAN_FASTA" "$HUMAN_IDS"
extract_ids_from_fasta "$TARGET_FASTA" "$TARGET_IDS"

cat "$HUMAN_FASTA" "$TARGET_FASTA" > "$COMBINED_FASTA"

echo "[1/5] Building DIAMOND database"
diamond makedb --in "$COMBINED_FASTA" -d "$DB_PATH"

echo "[2/5] Running all-vs-all DIAMOND blastp"
diamond blastp \
  -d "$DB_PATH" \
  -q "$COMBINED_FASTA" \
  --evalue 1e-5 --max-target-seqs 100 --threads 16 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen qstart qend sstart send \
  -o "$HITS_TSV"

echo "[3/5] Running RBH+in-paralog analysis"
"$PYTHON_BIN" "$REPO_ROOT/src/rbh_plus_orthologs.py" \
  --hits "$HITS_TSV" \
  --species H "$HUMAN_IDS" \
  --species T "$TARGET_IDS" \
  --pairs H:T \
  --tau 0.7 --min_cov 0.6 --len_ratio 0.7 \
  --evalue 1e-5 --min_align 50 \
  --out "$OUTDIR"

echo "[4/5] Converting mapping for PathLift"
CONVERT_ARGS=(
  --mapping "$RBH_MAPPING"
  --source-fasta "$HUMAN_FASTA"
  --out "$PATHLIFT_MAPPING"
)
if [[ -n "$TARGET_ID_MAP" ]]; then
  CONVERT_ARGS+=(--target-id-map "$TARGET_ID_MAP")
fi
"$PYTHON_BIN" "$REPO_ROOT/src/convert_rbh_to_pathlift.py" "${CONVERT_ARGS[@]}"

echo "[5/5] Done"
echo "RBH mapping table: $RBH_MAPPING"
echo "PathLift table: $PATHLIFT_MAPPING"
