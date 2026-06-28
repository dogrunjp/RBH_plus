#!/usr/bin/env bash
# Coffea arabica parental origin classification from BLAST hits
# Coffea arabicaのBLASTヒットから親由来を分類するスクリプト相手はCoffea canephoraとCoffea eugenioidesの2種
set -euo pipefail

cd "$(dirname "$0")"

TARGET_CSV="${1:-../caffeine_pathway_gene.csv}"
ARABICA_FASTA="../resource/GCF_036785885.1/rna.fna"
CANEPHORA_FASTA="../resource/DH200_94_v2/Coffea_canephora_DH200-94.cds.fna"
EUGENIOIDES_FASTA="../resource/GCF_003713205.1/rna.fna"

WORK_DIR="../output/parental_origin_caffeine"
TARGET_FASTA="${WORK_DIR}/C_arabica_caffeine_targets.fna"
PARENT_FASTA="${WORK_DIR}/parents_canephora_eugenioides.fna"
PARENT_DB="${WORK_DIR}/parents_nucl"
HITS_TSV="${WORK_DIR}/C_arabica_targets_vs_parents.blastn.tsv"
RESULT_TSV="${WORK_DIR}/C_arabica_parental_origin.tsv"

command -v makeblastdb >/dev/null 2>&1 || {
  echo "Error: makeblastdb was not found. Run: conda activate coffea_diamond" >&2
  exit 1
}

command -v blastn >/dev/null 2>&1 || {
  echo "Error: blastn was not found. Run: conda activate coffea_diamond" >&2
  exit 1
}

mkdir -p "${WORK_DIR}"

python3 extract_fasta_by_csv_ids.py \
  --csv "${TARGET_CSV}" \
  --id-column transcript_id \
  --fasta "${ARABICA_FASTA}" \
  --out "${TARGET_FASTA}"

awk '/^>/ {sub(/^>/, ">canephora|")} {print}' "${CANEPHORA_FASTA}" > "${PARENT_FASTA}"
awk '/^>/ {sub(/^>/, ">eugenioides|")} {print}' "${EUGENIOIDES_FASTA}" >> "${PARENT_FASTA}"

makeblastdb \
  -in "${PARENT_FASTA}" \
  -dbtype nucl \
  -out "${PARENT_DB}"

blastn \
  -task megablast \
  -db "${PARENT_DB}" \
  -query "${TARGET_FASTA}" \
  -evalue 1e-5 \
  -max_target_seqs 20 \
  -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore qlen slen qstart qend sstart send" \
  -out "${HITS_TSV}"

python3 classify_parental_origin.py \
  --targets "${TARGET_CSV}" \
  --hits "${HITS_TSV}" \
  --out "${RESULT_TSV}" \
  --min-cov 0.7 \
  --min-pident 85 \
  --min-delta-bitscore 20 \
  --min-ratio 1.02

echo "Wrote: ${RESULT_TSV}"
