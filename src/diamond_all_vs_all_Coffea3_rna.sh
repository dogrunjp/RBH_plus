#!/usr/bin/env bash
set -euo pipefail

# Coffea arabica, Coffea canephora, and Coffea eugenioides all-vs-all BLASTN
# 0) 環境設定
# ./setup_conda_env.sh coffea_diamond
# conda activate coffea_diamond

cd "$(dirname "$0")"

command -v makeblastdb >/dev/null 2>&1 || {
  echo "Error: makeblastdb was not found. Run: conda activate coffea_diamond" >&2
  exit 1
}

command -v blastn >/dev/null 2>&1 || {
  echo "Error: blastn was not found. Run: conda activate coffea_diamond" >&2
  exit 1
}

mkdir -p ../output ./data/rbh/out_rbh_inparalog_coffea3

# 1) RNA/CDS FASTA を結合
#cd ./src # srcディレクトリに移動
#cat ../resource/*/*.fna > ../resource/merged.fna

# 2) DB 作成
makeblastdb \
  -in ../resource/merged.fna \
  -dbtype nucl \
  -out ../resource/merged_nucl

# 3) all-vs-all 検索 (DBとfastaは同じものを使用)
blastn \
  -db ../resource/merged_nucl \
  -query ../resource/merged.fna \
  -evalue 1e-5 -max_target_seqs 100 -num_threads 16 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore qlen slen qstart qend sstart send" \
  -out ../output/Coffea3_all_2026_6_hits_rna.tsv


# 4) RBH + in-paralog 拡張
# TODO: idsファイルは、各種の遺伝子IDを1行ずつ記載したテキストファイルで、species A, B, T の順に指定する
# extract_id_cols.sh で作成したものを使用することを想定している
python3 rbh_plus_orthologs.py \
  --hits ../output/Coffea3_all_2026_6_hits_rna.tsv \
  --species T ./data/C_arabica_rna.ids \
  --species A ./data/C_canephora_cds.ids \
  --species B ./data/C_eugenioides_rna.ids \
  --pairs A:T B:T \
  --tau 0.7 --min_cov 0.6 --len_ratio 0.7 \
  --evalue 1e-5 --min_align 50 \
  --out ./data/rbh/out_rbh_inparalog_coffea3
