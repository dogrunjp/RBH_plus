# ヒト、ウシ、マウス、ラット、ボモのプロテオームに線虫とショウジョウバエを追加し、追加分のRBH + in-paralog を実行するスクリプト
# condsaのdiamond環境で実行することを想定
# 0) 3種のプロテオームを結合
cat /Users/oec/Desktop/docs/qp/data/2025_orthologs_analysis/GCF_000002985.6_WBcel235_protein.faa /Users/oec/Desktop/docs/qp/data/2025_orthologs_analysis/GCF_000001215.4_iso1_protein.faa /Users/oec/Desktop/docs/qp/data/2025_orthologs_analysis/Bomo_gene_models_prot.fa > ./data/Bomo_all_2025_12_v2.faa

# 1) DB 作成
diamond makedb -d Bomo_all_2025_12_v2 --in ./data/Bomo_all_2025_12_v2.faa

# 2) all-vs-all 検索 (DBとfastaは同じものを使用)
diamond blastp \
  -d Bomo_all_2025_12_v2 \
  -q ./data/Bomo_all_2025_12_v2.faa \
  --evalue 1e-5 --max-target-seqs 100 --threads 16 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen qstart qend sstart send \
  --header \
  -o ./data/Bomo_all_2025_12_v2_hits.tsv

# TODO: 対象生物が増えた場合の対応方法を検討すること

# 3) RBH + in-paralog 拡張
python3 rbh_plus_orthologs.py \
  --hits ./data/Bomo_all_2025_12_v2_hits.tsv \
  --species A ./data/tmp/cel_2025_12.ids \
  --species B ./data/tmp/fly_2025_12.ids \
  --species T ./data/tmp/bomo_2025_12.ids \
  --pairs A:T B:T \
  --tau 0.7 --min_cov 0.6 --len_ratio 0.7 \
  --evalue 1e-5 --min_align 50 \
  --out ./data/rbh/out_rbh_inparalog_bomo_v2