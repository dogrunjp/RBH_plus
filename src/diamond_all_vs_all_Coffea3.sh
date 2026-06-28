# Coffea arabica, Coffea canephora, and Coffea eugenioides all-vs-all BLASTP
# 0) 環境設定
# ./setup_conda_env.sh coffea_diamond
# conda activate coffea_diamond
# 1) プロテオームを結合
cd ./src # srcディレクトリに移動
cat ../resource/*/protein.faa > ../resource/merged.faa

# 2) DB 作成
diamond makedb -d ../resource/merged --in ../resource/merged.faa

# 3) all-vs-all 検索 (DBとfastaは同じものを使用)
diamond blastp \
  -d ../resource/merged \
  -q ../resource/merged.faa \
  --evalue 1e-5 --max-target-seqs 100 --threads 16 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen qstart qend sstart send \
  --header \
  -o ../output/Coffea3_all_2026_6_hits.tsv


# 4) RBH + in-paralog 拡張
# TODO: idsファイルは、各種の遺伝子IDを1行ずつ記載したテキストファイルで、species A, B, T の順に指定する
# extract_id_cols.sh で作成したものを使用することを想定している
python3 rbh_plus_orthologs.py \
  --hits ./output/Coffea3_all_2026_6_hits_rna.tsv \
  --species T ./data/C_arabica_rna.ids \
  --species A ./data/C_canephora_rna.ids \
  --species B ./data/C_eugenioides_rna.ids \
  --pairs A:T B:T \
  --tau 0.7 --min_cov 0.6 --len_ratio 0.7 \
  --evalue 1e-5 --min_align 50 \
  --out ./data/rbh/out_rbh_inparalog_coffea3_rna_v2


# プロテオームだと発現解析の結果と紐付けるのが難しいので、mRNAの配列でやり直す（C.canephoraはrna.fnaが含まれないため別ソースから取り直す）