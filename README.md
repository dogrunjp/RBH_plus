# RBH plus

DIAMOND `blastp` の all-vs-all 検索結果から、レシプロカルベストヒット
(RBH) と in-paralog を含む簡易的な ortholog/paralog テーブルを書き出すための
作業用ツールです。

## 概要

基本的な流れは次の通りです。

1. 比較したい生物種のプロテオーム FASTA を結合する
2. DIAMOND で protein DB を作成する
3. 同じ FASTA/DB に対して `diamond blastp` の all-vs-all 検索を行う
4. `rbh_plus_orthologs.py` で RBH と in-paralog を推定する

現状は汎用パイプラインというより、`src/diamond_all_vs_all_bomo2.sh` を雛形にして
プロテオーム、ID リスト、出力先を置き換えながら使う想定です。

## セットアップ

Conda または mamba が使える環境で、必要なツールをまとめて入れます。

```bash
./setup_conda_env.sh
conda activate rbh_plus_diamond
```

主な依存関係は DIAMOND、Python、pandas、numpy です。

## 使い方

まず、`src/diamond_all_vs_all_bomo2.sh` をコピーするなどして、以下を自分の解析対象に
合わせて編集します。

- 結合するプロテオーム FASTA
- `diamond makedb` の DB 名
- `diamond blastp` の入力 FASTA と出力 TSV
- 生物種ごとの protein ID リスト
- 比較したいペア指定
- `--out` の出力ディレクトリ

DIAMOND の出力形式は、少なくとも次の列を含む必要があります。

```text
qseqid sseqid pident length evalue bitscore qlen slen qstart qend sstart send
```

`rbh_plus_orthologs.py` の実行例:

```bash
python3 src/rbh_plus_orthologs.py \
  --hits data/all_vs_all_hits.tsv \
  --species A data/species_a.ids \
  --species B data/species_b.ids \
  --species T data/target.ids \
  --pairs A:T B:T \
  --tau 0.7 --min_cov 0.6 --len_ratio 0.7 \
  --evalue 1e-5 --min_align 50 \
  --out output/rbh_result
```

## 出力

指定した `--out` ディレクトリに、ペアごとの RBH + in-paralog テーブルが出力されます。

- `{A}:{B}.rbh_inparalog.tsv`: ペアごとの ortholog group
- `map_{src}_to_T.tsv`: `T` ラベルを指定した場合の target への対応表

## 注意

- `resource/` と `output/` は Git 管理対象外です。
- 大きな FASTA、DIAMOND DB、検索結果 TSV はローカルに置いて使う想定です。
- スクリプト内のパスは解析環境に合わせて書き換えてください。
- in-paralog の判定は bitscore 比 `--tau`、coverage、length ratio などのしきい値に依存します。
