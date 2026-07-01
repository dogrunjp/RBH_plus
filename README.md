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

## ヒトからハダカデバネズミへの対応表の作成例

ヒトとハダカデバネズミ（Taxonomy ID: 10181, 学名: Heterocephalus glaber）のプロテオーム FASTA がある場合、次のように実行すると対応表を作成できます。

```bash
chmod +x src/run_human_to_hglaber.sh
./src/run_human_to_hglaber.sh \
  /path/to/human_proteins.faa \
  /path/to/hglaber_proteins.faa \
  ./output/human_hglaber
```

出力先の [output/human_hglaber](output/human_hglaber) 配下に、以下のファイルが生成されます。

- `map_H_to_T.tsv`: ヒト ID からハダカデバネズミ ID への対応表
- `pathlift_H_to_T.tsv`: PathLift の `provided_table` として読める対応表
- `human_hglaber_hits.tsv`: DIAMOND の all-vs-all hits
- `human.ids`, `hglaber.ids`: 各種の ID リスト

必要に応じて、同じ手順で別の種間でも再利用できます。

### この作成例のために行った内部実装の変更

今回の作成例では、次のコード変更と追加を行っています。

- `src/rbh_plus_orthologs.py`: DIAMOND の検索結果全体を pandas の
  `DataFrame` に読み込む実装から、TSV を1行ずつ読み込み、条件を満たしたヒットの
  ベストヒットと in-paralog 展開用インデックスだけを保持する実装に変更しました。
  ヒトとハダカデバネズミの all-vs-all 検索結果は大きくなるため、解析時のメモリ消費を
  抑える必要があったためです。これに伴い、このスクリプトでの pandas／NumPy 依存も
  なくし、結果TSVはPython標準ライブラリで出力するようにしました。
- `src/run_human_to_hglaber.sh`: FASTAからのIDリスト作成、FASTAの結合、DIAMOND DBの
  作成と検索、RBH + in-paralog解析、PathLift用TSVへの変換を順番に実行するラッパーを
  追加しました。手作業による入力ファイルや種ラベル、出力先の指定間違いを避け、同じ
  条件で対応表を再作成できるようにするためです。
- `src/convert_rbh_to_pathlift.py`: `map_H_to_T.tsv` のカンマ区切りの対応候補を1候補1行に
  展開し、ヒトFASTAヘッダーから遺伝子記号を取得して、PathLiftの
  `target_id`, `source_symbol`, `source_pid` 形式に変換する処理を追加しました。
  RBHの元の出力形式をPathLiftへ直接入力できないために必要です。また、任意の
  protein-to-gene TSVを使ってハダカデバネズミのprotein accessionをgene IDへ変換し、
  複数アイソフォームを遺伝子単位にまとめられるようにしました。

## PathLift 用TSVへの変換

`convert_rbh_to_pathlift.py` は `map_*_to_T.tsv` のカンマ区切り候補を1候補1行に展開し、
source FASTA のヘッダーから遺伝子記号を補完します。出力列は
`target_id`, `source_symbol`, `source_pid` です。

```bash
python3 src/convert_rbh_to_pathlift.py \
  --mapping output/human_hglaber/map_H_to_T.tsv \
  --source-fasta data/human_proteins.faa \
  --out output/human_hglaber/pathlift_H_to_T.tsv
```

この実行では対象ID（`XP_...` / `NP_...`）をprotein accessionのまま出力します。
PathLiftで遺伝子単位に集約する場合は、次のヘッダーを持つprotein-to-gene TSVを渡します。

```text
protein_id\tgene_id
XP_004875082.1\t100123456
```

```bash
python3 src/convert_rbh_to_pathlift.py \
  --mapping output/human_hglaber/map_H_to_T.tsv \
  --source-fasta data/human_proteins.faa \
  --target-id-map data/hglaber_protein_to_gene.tsv \
  --out output/human_hglaber/pathlift_H_to_T.tsv
```

`run_human_to_hglaber.sh`では変換も自動実行します。protein-to-gene TSVがある場合は
第4引数に指定できます。

```bash
./src/run_human_to_hglaber.sh \
  data/human_proteins.faa \
  data/hglaber_proteins.faa \
  output/human_hglaber \
  data/hglaber_protein_to_gene.tsv
```

PathLift recipeの列指定は次のとおりです。

```yaml
provided_table:
  path: /path/to/pathlift_H_to_T.tsv
  format: tsv
  columns:
    target_id: target_id
    source_symbol: source_symbol
    source_pid: source_pid
```

## 注意

- `resource/` と `output/` は Git 管理対象外です。
- 大きな FASTA、DIAMOND DB、検索結果 TSV はローカルに置いて使う想定です。
- スクリプト内のパスは解析環境に合わせて書き換えてください。
- in-paralog の判定は bitscore 比 `--tau`、coverage、length ratio などのしきい値に依存します。
