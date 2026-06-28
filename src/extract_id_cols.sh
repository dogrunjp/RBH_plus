# blastに利用した元ファイルからID列を抽出するスクリプト
# ファイルパスと出力するファイル名は適宜変更してください
awk '/^>/{sub(/^>/,""); print $1}' ../resource/DH200_94_v2/Coffea_canephora_DH200-94.cds.fna > ./data/C_canephora_cds.ids
awk '/^>/{sub(/^>/,""); print $1}' ../resource/GCF_036785885.1/rna.fna > ./data/C_arabica_rna.ids
awk '/^>/{sub(/^>/,""); print $1}' ../resource/GCF_003713205.1/rna.fna > ./data/C_eugenioides_rna.ids