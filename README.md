# Vanilla-CKKS
修論のためにCKKSの暗号化、復号を理解する必要があったため作成しました。私の英語力は赤ちゃん未満なので間違ったコメントをつけている可能性もあります。ご了承ください。

## このページで分かること
・CKKSの暗号化、復号の最適化がされていない実装

・CKKSの各処理における入出力

・CKKSのパラメータの例（研究目的のため、1パターンしか扱いません）

・実装に必要なざっくりとした説明

## このページで触れないこと
・数学的な証明および論文等の解説

・既存ライブラリで行われているような最適化

・パラメータの厳密な設定方法（どこにも書いてなかったものもある）

## 参考文献
CKKS explained: Part 1, Vanilla Encoding and Decoding(https://openmined.org/blog/ckks-explained-part-1-simple-encoding-and-decoding/)

準同型暗号CKKSその1 多項式環(https://zenn.dev/herumi/articles/ckks-ring-iso)
