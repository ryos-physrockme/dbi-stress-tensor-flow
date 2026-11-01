# DBI 流れの二分枝連続族と Maxwell 極限

2026-11-01。日本語の自己完結型ノートと小規模な厳密検算。

## 成果物

- `note.tex`: 編集可能な日本語 XeLaTeX ソース。
- `note.pdf`: 読むための 4 頁 PDF。最終版は全頁の画像確認済み。
- `check.py`: SymPy による主検算。
- `results.json`: 同検算の保存結果。成分ごとの厳密零残差 50 件、別途区間正値性の assertion。
- `independent_check.py`: 別実装による端点代入前の計量変分・固定場微分と成分別の検査。
- `independent_check.txt`: 別実装の保存結果。成分ごとの厳密零残差 43 件。
- `LITERATURE.md`: 一次文献との適用範囲比較。
- `qa/`: 非公開の取得済み原稿・前段ノートと、PDF のビルド／描画／確認記録。公開パッケージには含めない。

## 結論の範囲

共通の固定単位と結合定数規格化のもと、時間的・空間的スカラー勾配の二つの軸方向電場配置が、同じ完全混合応力 `diag(-rho,-s,q,q)` と同じ `L=q` を持ち、異なる固定場の `Q=partial_lambda L` を持つ連続族を与える。両者は同じ非零 Maxwell 背景へ結合定数を零にして近づく。電場を零にする極限ではなく、両方の Q は同じ q²/2 に収束するため不連続性も主張しない。従って、両セクターを含むその点の全場・結合定数近傍に単一の応力のみの自律関数は存在しない。

これは `Q=f(T;lambda)`、固定した `X` 符号の一般領域、一般の単純スペクトル、安定性に対する否定ではない。epsilon に沿う弱結合極限は場も変え、Q は各 epsilon の配置で別々に固定場微分する。新規性は未確認。

## 再生成

Noto Serif CJK JP と Noto Sans CJK JP を備えた XeLaTeX 環境で、このディレクトリから次を実行する。

```sh
mkdir -p qa/build qa/render
xelatex -interaction=nonstopmode -halt-on-error -output-directory=qa/build note.tex
xelatex -interaction=nonstopmode -halt-on-error -output-directory=qa/build note.tex
cp qa/build/note.pdf note.pdf
pdftoppm -png -r 120 note.pdf qa/render/page
```

数学的検算の再実行は次のとおり。

```sh
python check.py
python independent_check.py
```

第1コマンドは `results.json` の経過時間を更新する。既存記録は Python 3.12.14、SymPy 1.14.0、0.363 秒。第2コマンドは `independent_check.txt` を更新する。2つの検算は重複を含み、50件と43件は独立な物理恒等式の数ではない。ノート作成工程では検算を再実行していない。

## 固定出典

- 取得基点: `ryos-physrockme/dbi-stress-tensor-flow`, commit `45c11bc3f2396bcb523a6e4906bba89316e16ce9`。
- 取得済み原稿: `qa/paper.tex`。
- 取得済み前段ノート: `qa/note.tex`。
- 非公開の `qa/source_manifest.json` に取得した3ファイルの Git blob SHA を記録。
- 文献比較は範囲の確認に限定し、新規性を確定しない。

## 品質確認

最終ビルドは警告なし。120 dpi の全 4 頁 PNG を開き、欠字、数式の崩れ、重なり、はみ出しがないことを確認した。非公開の確認記録は `qa/visual_qa.md` に保存してある。
