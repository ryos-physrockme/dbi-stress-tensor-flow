# 空間的勾配セクター内における DBI 流れの多価性

## 結論

四次元の単一実スカラーと電磁場を含む DBI 系について、同じ厳密な空間的勾配セクター `X=vᵀg⁻¹v>0` の内部で、完全な混合応力と L が一致しても、固定場微分 `Q=∂λL|g,Fcov,v` が異なる二つの実配置を構成した。X の数値そのものを一致させる主張ではない。

- 共通応力：Lorentz 整列後の `diag(-3,-2,1/2,7/10)`。四固有値はすべて単純。
- 共通 L：`1/2`。
- 配置 A：`λ=23/21`, `X=17/42`, `Δ=19/42`, `Q=21/20`。
- 配置 B：`λ=13/14`, `X=31/28`, `Δ=15/28`, `Q=7/10`。
- 差：`Q_A−Q_B=7/20`。

有理的な場の再構成により、目標固有値 `(-ρ,-s,p,q)` の `(3,2,1/2,7/10)` 近傍で、両分枝が同時に実現される。`L=p`, `Q_A=ρq/2`, `Q_B=sq/2`。Vandermonde Jacobian が非零なので、この重なりは独立な応力不変量の開近傍でもある。

各配置について場を固定し、λ=0 から端点までの全区間で正の DBI 平方根と Lorentz 計量 h を確認する。二配置は同じ連結な実正則 X>0 領域に属するが、それらを結ぶ経路は応力固有値の退化を通り得る。

**適用範囲**：分枝を指定しない単一の自律形 `Q=f(T)`、さらには `Q=f(T,L)` への障害である。λ が異なる比較なので `Q=f(λ,T,L)` や同一固定 λ での閉包は否定しない。各根を指定した局所表示とも両立する。安定性、追加の因果性条件、X<0 の一般分類、絶対的優先権は主張しない。

## ファイル

- `note.pdf`：自己完結した日本語科学ノート、A4 全7頁。
- `note.tex`：編集可能な XeLaTeX ソース。直接行列計算、一般再構成、開近傍の証明、明示的 Lorentz フレーム、変換後の場、文献上の範囲比較を含む。
- `check.py` / `results.json`：二配置の厳密な直接検算。記録済みの零残差72件。
- `reconstruction.py` / `reconstruction_results.json`：一般の再構成恒等式。零恒等式16件、基準点の厳密正値性、二つの非零 `K_Q`。
- `independent_check.py` / `independent_results.txt`：別実装による作用密度の計量10成分直接変分、固定場微分、因果型、実経路の検査。
- `LITERATURE.md`：一次文献との詳細な比較、調査範囲、評価、残る問い。
- `build_pdf.sh`：ローカル再ビルドとPNG描画。キャッシュと中間物は `qa/` 内。

再ビルド時に生成する PNG、TeX キャッシュ、ビルドログは品質確認用の中間物であり、数学的検算には不要。

## 再生成

XeLaTeX、Noto Serif CJK JP / Noto Sans CJK JP、DejaVu Sans Mono、Poppler の `pdftoppm` を用いる。Linux の TeX Live 標準パスを仮定した同梱スクリプトは、このフォルダ内だけに出力する。

```sh
./build_pdf.sh
```

他の TeX 環境では、フォントを導入したうえで通常のコマンドを使える。

```sh
xelatex -interaction=nonstopmode -halt-on-error note.tex
xelatex -interaction=nonstopmode -halt-on-error note.tex
pdftoppm -png -r 110 note.pdf qa/render/page
```

数学的検算には Python と SymPy が必要。既存記録は Python 3.12.14、SymPy 1.14.0。

```sh
python3 check.py
python3 reconstruction.py
python3 independent_check.py > independent_results.txt
```

最初の2コマンドは各 JSON と経過時間を更新する。72件と16件には互いに関連する成分恒等式が含まれ、独立な物理法則の個数ではない。これらの計算は文献上の優先権や摂動安定性を検証しない。ノートの作成工程では既存の数学的検算記録を読み、数式を導出・照合した。

## 入力と品質確認

基礎となる原稿と前段階のノートの固定版は `ryos-physrockme/dbi-stress-tensor-flow` の commit `1b546752fb96dc1ee65a3b7bb1827110b1d2eb20`。原稿 `paper/paper.tex` の Git blob は `4e78787d298956360ec2971bc06a1c5a973adc9a`。今回のノートはこの原稿を改変せず、追加の構成的結果として別途記した。

最終PDFは7頁。最終ビルドにエラー、LaTeX警告、overfull / underfull box はなく、110 dpi で描画した全7頁を目視確認した。欠字、数式の崩れ、重なり、切れ、過度なはみ出しは認めない。この作成工程でリモートファイルは変更していない。
