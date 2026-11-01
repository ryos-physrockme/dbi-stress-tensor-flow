# DBI CH pinned follow-up

2026-11-01 の追試。上流の式を新規結果として扱わず、計量の連鎖律と電場一成分の非零 null 族での平方根の枝を導出・検査する。

## Files

- `note.tex`: 日本語数式ノートの編集可能なソース。XeLaTeX で `note.pdf` を生成する。
- `check.py`: 固定ソースの Git blob SHA-1 を検証し、ベースラインと 18 件の追加厳密残差を実行する。
- `source_ch_dbi_sympy.py`: 上流 `code/ch_dbi_sympy.py` の固定コピー。末尾改行を含めて変更しない。
- `results.json`: 保存済みの実行結果。再実行すると上書きされる。

## Provenance

- Repository: `ryos-physrockme/dbi-stress-tensor-flow`
- Branch: `feature/ch-dbi-construction`
- Commit: `62e0d609e6d5109675834b0d5f14e5afdfc47550`
- Git blob SHA-1: `8f5056aa25bc51d10f03e1aea3d64df2ae435d37`
- [PR 5](https://github.com/ryos-physrockme/dbi-stress-tensor-flow/pull/5)
- [Pinned source](https://github.com/ryos-physrockme/dbi-stress-tensor-flow/blob/62e0d609e6d5109675834b0d5f14e5afdfc47550/code/ch_dbi_sympy.py)

文献の追加照合と評価は [LITERATURE.md](LITERATURE.md) を参照。

`check.py` は Git の `blob <byte count>\0<bytes>` 形式を SHA-1 に入力し、完全一致を要求する。外部リポジトリは変更しない。

## Run the exact checks

保存済み結果の環境は Python 3.12.14、SymPy 1.14.0、Linux 6.18.57 x86_64、glibc 2.41。`/proc/cpuinfo` で観測された CPU モデル名は `AMD EPYC 9V74 80-Core Processor` だが、割り当て CPU 数はこの名称からは推定しない。適切な Python 環境から実行する。

```bash
cd path/to/this-directory
python check.py
```

依存は SymPy。ネットワークは不要。保存済み結果は `status: all_passed`、ベースライン表示 7 群の PASS と、それとは別の追加残差 18 件がすべて厳密零、記録時間 5.204 秒。CPU 上の記号計算であり、時間は環境依存の単一実行記録である。

## Build the note

使用エンジンは XeTeX 0.999996（TeX Live 2025/dev/Debian）。XeLaTeX と `fontspec`、`amsmath`、`amssymb`、`mathtools`、`geometry`、`booktabs`、`hyperref`、フォント `Noto Serif CJK JP`、`Noto Sans CJK JP`、`DejaVu Sans Mono` を用いる。

```bash
xelatex -interaction=nonstopmode -halt-on-error note.tex
xelatex -interaction=nonstopmode -halt-on-error note.tex
pdfinfo note.pdf
pdftoppm -r 120 -png note.pdf note-page
```

日本語には XeLaTeX の `fontspec` と `\XeTeXlinebreaklocale "ja"` を使用する。`xeCJK` は不要。今回の環境には `luatexja` がないため LuaLaTeX は採用していない。配布前には全ページ画像を目視し、参照番号、数式、改ページ、欠字、はみ出しを確認する。

## Scope and limitations

- 一般の Jacobi/トレース恒等式は、非退化行列と正の行列式比の仮定の下で導出する。方向収縮 `v^T T v` を含むため、original-EMT 不変量だけの閉形式の証明ではない。
- null 根式は固定ソースに既存。本追試の独立再導出は `v=(a,a,0,0)`、`F02=e` のみ、`lambda,a,e>0`、`1-lambda*e^2-lambda^2*a^2*e^2>0` に限定する。
- `q=1` でも誘導計量は `lambda` に依存する。計量補正は `a^2 e^2/(2D)` で、一般に非零。
- 対照の整列電場 `F01=e` では、同じ非零 null 勾配でも方向収縮と計量補正は零。非零補正をすべての null 配置へ一般化しない。横方向族の弱結合極限は `a^2 e^2/2+e^4/8`。
- 次の課題は独立ゲージ成分二つの null 族、その後に一般 null 不変量恒等式と枝の符号の証明。今回それらは実行していない。
- 有限個の配置や残差の通過は任意配置の定理や物理的安定性の証明ではない。元の EMT だけによる有理流れは固定ソースに既存。未解決なのはその存在ではなく、CH 連鎖律からの独立な方向収縮の一般消去と、lambda を陽に含まない自律形の構成である。`note.tex` 自体は論文本体の照合や文献新規性評価を目的とせず、制限族の計算を記録する。追加の文献照合は別紙に分離した。
