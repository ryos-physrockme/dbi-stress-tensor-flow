# Two-field null DBI and the CH bridge

2026-11-01 の二段階目の検算。非零 null 勾配 `v=(a,a,0,0)` と `F02=e,F12=b` の二成分族を直接検査し、既知の一般行列式・特性多項式恒等式を再利用して一般 null の CH 連鎖律から方向収縮を消去する。新公式や全行列代数の独立証明は主張しない。前段階の成果物は変更しない。

## Files and provenance

- `note.tex`: 日本語の導出ノート。XeLaTeX で `note.pdf` を生成する。
- `check.py`: 二成分行列の厳密検査、六成分の低次数多項式検査3件、既知恒等式の形式代入、極限・対照の検査。
- `results.json`: 実行結果。再実行すると上書きされる。
- Repository: `ryos-physrockme/dbi-stress-tensor-flow`
- Pinned commit: `c304d28526e7ed1b8deff80890eb1f7ac060f712`
- Manuscript Git blob SHA-1: `4e78787d298956360ec2971bc06a1c5a973adc9a`
- [Pinned manuscript](https://github.com/ryos-physrockme/dbi-stress-tensor-flow/blob/c304d28526e7ed1b8deff80890eb1f7ac060f712/paper/paper.tex)

文献の追加照合と評価は [LITERATURE.md](LITERATURE.md) を参照。

`check.py` は出典コミットと blob をメタデータとして保存するが、原稿ファイルの Git blob 検証や前段階ベースラインの再実行は行わない。実行にネットワークは不要。上流リポジトリは変更しない。

## Run the exact checks

```bash
cd path/to/this-directory
python check.py
```

依存は SymPy。記録環境は Python 3.12.14、SymPy 1.14.0、Linux 6.18.57 x86_64、glibc 2.41。`/proc/cpuinfo` で観測された CPU モデル名は `AMD EPYC 9V74 80-Core Processor` だが、割り当て CPU 数はこの名称から推定しない。最終実行は厳密残差27件がすべて零、所要時間44.454秒。25件は直接・特殊化・多項式検査、2件は既知恒等式の形式代入。各残差は `results.json` に記録する。CPU 上の記号計算であり、時間は性能保証ではない。

## Build the note

XeTeX 0.999996（TeX Live 2025/dev/Debian）、`fontspec`、`amsmath`、`amssymb`、`mathtools`、`geometry`、`booktabs`、`hyperref` と、`Noto Serif CJK JP`、`Noto Sans CJK JP`、`DejaVu Sans Mono` を用いる。

```bash
xelatex -interaction=nonstopmode -halt-on-error note.tex
xelatex -interaction=nonstopmode -halt-on-error note.tex
pdfinfo note.pdf
pdftoppm -r 120 -png note.pdf note-page
```

日本語には XeLaTeX の `fontspec` と `\XeTeXlinebreaklocale "ja"` を用いる。`xeCJK` は不要。今回の環境には `luatexja` がないため LuaLaTeX は採用しない。配布前には全ページ画像で数式、参照番号、欠字、はみ出し、改ページを確認する。

## Mathematical scope

- 記号対応: スクリプトの `P=Delta^2` をノートでは `R` と書き、電磁擬スカラー `P_em` と区別する。
- Two-field assumptions: `lambda>0,a>0,e,b real`, `P=1+lambda*(b^2-e^2)-lambda^2*k>0`, `k=a^2*(e-b)^2`, `Delta=positive sqrt(P)`.
- Direct two-field results: `Tvv=k/Delta`, metric correction `k/(2Delta)`, discriminant `lambda^4*k^2/Delta^2`.
- Controls: `e=b` gives `Delta=1,Q=C=0` even with nonzero fields; `e=-b` gives `P=1-4lambda^2*a^2*e^2`; `b=0` recovers the earlier transverse family; aligned `F01=e` has zero directional contraction and correction.
- Weak limit: `Q -> (b^2-e^2)^2/8+k/2`. Pure BI is recovered at `a=0`.
- General-null bridge: known manuscript determinant and characteristic-polynomial identities imply `Tvv=-Z/Delta`, `Z=-w^2<=0`, hence the positive discriminant root is `-lambda^2*Z/Delta`. The resulting invariant formula and starting trace-directional relation are already in the manuscript. Source labels and the imported inputs are explicit in the note.
- This closes the earlier general-null directional-elimination/branch issue **given those known identities**. It is not an independent derivation of all matrix identities, nor a general six-field inverse/EMT/discriminant symbolic test. The three delivered six-field checks cover only the canonical null Z polynomial, the amplitude determinant derivative, and t-linearity. Distinguish direct matrix residuals, low-degree polynomial tests, and formal substitutions of imported identities.
- Non-null CH invariant elimination independent of known formulas and a lambda-free autonomous form are the next bounded questions; neither is executed here. Literature comparison is separate, and no priority claim is made.
