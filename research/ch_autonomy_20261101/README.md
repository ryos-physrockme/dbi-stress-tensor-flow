# Local autonomy and global branch ambiguity in DBI

2026-11-01 の限定的な代数評価。既知の DBI 特性因数分解から lambda を消去した六次流れ多項式を得て、実際の実点での単純根を確認する。結論は四つの混合応力不変量を使う**分枝局所的**な自律形の存在であり、同一の完全混合応力と L を持つ二つの実 DBI 配置で Q が異なるため、時間的・空間的勾配を合わせた全物理領域上の単一の応力のみの自律関数は存在しない。固定した X 符号の領域や二つの二次不変量だけの CH 形への結論ではない。既知結果の代数的帰結であり、新規性は主張しない。前段階の成果物は変更しない。

## Files and provenance

- `note.tex`: 日本語の導出・限界・実点検査。XeLaTeX で `note.pdf` を生成する。
- `check.py`: 六対の因数分解と対称多項式、BI 重根、実点での単純根を確認する限定的 CPU 検査。
- `results.json`: 保存済み結果。再実行すると上書きされる。
- `metric_variation_check.py`: 対角計量の直接変分による Hilbert 応力対角成分、完全混合応力行列、固定場の直接微分 Q を別経路で照合する。
- Pinned commit: `c304d28526e7ed1b8deff80890eb1f7ac060f712`
- [Pinned manuscript](https://github.com/ryos-physrockme/dbi-stress-tensor-flow/blob/c304d28526e7ed1b8deff80890eb1f7ac060f712/paper/paper.tex)
- Known inputs: `eq:char-factor-app`, `eq:beta-app`, `eq:algebraic-relation`.

文献の追加照合と評価は [LITERATURE.md](LITERATURE.md) を参照。

スクリプトは原稿 blob を検証せず、以前のベースラインや追試も再実行しない。出典コミットはメタデータとして保存する。実行にネットワークは不要。上流リポジトリは変更しない。

## Run the check

```bash
cd path/to/this-directory
python check.py
python metric_variation_check.py
```

依存は SymPy。保存済み結果は Python 3.12.14、SymPy 1.14.0、CPU 上で 1.054 秒。**零残差19件**と、それとは別の**非零微分 assertion 2件**と、反例の区間正性などの assertions が通過した。一回の時間は性能保証ではない。前段階と同じ環境は Linux 6.18.57 x86_64、glibc 2.41、観測 CPU モデル `AMD EPYC 9V74 80-Core Processor`。モデル名から割り当て CPU 数を推定しない。

JSON の `H_lambda` と `K_Q` は多項式の文字列表現であり、微分ではない。微分の厳密値と近似値は `physical_point.H_derivative_nonzero`、`H_derivative_numeric`、`K_derivative_nonzero`、`K_derivative_numeric` に記録する。記録微分はそれぞれ約 `-1.78695607843948e-9`、`2.65874641025876e-13` だが、非零判定は厳密式で行う。数値条件の良さは主張しない。

## Build the note

XeTeX 0.999996（TeX Live 2025/dev/Debian）、`fontspec`、`amsmath`、`amssymb`、`mathtools`、`geometry`、`booktabs`、`hyperref` と、`Noto Serif CJK JP`、`Noto Sans CJK JP`、`DejaVu Sans Mono` を用いる。

```bash
xelatex -interaction=nonstopmode -halt-on-error note.tex
xelatex -interaction=nonstopmode -halt-on-error note.tex
pdfinfo note.pdf
pdftoppm -r 120 -png note.pdf note-page
```

日本語は XeLaTeX の `fontspec` と `\XeTeXlinebreaklocale "ja"` で処理する。`xeCJK` は不要。今回の環境には `luatexja` がないため LuaLaTeX は採用しない。配布前には全ページ画像で式、参照、欠字、改ページ、はみ出しを確認する。

## Scope and cautions

- `e_k` は混合応力 `Tmix` の基本対称式。`E_k` は `I-lambda*Tmix` の基本対称式であり、混同しない。スペクトル `t_i` は代数的には複素数を許し、全固有値の実数性や対角化を仮定しない。
- 定数項1の既知因子に属する対から `lambda=(t_i+t_j)/(t_i*t_j)`、`Q=-t_i*t_j/2` を得る。Q の積表示には積での除算は不要。lambda の再構成では `lambda!=0` と積の非零が除算条件。積が零なら両固有値が零で `Q=0`、この式で割らない。
- 六対の積で lambda が消えるが、対のラベルや物理的な根を選択する情報も落ちる。実際の実点で `dK/dQ!=0` なら実解析的陰関数定理により枝局所の `Q=f(e1,e2,e3,e4)` が存在する。これは全四不変量を使う結論であり、二不変量 CH NLED クラスや理論全体への一般主張ではない。
- 代入で生じる `lambda^6` は lambda=0 の物理分枝の証拠ではない。BI の重複スペクトルでは交差対の根が四重となり判定が失敗するが、既知の BI 自律流れは存在する。単純根条件は十分条件であって必要条件ではない。
- 実多項式の根だけでは DBI 配置の実現可能性を保証しない。`lambda>0`、`det(h)/det(g)>0`、`Delta>0` は必要であり、全物理領域の単一応力関数は下記の実反例で否定される。制限領域の分枝選択や継続は未検討。
- 実点は lambda=1、v=(1/3,1/4,0,1/5)、F01=1/5、F02=1/7、F12=1/6、F23=1/4。`X=-31/3600` で非 null。原稿の beta 比表示を null で割って使わない。
- 反例は同じ共通質量尺度で比較する。A: lambda=3/2, v=(1/2,0,0,0), F01=sqrt(3/8), Delta=1/4, Q=1/2。B: lambda=1, v=(0,1/sqrt(2),0,0), F01=sqrt(5)/2, Delta=1/2, Q=1/4。両方で Tmix=diag(-2,-1,1/2,1/2), L=1/2。各固定場の弱結合極限から端点まで h と G は非退化で行列式比は正。同じ lambda の衝突ではなく、固定場の力学の非一意性でもない。
- 数値では共通の固定質量 mu により lambda_hat=lambda*mu^4, v_hat=v/mu^2, F_hat=F/mu^2, T_hat=T/mu^4, L_hat=L/mu^4, Q_hat=Q/mu^8 として帽子を省略する。状態ごとの別単位を使わない。
- 次の制限領域の分枝分類・CH クラスの問いは今回実行していない。単純な次元勘定による不可能性も主張しない。

この比較は元の DBI 作用の Hilbert 応力、共通の結合定数規格化、固定 `(g,F,v)` での微分を用いる。改良応力や場に依存する結合定数の再定義についての主張ではない。実正則枝は `h` が Lorentz 計量で `Delta>0`、弱結合極限に連結する意味であり、動力学的・因果的安定性の判定ではない。局所的な一定場で反例が成立し、有限全エネルギーは要求しない。
