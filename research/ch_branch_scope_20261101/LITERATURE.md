# DBI 自律流れの分枝条件と先行研究（2026-11-01）

## 結論と今回の差分

既存の「同じ完全な混合 Hilbert 応力 T と L、異なる固定場微分 Q=∂λL」という二配置は、`T=diag(-ρ,-s,q,q), ρ>s>q>0` の連続族に拡張できる。時間的勾配の枝で `Q=ρq/2`、空間的勾配の枝で `Q=sq/2`。`ρ=q+2ε, s=q+ε` とすれば、両枝は同じ Maxwell 配置へ接近し、両方の λ が 0 に近づく。したがって非一意性を強結合領域だけの問題として除外できない。ここでは q>0 を固定するので非零 Maxwell 背景への弱結合極限であり、電磁場の振幅が零になる極限ではない。両方の Q は q²/2 に収束し、差も零に近づくため、不連続性を主張していない。これは ε に沿う全微分ではなく、各配置で場を固定した λ 微分の比較である。

この族は横方向に二重固有値を持つ。各因果型を固定した**今回の軸対称・純電場 ansatz の内部**では、応力から非零 λ が一意に再構成される。固定 X 符号の一般的な領域、四固有値が相異なる領域、あるいは追加の分枝条件を課す領域の非一意性は示していない。

## 一次本文との照合

| 一次資料 | 確認箇所と適用域 | 今回の結果との関係 |
|---|---|---|
| Ferko, Kuzenko, Smith, Tartaglino-Mazzucchelli (2023), *Duality-Invariant Non-linear Electrodynamics and Stress Tensor Flows*, [arXiv:2309.04253v4](https://arxiv.org/html/2309.04253v4) | §5, Theorem 3, 式(5.9)–(5.12)。純電磁の双対不変族。局所逆関数で得る流れは明示的な結合定数を含み得る。 | 同定理は全領域の λ 自由な自律性を保証しない。実スカラーを含む今回の模型も直接の対象外。 |
| Pezhman, Babaei Velni, Babaei-Aghbolagh (2026), *Exact Stress Tensor flows in Causal Electrodynamics: Courant-Hilbert Approach*, [arXiv:2610.33887v1](https://arxiv.org/html/2610.33887v1) | §2.2–2.3。一尺度の自己双対な純 NLED において、不変量比を反転できる分枝で自律流れを構成。 | 反転可能な分枝という条件を落としてスカラー＋電磁 DBI 全域に拡張することはできない。 |
| Pezhman, Babaei-Aghbolagh, Babaei Velni (2026), *Exact Stress Tensor Flows in Causal Electrodynamics: Auxiliary Field Formalism*, [arXiv:2610.35135v1](https://arxiv.org/html/2610.35135v1) | §4, 式(42)–(44)。純 NLED で因果性と大域単射性を区別。同じ正の不変量比に異なる補助変数が対応すると核の値が異なることを説明。 | 局所／大域の区別と単射性の障害は既知。今回の差分は特定の4D DBIにおける完全なTの一致と、同じ弱結合配置へ接近する実現可能な連続族。 |
| Brizio, Kade, Sfondrini, Sorokin (2026), *More on TbarT-like deformations in higher dimensions*, [arXiv:2602.24058v2](https://arxiv.org/html/2602.24058v2) | §5.2末尾、式(119)の後。4D単一スカラーの場合に流れを「T と λ」によって表せるかを将来課題として明記。§6も同課題を挙げる。 | 文献の問いは明示的 λ 依存を許す。今回の λ 自由な大域関数の反証は、この期待の否定ではない。 |
| Blair, Lahnsteiner, Obers, Yan (2025), *Matrix Theory Reloaded: A BPS Road to Holography*, [arXiv:2410.03591](https://arxiv.org/html/2410.03591) | §6, 式(212)–(217)。スカラーとゲージ場の混合系で、4Dでは電磁場不変量を含む摂動的流れを提示。 | この表式で追加不変量が現れることだけから、別の応力表示が存在しないとは結論できない。今回の厳密な一致配置は、より限定された大域自律性の問いを直接判定する。 |
| 既存 DBI 原稿・検算ノート | [固定原稿](https://github.com/ryos-physrockme/dbi-stress-tensor-flow/blob/45c11bc3f2396bcb523a6e4906bba89316e16ce9/paper/paper.tex)、[前段階ノート](https://github.com/ryos-physrockme/dbi-stress-tensor-flow/blob/45c11bc3f2396bcb523a6e4906bba89316e16ce9/research/ch_autonomy_20261101/note.tex)。元の Hilbert 応力と作用を継承。 | 前段階で個別反例と単純根近傍の局所自律形を確認済み。今回の連続族は既知作用の小行列計算であり、新しい作用や独立な一般閉包定理ではない。 |

## 新規性についての調査範囲

2026-11-01 に上記4近接論文を本文で再照合し、2602.24058 の引用先 2410.03591 の §6 まで遡った。検索語は `"Dirac-Born-Infeld" "autonomous" "stress" scalar`、`"Born-Infeld" "same" "stress tensor" "flow"`、`"Dirac-Born-Infeld" "autonomous" "counterexample"`、`"DBI" "stress" "timelike" "spacelike" flow`、`"Born-Infeld" "same stress"`。さらに `Born Infeld scalar electromagnetic energy momentum tensor inverse reconstruction Rainich`、`Dirac Born Infeld stress tensor eigenvalues scalar gradient reconstruction`、`"Born-Infeld" "stress tensor" "non-uniqueness"` で逆問題としての同等表現も検索した。これらの追加検索では直接対応する一次結果を特定できなかった。検索結果の第三者要約は主張の根拠に用いていない。

今回確認した本文には、この時間的／空間的二配置の厳密な一致族または同じ弱結合点への接近による同等の反証を確認しなかった。ただし、限定的な検索であり、先行結果が存在しない証明ではない。**新規性・優先権は未確定**である。一般的な逆関数・分枝・単射性の問題を新しい発見と主張しない。

## 適用範囲の分類

| 問い | 現状 | 根拠・未解決事項 |
|---|---|---|
| 元の DBI 族で Q=f(T;λ) と表せるか | 既存原稿の結果を維持 | 今回は同じ λ の衝突を示していない。 |
| 選択した単純根の近傍で Q=f(e1,…,e4) と表せるか | 前段階で十分条件を確認 | 陰関数定理。根の選択を含む局所的な主張。 |
| 両因果型を含む実正則 DBI 全域で単一値 Q=f(T) が存在するか | 否定 | 同じ完全なT、異なるQ。Lを追加しても同じ。 |
| この障害を λ が十分小さい領域から除けるか | 両因果型を含む Maxwell 点の全近傍では除けない | 任意の ε>0 で差が残り、ε→0で両λ→0、場も同じ極限へ接近。 |
| 今回の軸対称 ansatz の一方の因果型のみ | λ とQが応力から一意 | それぞれ異なる相反固有値対の選択で決まる。 |
| 固定 X 符号の一般の4D実配置 | 未解決 | 今回のansatzを超える場の向き・電磁成分を含む分類が必要。 |
| 四固有値が相異なる領域 | 未解決 | 反例の横方向固有値qは二重。単純根と単純スペクトルも同じ条件ではない。 |

実正則とは正の結合定数、Lorentz 型の誘導計量、正の Δ、および自由極限までの非退化性を意味する。動力学的・因果的安定性は検証していない。改良応力、場依存の結合定数再定義、有限全エネルギーの要請は対象外。

## 五軸による再評価

| 軸 | 評価 | 根拠と前段階からの変化 |
|---|---|---|
| 新規性 | 判断保留 | 正確な近接定理との適用域の違いを明示できたが、優先権は未確認。 |
| 進歩性 | 限定的だが具体的 | 小行列の標準的計算による一般化。単独の数値例に依存せず、任意に弱い結合で障害が残る条件を明確化した。 |
| 実現可能性 | 今段階は完了 | 数値探索なし。SymPy厳密検査50件と区間正性の証明で閉じる。 |
| 物理的な面白さ | 変形の枝選択に関して有用 | 同じMaxwell極限と同じ有限結合の応力値でも、両因果型を統合する単一の変形方向は得られない。固定場の時間発展の非一意性ではない。 |
| 筋の良さ | 元原稿の補足として適切 | 明示的λ流れを保持し、自律化の適用域を切り分ける。独立した大規模な一般定理としては未成熟。 |

## 今回の終了点と次の未解決事項

連続族・弱結合極限・各制限ansatz内の一意性までで終了する。一般の固定 X 符号領域を扱うには、固有値対の代数条件に加えて実 DBI 場としての実現可能性を分類する必要がある。その分類や広い探索は今回開始していない。
