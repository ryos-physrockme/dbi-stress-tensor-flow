# Detailed calculation notes

This directory is for the detailed Japanese calculation notes underlying the manuscript. The notes are intentionally more explicit than the paper and retain intermediate algebra, conventions, and reproducibility details needed to audit the derivation.

The current revision has been synchronized with the final manuscript. In particular, it includes:

- the corrected five-dimensional to four-dimensional reduction with distinct deformation parameters, compactification length `L_y`, and canonical zero-mode rescaling;
- the assumptions and continuous extension used in the characteristic-polynomial factorization;
- the null-gradient branch and pure Born--Infeld limit;
- the neural-network map, sampling prescription, loss, rational SVD readout, and independent validation;
- discussion of stress-tensor characterization of DBI theory and possible extensions.
- the Courant--Hilbert chapter in `ch_dbi_construction_ja.tex`, included by `calculation_note_ja.tex`; it derives the induced-metric Born--Infeld decomposition, the stress-tensor relation, and the parameter-dependent metric chain rule.

Compile the current TeX source from the `notes/` directory using `latexmk -lualatex calculation_note_ja.tex`. The tracked `calculation_note_ja.pdf` has not yet been regenerated on this source-only branch. Exact checks are in `code/ch_dbi_sympy.py`.
