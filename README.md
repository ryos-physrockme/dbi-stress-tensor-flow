# Stress-tensor flows in four-dimensional single-scalar DBI theory

Public research repository for a study of stress-tensor deformation flows in four-dimensional Dirac--Born--Infeld theory with one real scalar field and one Abelian gauge field.

The project combines an exploratory neural-network search with analytic derivations of the resulting stress-tensor flow. The scientific claims in the manuscript are established analytically; the machine-learning calculation is used to identify candidate structure and guide the search.

## Repository structure

- `paper/` — English manuscript (`paper.tex`, `paper.pdf`) and Japanese companion version (`paper_ja.tex`, `paper_ja.pdf`).
- `notes/` — detailed Japanese calculation notes (`calculation_note_ja.tex`, `calculation_note_ja.pdf`) with intermediate algebra and background explanations.
- `code/` — neural-network exploration, exact symbolic checks, and independent numerical validation.
- `data/` — numerical outputs used to validate the exploratory and analytic calculations.

## Reproducibility

The repository contains code for:

- the neural-network regression from stress-tensor invariants to the deformation generator;
- rational-function readout of the numerical data;
- exact symbolic checks of the four-dimensional matrix identities;
- independent numerical validation on generic field configurations;
- the five-dimensional Born--Infeld to four-dimensional DBI reduction checks.

By default, the scripts write outputs to the repository-level `data/` directory. This can be overridden with the `DBI_OUTPUT_DIR` environment variable.

## Status

Work in progress. The English manuscript and the detailed calculation notes are publicly available in this repository and are being prepared for preprint/submission.

## Author

Ryosuke Suda  
Graduate School of Science and Engineering, Saitama University  
255 Shimo-Okubo, Sakura-ku, Saitama 338-8570, Japan

## Licensing

This repository uses a mixed licensing scheme. See `LICENSE` for details. Source code under `code/` is released under the MIT License. Manuscript and calculation-note text remain copyright Ryosuke Suda unless a separate license is stated in a future preprint release.
