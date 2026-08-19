# Stress-tensor flows in four-dimensional single-scalar DBI theory

Public research repository for a study of stress-tensor deformation flows in four-dimensional Dirac--Born--Infeld theory with one real scalar field and one Abelian gauge field.

The project combines an exploratory neural-network search with analytic derivations of the resulting stress-tensor flow. The scientific claims in the manuscript are established analytically; the machine-learning calculation is used to identify candidate structure and guide the search.

## Repository structure

- `paper/` — manuscript status and submission-related material. The current Japanese manuscript source will be synchronized after the public-review draft is frozen.
- `notes/` — detailed calculation-note material and documentation.
- `code/` — neural-network exploration, exact symbolic checks, and independent numerical validation.
- `data/` — numerical outputs used to validate the exploratory and analytic calculations.

## Reproducibility

The repository contains code for:

- the neural-network regression from stress-tensor invariants to the deformation generator;
- rational-function readout of the numerical data;
- exact symbolic checks of the four-dimensional matrix identities;
- independent numerical validation on generic field configurations;
- the five-dimensional Born--Infeld to four-dimensional DBI reduction checks.

## Status

Work in progress. The Japanese manuscript is undergoing public review before preparation of an English preprint/submission version.

## Author

Ryosuke Suda  
Graduate School of Science and Engineering, Saitama University  
255 Shimo-Okubo, Sakura-ku, Saitama 338-8570, Japan

## Licensing

This repository uses a mixed licensing scheme. See `LICENSE` for details. In brief, source code under `code/` is released under the MIT License. Manuscript and calculation-note text remain copyright Ryosuke Suda unless a separate license is stated in a future preprint release.
