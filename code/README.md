# Code

This directory contains the reproducibility and exploratory code used in the four-dimensional single-scalar DBI stress-tensor-flow project.

## Files

- `nn_closure_search.py` — neural-network scout for the map from stress-tensor invariants to the dimensionless deformation generator, rational SVD readout, and independent numerical validation.
- `coordinate_free_sympy.py` — exact symbolic checks of the four-dimensional matrix identities used in the analytic derivation.
- Additional reduction/validation scripts will be added as the manuscript is frozen.

## Role of the neural network

The neural-network calculation is exploratory. It tests whether the deformation generator can be predicted from stress-tensor invariants and was used to motivate the search for a closed analytic expression. The scientific claims in the manuscript are established analytically and are intended to be independently verifiable without the trained network.

## Reproducibility

Python dependencies are listed in the repository-level `requirements.txt`. Numerical output files are stored under `data/` where appropriate. The repository will be made public together with the preprint, at which point a versioned archive/DOI can be added.
