# Reproducibility code

This directory contains the exploratory and verification code used in the four-dimensional single-scalar DBI stress-tensor-flow project.

## Files

- `nn_closure_search.py` — neural-network regression from stress-tensor invariants to the dimensionless deformation generator, affine-over-affine SVD readout, and numerical validation.
- `coordinate_free_sympy.py` — exact symbolic checks of the four-dimensional matrix identities.
- `coordinate_free_validation.py` — independent float64 validation in timelike, spacelike, null, and pure Born--Infeld sectors.
- `reduction_5d_to_4d_sympy.py` — exact and numerical checks of the five-dimensional Born--Infeld to four-dimensional DBI reduction, including the compactification-length normalization.

## Role of the neural network

The neural-network calculation is exploratory. It tests whether the deformation generator can be predicted from stress-tensor invariants and motivated the search for a closed analytic expression. The scientific claims in the manuscript are established analytically and do not depend on the trained network.

## Running the scripts

Python dependencies are listed in the repository-level `requirements.txt`. By default, generated tables, summaries, and figures are written to the repository-level `data/` directory. To use another output directory, set the environment variable `DBI_OUTPUT_DIR`, for example

```bash
DBI_OUTPUT_DIR=/tmp/dbi-output python code/coordinate_free_sympy.py
```

The scripts use fixed random seeds where stochastic sampling is involved.

## License

Unless otherwise stated, source code in this directory is released under the MIT License; see the repository-level `LICENSE` file.
