# Reproducibility code

This directory contains the exploratory and verification code used in the four-dimensional single-scalar DBI stress-tensor-flow project.

## Files

- `nn_closure_search.py` — neural-network regression from stress-tensor invariants to the dimensionless deformation generator, rational SVD readout, and independent numerical validation.
- `coordinate_free_sympy.py` — exact symbolic checks of four-dimensional matrix identities.
- Additional reduction/validation scripts are kept here as the public manuscript is synchronized with the repository.

## Role of the neural network

The neural-network calculation is exploratory. It tests whether the deformation generator can be predicted from stress-tensor invariants and motivated the search for a closed analytic expression. The scientific claims in the manuscript are established analytically and do not depend on the trained network.

## Reproducibility

Python dependencies are listed in the repository-level `requirements.txt`. Numerical output files are stored under `data/`.

## License

Unless otherwise stated, source code in this directory is released under the MIT License; see the repository-level `LICENSE` file.
