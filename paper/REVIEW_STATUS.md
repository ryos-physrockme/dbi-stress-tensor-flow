# Referee-style manuscript checklist

## Addressed in the current Japanese draft

- Define the deformation derivative as the derivative of the DBI Lagrangian with respect to the deformation parameter at fixed metric, field strength, and scalar gradient.
- Use `次元還元` consistently in Japanese prose.
- State explicitly that the 5D-to-4D reduction is a zero-mode truncation on a product background with no Kaluza--Klein vector or radion fluctuations.
- Define the 5D action, determinant factor, and Hilbert stress tensor before using them.
- Show the 5D block inverse used in reducing the stress-tensor flow.
- Give an explicit derivation of the characteristic-polynomial factorization used for the null-gradient formula and the algebraic relation among stress-tensor invariants.
- State the convention map to Blair--Lahnsteiner--Obers--Yan when comparing the weak-deformation limit.
- Describe the neural-network map, loss, sampling region, rejection cuts, architecture, and the rational SVD readout.
- Remove statements about SymPy residuals from the manuscript; symbolic checks remain available in the repository.
- Explain the non-uniqueness of the algebraic representation of the flow operator modulo the invariant constraint.
- Add author/affiliation front matter, keywords, a generative-AI declaration, and a data/code availability statement.

## Checks still worth doing before English submission

- Re-read the 5D determinant/trace reduction line by line in an independent notation.
- Re-check the null-gradient branch and the sign choice of the square root directly from field variables.
- Verify the Blair et al. convention dictionary against their published equations after the final notation is fixed.
- Make the repository public and freeze a tagged version or Zenodo archive when the preprint is released.
- Perform one final citation-chain search for a later paper giving an equivalent 4D single-scalar DBI stress-tensor-only flow.
