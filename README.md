# 'Vanilla' Model Scripts

---

A collection of scripts for Bayesian modeling and analysis using the 'vanilla' model, including plotting for research papers.

## Papers:
- [Investigating the hard state of MAXI J1820 + 070: a comprehensive Bayesian approach to black hole spin and accretion properties](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529.1752D/abstract)
- Bayesian X-ray Spectroscopy of GRO J1655-40 Across Accretion States (in prep)

---

## Folder Structure
Folders are separated as follows:

**GROj1655-40_Paper:** scripts to generate tables and figures for the paper applying the 'vanilla' model to GRO J1655-40.

**MCMC:** XSPEC scripts and code for implementing Markov-chain Monte Carlo in selected examples, including setup, execution, and convergence diagnostics. 

**Nested_Sampling:** code to perform Nested Sampling analysis (including check the reliability of step sampler results) across four examples, with custom truncated priors for gamma, isotropic (sin), and lognormal distributions.

**NS_vs_MCMC:** scripts to compare the results from the Nested Sampling and Markov-chain Monte Carlo algorithms.

**System_Parameters:** code to calculate and constrain system parameters such as distance, mass, and mass accretion rates.