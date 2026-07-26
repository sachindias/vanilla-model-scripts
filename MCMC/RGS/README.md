# Notes On Markov-chain Monte Carlo - RGS

---

## What is in this Folder?

The contents of the subfolder are concerned with the running of MCMC scripts RGS data.
With the aim of getting constraints on ISM parameters to use for priors in the EPIC (+ NuSTAR) analysis.

There are two types of file in this subfolder:
1. **The models** - These are example RGS XSPEC models, for ```ISMabs * (SIMPL * KERRBB)```.
Files are named as ```<source_shorthand>_RGS_<date>.xcm```.
2. **The run scripts** - These are example XSPEC scripts for running MCMC on the RGS data & model.
Files are named as ```<source_shorthand>_RGS_run.txt```.

For both```source_shorthand``` is either ```MAXI``` for MAXI J1820+070 or ```GRO``` for GRO J1655-40.