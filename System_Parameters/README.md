# Notes On System Parameters

## What is in this Folder?

- **distance_priors.py** - For constructing a Gamma prior on distance for the MCMC and NS analyses.
Distances are sampled for each of the systems using parallax measurements and the method defined in Bailer-Jones (2015). 

- **mass_priors.py** - For constructing a shifted log-normal prior on mass for the MCMC and NS analyses. 
Black hole masses are sampled for each of the systems using the method defined in Farr et al. (2011), 
i.e. using priors on the mass function, mass ratio and inclination.
- **mass_accretion_rate_limits.py** - Calculates upper and lower limits on the mass accretion rate for each of the systems,
using the 10th and 90th centiles from the distribution of black hole masses (see above).

Justification behind values used in the code can be found in: 
  - [Dias et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529.1752D/abstract) - for MAXI J1820=070
  - Dias et al. (in prep) - for GRO J1655-40
  - [Dias (2024)](https://figshare.le.ac.uk/articles/thesis/A_Bayesian_Approach_to_Black_Hole_Spin_and_Accretion_Properties/26251382?file=47583947) - for all sources, including the above two.