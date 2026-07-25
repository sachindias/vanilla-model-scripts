# Notes on the GROJ1655-40 Paper Folder

## What is in this Folder?

This folder is concerned with scripts used for the GRO J1655-40 paper (currently in prep).
In order of appearance in the paper, the scripts are as follows:

1. **RGS_model_fit.py** - Plots RGS counts and ratios for GRO2, GRO3 and GRO4.
2. **GRO_EPIC_wind_model_fit.py** - Plots the EPIC counts and ratios for GRO5, with and without model components to account for the wind.
3. **MCMC_NS_GRO_violin_plots.py** - Violin plots for MCMC vs NS for spin, inclination, mass and the (base-10) logarithm of the mass accretion rate.
4. **spin_fcol_corner_plot.py** - MCMC Corner plot for spin vs fcol for GRO2-5, demonstrating the $R_{in}$ and $f^2_{col}$ relationship. 
5. **my_corner.py** - My own corner plotting script, used in the above script.
6. **GRO_best_fits.py** - Plots the EPIC counts and ratios for each of the observations. 
For GRO1 an additional plot is included to show the POWERLAW + GAUSSIAN model over 4 - 8 keV.
7. **GRO_automated_table_MCMC.py** - Creates the LaTeX table for the MCMC results in the appendix.
8. **GRO_automated_table_NS.py** - Creates the LaTeX table for the Nested Sampling results in the appendix.