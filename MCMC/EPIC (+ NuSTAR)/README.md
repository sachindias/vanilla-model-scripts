# Notes On Markov-chain Monte Carlo - EPIC (+ NuSTAR)

---

## What is in this Folder?

Scripts are used for monitoring, determining convergence on and extracting data from the MCMC chain:

- **trace_plots.py** - For visualising the evolution of chains (the trace) per parameter. 
Useful in determining if chains are clustering, hitting the bounds of the model, or approaching convergence.
There are two options:
  1. Plotting each chain as distinct.
  2. Plotting the mean and standard deviation of the bulk of chains.
- **per_param_autocorrelation.py** - For calculating (or extrapolating to estimate) the autocorrelation time of a single parameter.
Due to the extreme length of chains (all of which is needed for the calculation), 
this may need to be run on a server with sufficient memory to hold the entire chain in RAM.
- **mean_v_width.py** - For calculating the "measure of the change in the mean relative 
to the width of the distribution" between different points ($i$ and $j$) in the chain.
The formula:

$$
\sqrt{2} \frac{\mu_i - \mu_j}{(S_i^2 + S_j^2)^{1/2}}
$$

comes from [Dias et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529.1752D/abstract) and is a useful diagnostic for determining convergence.
- **min_fit_stat_param_values.py** - For extracting the parameter values from the chain corresponding to the minimum fit statistic (chi-squared) value.
- **average_fit_stat.py** - For calculating the average fit statistic (chi-squared) value across the final sample.
Useful when comparing to the typical set.

### What is in the XSPEC Folder
The contents of the subfolder are concerned with the actual running of MCMC scripts for EPIC (and in the case of MAXI J1820+070, NuSTAR) data.

There are python scripts and a further three additional subfolder (GRO5, MAXI_Obs1, MAXI_Obs2) each containing:
- An example of the model used - `<obs>_<date>_<spin>.xcm`, 
with obs being the XSPEC observation name, e.g. rev0966 for GRO5.
- An initialisation script for setting up the runs in XSPEC - `<obs>_<spin>_initialisation.txt`, 
with obs being the `MAXI_Obs<X>` for MAXI J1820+070 and the `GRO<X>` for GRO J1655-40,
where `X` is the observation number in either case.
- A script to run the MCMC analysis in XSPEC - `<obs>_<spin>_run.txt`, following the same naming convention as the previous point.

The scripts in this subfolder are as follows:
- **initial_walker_selection.py** - For ensuring that the initial walkers are not outside the bounds of the model.
- **additional_batch_selection.py** - For setting the next batch to start from the end of the previous batch,
i.e. replacing the initial walker positions of the new batch with the final walker positions from the last.
- **remove_clustering.py** - For removing clustering behaviour in the walkers during the warm-up period.

## FITS File Structure

### Explaining the names for each FITS files

For those who are able to access the data, it consists of several ~ 400MB files. Due to the complexities during warm-up and having to guide the walkers (to prevent them getting stuck in local minima) files have complicated names. These are detailed here for ease of understanding.

**The naming convention isn't perfect in every case, but should give a rough idea as to what the filenames refer to.**

File names are according to the following format:
`<XMM_revolution>_<date>_<spin_value>_<chain_length>_<additional_information>.fits`

`<XMM_revolution>` refers to the XMM resolution the observation was taken from.
`<date>` refers to when the first run of the chain began.
`<spin_value>` refers to the spin value for that run (a = 0, a = 0.998 or a is a free parameter).
`<chain_length>` refers to the length of the chain as set in XSPEC (10,000 or 2,000,000).
`<additional_information>` refers to additional information related to the order of runs from continuations (or extensions) or due to long and complicated warm-up periods.

### XMM_revolution:
###### MAXI J1820+070:
- rev3531: _Obs1_
- rev3533: _Obs2_
- rev3623: _Obs3_
###### GRO J1655-40:
- rev0956: _GRO1_
- rev0964a: _GRO2_
- rev0964b: _GRO3_
- rev0965: _GRO4_
- rev0966: _GRO5_
- rev0970: _GRO6_

### spin_value:
- afree: _spin is a free parameter_
- a0: _spin is fixed at a=0_
- amax: _spin is fixed at a=0.998_

### chain_length:
- 10K: _chain length of 10,000. These 10K were primarily run as a test before the main bulk._ 
- 2M: _chain length of 2,000,000_

### additional_information:
- P_: _this run was started from values in the best 25% chi-sq values of the previous run (to avoid clustering behaviour). For each of these P is followed by a number, this is just a rough guide for how many times this has occured. P2 meaning once, P3 meaning twice, etc..._
- ext_: _this run was started from the final values of the previous run, effectively extending it._

-----
### Filenames:

#### Obs1:

###### afree

rev3531_1_9_22_afree_10K.fits<br>
rev3531_1_9_22_afree_2M.fits<br>
rev3531_1_9_22_afree_2M_P2.fits<br>
rev3531_1_9_22_afree_2M_P2_ext.fits<br>
rev3531_1_9_22_afree_2M_P2_ext2.fits<br>
↓<br>
rev3531_1_9_22_afree_2M_P2_ext36.fits<br>

###### a0

rev3531_1_9_22_a0_10K.fits<br>
rev3531_1_9_22_a0_2M.fits<br>
rev3531_1_9_22_a0_2M_P2.fits<br>
rev3531_1_9_22_a0_2M_P2_ext.fits<br>
rev3531_1_9_22_a0_2M_P2_ext2.fits<br>
↓<br>
rev3531_1_9_22_a0_2M_P2_ext36.fits<br>

###### amax

rev3531_1_9_22_amax_10K.fits<br>
rev3531_1_9_22_amax_2M.fits<br>
rev3531_1_9_22_amax_2M_P2.fits<br>
rev3531_1_9_22_amax_2M_P2_ext.fits<br>
rev3531_1_9_22_amax_2M_P2_ext2.fits<br>
↓<br>
rev3531_1_9_22_amax_2M_P2_ext38.fits<br>


#### Obs2:

###### afree

rev3533_1_9_22_afree_10K.fits<br>
rev3533_1_9_22_afree_2M.fits<br>
rev3533_1_9_22_afree_2M_P2.fits<br>
rev3533_1_9_22_afree_2M_P2_ext.fits<br>
rev3533_1_9_22_afree_2M_P2_ext2.fits<br>
↓<br>
rev3533_1_9_22_afree_2M_P2_ext31.fits<br>

###### a0

rev3533_1_9_22_a0_10K.fits<br>
rev3533_1_9_22_a0_2M.fits<br>
rev3533_1_9_22_a0_2M_P2.fits<br>
rev3533_1_9_22_a0_2M_P2_ext.fits<br>
rev3533_1_9_22_a0_2M_P2_ext2.fits<br>
↓<br>
rev3533_1_9_22_a0_2M_P2_ext24.fits<br>

###### amax

rev3533_1_9_22_amax_10K.fits<br>
rev3533_1_9_22_amax_2M.fits<br>
rev3533_1_9_22_amax_2M_P2.fits<br>
rev3533_1_9_22_amax_2M_P2_ext.fits<br>
rev3533_1_9_22_amax_2M_P2_ext2.fits<br>
↓<br>
rev3533_1_9_22_amax_2M_P2_ext23.fits<br>

#### Obs3:

###### afree

rev3623_1_9_22_afree_10K.fits<br>
rev3623_1_9_22_afree_2M.fits<br>
rev3623_1_9_22_afree_2M_ext.fits<br>
rev3623_1_9_22_afree_2M_ext2.fits<br>
↓<br>
rev3623_1_9_22_afree_2M_ext23.fits<br>

###### a0

rev3623_1_9_22_a0_10K.fits<br>
rev3623_1_9_22_a0_2M.fits<br>
rev3623_1_9_22_a0_2M_P2.fits<br>
rev3623_1_9_22_a0_2M_P3.fits<br>
rev3623_1_9_22_a0_2M_P3_ext.fits<br>
rev3623_1_9_22_a0_2M_P3_ext2.fits<br>
↓<br>
rev3623_1_9_22_a0_2M_P3_ext25.fits<br>

###### amax

rev3623_1_9_22_amax_10K.fits<br>
rev3623_1_9_22_amax_2M.fits<br>
rev3623_1_9_22_amax_2M_P2.fits<br>
rev3623_1_9_22_amax_2M_P2_ext.fits<br>
rev3623_1_9_22_amax_2M_P2_ext2.fits<br>
↓<br>
rev3623_1_9_22_amax_2M_P2_ext25.fits<br>

#### GRO1:

rev0956_27_9_23_afree_10K.fits<br>
rev0956_27_9_23_afree_2M.fits<br>
rev0956_27_9_23_afree_2M_ext.fits<br>
rev0956_27_9_23_afree_2M_ext2.fits<br>
↓<br>
rev0956_27_9_23_afree_2M_ext162.fits<br>

#### GRO2:

rev0964a_15_9_23_afree_10K.fits<br>
rev0964a_15_9_23_afree_2M.fits<br>
rev0964a_15_9_23_afree_2M_P2.fits<br>
rev0964a_15_9_23_afree_2M_P2_ext.fits<br>
rev0964a_15_9_23_afree_2M_P2_ext2.fits<br>
↓<br>
rev0964a_15_9_23_afree_2M_P2_ext32.fits<br>
rev0964a_15_9_23_afree_2M_P2_ext32_P2.fits<br>
rev0964a_15_9_23_afree_2M_P2_ext32_P2_ext.fits<br>
rev0964a_15_9_23_afree_2M_P2_ext32_P2_ext2.fits<br>
↓<br>
rev0964a_15_9_23_afree_2M_P2_ext32_P2_ext119.fits<br>

#### GRO3:

rev0964b_27_9_23_afree_10K.fits<br>
rev0964b_27_9_23_afree_2M.fits<br>
rev0964b_27_9_23_afree_2M_P2.fits<br>
rev0964b_27_9_23_afree_2M_P2_ext.fits<br>
rev0964b_27_9_23_afree_2M_P3.fits<br>
rev0964b_27_9_23_afree_2M_P3_ext.fits<br>
rev0964b_27_9_23_afree_2M_P3_ext2.fits<br>
↓<br>
rev0964b_27_9_23_afree_2M_P3_ext162.fits<br>

#### GRO4:

rev0965_27_9_23_afree_10K.fits<br>
rev0965_27_9_23_afree_2M.fits<br>
rev0965_27_9_23_afree_2M_ext.fits<br>
rev0965_27_9_23_afree_2M_ext2.fits<br>
↓<br>
rev0965_27_9_23_afree_2M_ext19.fits<br>
rev0965_27_9_23_afree_2M_ext19_P2.fits<br>
rev0965_27_9_23_afree_2M_ext19_P2_ext.fits<br>
rev0965_27_9_23_afree_2M_ext19_P2_ext2.fits<br>
↓<br>
rev0965_27_9_23_afree_2M_ext19_P2_ext136.fits<br>

#### GRO5:

rev0966_27_9_23_afree_10K.fits<br>
rev0966_27_9_23_afree_2M.fits<br>
rev0966_27_9_23_afree_2M_P2.fits<br>
rev0966_27_9_23_afree_2M_P2_ext.fits<br>
rev0966_27_9_23_afree_2M_P3.fits<br>
rev0966_27_9_23_afree_2M_P3_ext.fits<br>
rev0966_27_9_23_afree_2M_P3_ext2.fits<br>
↓<br>
rev0966_27_9_23_afree_2M_P3_ext161.fits<br>

#### GRO6:

rev0970_27_9_23_afree_10K.fits<br>
rev0970_27_9_23_afree_2M.fits<br>
rev0970_27_9_23_afree_2M_P2.fits<br>
rev0970_27_9_23_afree_2M_P2_ext.fits<br>
rev0970_27_9_23_afree_2M_P3.fits<br>
rev0970_27_9_23_afree_2M_P3_ext.fits<br>
rev0970_27_9_23_afree_2M_P3_ext2.fits<br>
↓<br>
rev0970_27_9_23_afree_2M_P3_ext42.fits<br>
rev0970_27_9_23_afree_2M_P3_ext42_P2.fits<br>
rev0970_27_9_23_afree_2M_P3_ext42_P2_ext.fits<br>
rev0970_27_9_23_afree_2M_P3_ext42_P2_ext2.fits<br>
↓<br>
rev0970_27_9_23_afree_2M_P3_ext42_P2_ext103.fits<br>