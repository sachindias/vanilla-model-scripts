## Notes On Nested Sampling

### What is in this Folder?

**Custom_Priors:** These are probability distributions (not previously included in PyXspec) for priors needed in the Nested Sampling analysis. 

**NS_Run_Scripts:** Here are some example scripts needed to run the Nested Sampling analysis. 

**Comparison Scripts:** Due to unreliable results with the smaller step sizes, these scripts give examples for how to check for consistency in both parameter values and the evidence.

### Running NS Scripts

To run the scripts in **NS_Run_Scripts**, you will need to follow these steps:

1. initialise **heasoft**.
2. activate **conda** environment.
3. `cd` into the correct directory.
4. Set **OpenMP** to use a single thread.
5. Run the script using `mpiexec` for parallel comupting.

e.g.
```
. /home/heasoft-6.30.1/x86_64-pc-linux-gnu-libc2.34/headas-init.sh
source ~/miniconda3/bin/activate nested

cd /to/your/Nested_Sampling/NS_Run_Scripts

export OMP_NUM_THREADS=1
mpiexec -n 4 python EXAMPLE_SCRIPT.py
```