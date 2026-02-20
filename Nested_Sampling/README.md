## Notes On Nested Sampling

### Running Scripts

1. initialise **heasoft**.
2. activate **conda** environment.
3. `cd` into the correct directory.
4. Set **OpenMP** to use a single thread.
5. Run the script using `mpiexec` for parallel comupting.

e.g.
```
. /home/heasoft-6.30.1/x86_64-pc-linux-gnu-libc2.34/headas-init.sh
source ~/miniconda3/bin/activate nested

cd post_processing/nested_sampling_scripts

export OMP_NUM_THREADS=1
mpiexec -n 4 python EXAMPLE_SCRIPT.py
```