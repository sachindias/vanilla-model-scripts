import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import datetime
import os
from scipy.optimize import minimize
import celerite
from celerite import terms

# TIMER
print("\n--------------------------------------------")
start_timer = datetime.datetime.now()
print(start_timer)
print("--------------------------------------------")


# -----------------------------------------
# FUNCTIONS

def extender(new_x, new, N):
    # FIT A POWER LAW TO THE LAST TWO POINTS AND EXTRAPOLATE TO new_x
    m = np.log10(new[-1] / new[-2]) / np.log10(N[-1] / N[-2])
    b = new[-1] / (N[-1] ** m)
    y = b * new_x ** m

    return y


def next_pow_two(n):
    # ROUND n UP TO THE NEXT POWER OF TWO (FOR EFFICIENT FFT LENGTH)
    i = 1
    while i < n:
        i = i << 1
    return i


def autocorr_func_1d(x, norm=True):
    # COMPUTE THE 1D AUTOCORRELATION FUNCTION OF A SINGLE CHAIN VIA FFT
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # COMPUTE THE FFT AND THEN (FROM THAT) THE AUTO-CORRELATION FUNCTION
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # OPTIONALLY NORMALIZE
    if norm:
        acf /= acf[0]

    return acf


# AUTOMATED WINDOWING PROCEDURE FOLLOWING SOKAL (1989)
def auto_window(taus, c):
    # FIND THE SMALLEST WINDOW SIZE WHERE M < c * tau(M), OR USE THE FULL LENGTH IF NONE FOUND
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1


# FOLLOWING THE SUGGESTION FROM GOODMAN & WEARE (2010)
def autocorr_gw2010(y, c=5.0):
    # AVERAGE CHAINS FIRST, THEN COMPUTE AUTOCORRELATION AND INTEGRATE TO GET TAU
    f = autocorr_func_1d(np.mean(y, axis=0))
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]


def autocorr_new(y, c=5.0):
    # COMPUTE AUTOCORRELATION PER CHAIN THEN AVERAGE (IMPROVED ESTIMATOR VS GW2010)
    f = np.zeros(y.shape[1])
    for yy in y:
        f += autocorr_func_1d(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]


def autocorr_ml(y, thin=1, c=5.0):
    # COMPUTE THE INITIAL ESTIMATE OF TAU USING THE STANDARD METHOD
    init = autocorr_new(y, c=c)
    z = y[:, ::thin]
    N = z.shape[1]

    # BUILD THE GP MODEL
    tau = max(1.0, init / thin)
    kernel = terms.RealTerm(
        np.log(0.9 * np.var(z)),
        -np.log(tau),
        bounds=[(-5.0, 5.0), (-np.log(N), 0.0)],
    )
    kernel += terms.RealTerm(
        np.log(0.1 * np.var(z)),
        -np.log(0.5 * tau),
        bounds=[(-5.0, 5.0), (-np.log(N), 0.0)],
    )
    gp = celerite.GP(kernel, mean=np.mean(z))
    gp.compute(np.arange(z.shape[1]))

    # DEFINE THE OBJECTIVE
    def nll(p):
        # UPDATE THE GP MODEL
        gp.set_parameter_vector(p)

        # LOOP OVER THE CHAINS AND COMPUTE LIKELIHOODS
        v, g = zip(*(gp.grad_log_likelihood(z0, quiet=True) for z0 in z))

        # COMBINE THE DATASETS
        return -np.sum(v), -np.sum(g, axis=0)

    # OPTIMIZE THE MODEL
    p0 = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()
    soln = minimize(nll, p0, jac=True, bounds=bounds)
    gp.set_parameter_vector(soln.x)

    # COMPUTE THE MAXIMUM LIKELIHOOD TAU
    a, c = kernel.coefficients[:2]
    tau = thin * 2 * np.sum(a / c) / np.sum(a)
    return tau


# -----------------------------------------
# MAIN SCRIPT - DATA SORTING & CALCULATIONS
# THE EXAMPLE HERE IS SPIN FOR GRO1 (rev0956)

# INITIAL SETTINGS
param = "a__37"
n_walkers = 200
save_filename = "GRO1_AutoCorrelation_Spin"

# FINDS ALL FILES WITH DATA FOR A SPECIFIC MCMC RUN
# THIS WILL BE DIFFERENT FOR EACH CHAIN
# AS IT MUST INCLUDE EVERY BATCH
section_time = datetime.datetime.now()
FILENAME = "rev0956_27_9_23_afree_2M"

hdus = [
    '',
    "_ext",
]

ext_cond = 0
ext_count = 1

while (ext_cond < 1):
    try:
        ext_count = ext_count + 1

        if (os.path.getsize("%s_ext%s.fits" % (FILENAME, ext_count)) > 0):
            print("%s_ext%s.fits" % (FILENAME, ext_count))
            hdus.append("_ext%s" % ext_count)
    except Exception:
        ext_cond = 2

print("\nDATA SORTED\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))
print("--------------------------------------------")

# EXTRACTS ALL VALUES OF THE PARAMETER FROM ALL BATCHES AND PUTS INTO 1D ARRAY
section_time = datetime.datetime.now()
min_iteration = 1
hdu_param_values = []

for hdu in hdus:
    hdul = fits.open('%s%s.fits' % (FILENAME, hdu))
    cols = hdul[1].columns
    data = hdul[1].data

    hdu_param_values.append(data[param])
    hdul.close()

flat_param_values = np.concatenate(hdu_param_values, axis=0)

print("DATA EXTRACTED\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))
print("--------------------------------------------")

# RESHAPES THE 1D ARRAY INTO 2D ARRAY ORGANISED BY ITERATION& WALKER
section_time = datetime.datetime.now()
orgnaised_param_values = []

for n in range(len(flat_param_values) // n_walkers):
    chain_rev3623_iteration = []
    for i in range(n_walkers):
        chain_rev3623_iteration.append(flat_param_values[n * n_walkers + i])
    orgnaised_param_values.append(np.array(chain_rev3623_iteration))

transposed_param_values = np.array(orgnaised_param_values).T  # shape: (walkers, iterations)

print("DATA RESHAPED\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))
print("--------------------------------------------")

# COMPUTES THE ESTIMATORS FOR A FEW DIFFERENT CHAIN LENGTHS
section_time = datetime.datetime.now()
N = np.exp(np.linspace(np.log(100), np.log(transposed_param_values.shape[1]), 10)).astype(int)
gw2010 = np.empty(len(N))
new = np.empty(len(N))
for i, n in enumerate(N):
    gw2010[i] = autocorr_gw2010(transposed_param_values[:, :n])
    new[i] = autocorr_new(transposed_param_values[:, :n])

print("GW & NEW METHODS DONE\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))
print("--------------------------------------------")

# COMPUTE ML AUTOCORRELATION ESTIMATES FOR A RANGE OF CHAIN LENGTHS (SKIPPING THE FIRST N)
section_time = datetime.datetime.now()

try:
    ml = np.empty(len(N))
    ml[:] = np.nan
    for j, n in enumerate(N[1:10]):
        i = j + 1
        thin = max(1, int(0.05 * new[i]))
        ml[i] = autocorr_ml(transposed_param_values[:, :n], thin=thin)

    print("ML METHOD DONE\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
    datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))
except Exception:
    print("ML METHOD FAILED\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
    datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))
print("--------------------------------------------")

# EXTRAPOLATE THE 'NEW' AND 'ML' ESTIMATORS OUT TO A MUCH LARGER SAMPLE SIZE
section_time = datetime.datetime.now()
print("EXTENDING: %s" % param)

# EXTRAPOLATE 'NEW' ESTIMATOR TO new_x SAMPLES
new_x = 1e7
new_y = extender(new_x, new, N)

N1 = np.array(list(N) + [new_x])
new1 = np.array(list(new) + [new_y])

# BUILD DENSE ARRAY BETWEEN LAST COMPUTED N AND new_x TO SEARCH FOR INTERSECTION WITH N/50
cross_array = np.logspace(np.log10(N[-1]), np.log10(new_x), 10000)
cross_array_50 = cross_array / 50

# EVALUATE THE EXTRAPOLATED 'NEW' CURVE ACROSS THAT DENSE RANGE
cross_y_array = []
for n in range(len(cross_array)):
    cross_y = extender(cross_array[n], new, N)
    cross_y_array.append(cross_y)

cross_y_array = np.array(cross_y_array)
minus_array = np.subtract(cross_y_array, cross_array_50)

# FIND WHERE THE EXTRAPOLATED 'NEW' CURVE CROSSES N/50 (CLOSEST POINT TO ZERO DIFFERENCE)
if (len(np.where(minus_array < 0)[0]) != 0):
    mag_array = abs(minus_array)
    mag_array_sorted = np.sort(mag_array)
    mag_low = mag_array_sorted[0]
    min_pos = np.where(mag_array == mag_low)[0][0]
    print("EXTENDED 'new' : MAGNITUDE: %.4f ; TAU ESTIMATE: %.0f" % (mag_array[min_pos], cross_y_array[min_pos]))
else:
    print("EXTENDED 'new' : MAGNITUDE: N/A ; TAU ESTIMATE: NO INTERSECTION")
    min_pos = len(cross_array) - 1

# REPEAT THE SAME EXTRAPOLATION + INTERSECTION SEARCH FOR THE 'ML' ESTIMATOR
try:
    new_y_ml = extender(new_x, ml, N)
    new1_ml = np.array(list(ml) + [new_y_ml])

    cross_array = np.logspace(np.log10(N[-1]), np.log10(new_x), 10000)
    cross_array_50 = cross_array / 50

    cross_y_ml_array = []
    for n in range(len(cross_array)):
        cross_y_ml = extender(cross_array[n], ml, N)
        cross_y_ml_array.append(cross_y_ml)

    cross_y_ml_array = np.array(cross_y_ml_array)
    minus_array_ml = np.subtract(cross_y_ml_array, cross_array_50)

    if (len(np.where(minus_array_ml < 0)[0]) != 0):
        mag_array_ml = abs(minus_array_ml)
        mag_array_ml_sorted = np.sort(mag_array_ml)
        mag_low_ml = mag_array_ml_sorted[0]
        min_pos_ml = np.where(mag_array_ml == mag_low_ml)[0][0]
        print("EXTENDED 'ml' : MAGNITUDE: %.4f ; TAU ESTIMATE: %.0f" % (
        mag_array_ml[min_pos_ml], cross_y_ml_array[min_pos_ml]))
    else:
        print("EXTENDED 'ml' : MAGNITUDE: N/A ; TAU ESTIMATE: NO INTERSECTION")
        min_pos_ml = 0
except Exception:
    print("ML EXTEND FAILED")

print("EXTENSION DONE\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))


# -----------------------------------------
# PLOTTING
section_time = datetime.datetime.now()

# SET UP FIGURE
plt.figure(figsize=(10, 10))

# PLOT THE COMPARISONS (G&W 2010 vs 'NEW' METHOD, INCLUDING EXTENDED 'NEW' RANGE)
plt.loglog(N, gw2010, "o-", color="tab:blue", label="G&W 2010", linewidth=2, markersize=8)
plt.loglog(N, new, "o-", color="tab:orange", label="New", linewidth=2, markersize=8)
plt.loglog(N1[9:], new1[9:], "o--", color="tab:orange", label="New (extended)", linewidth=2, markersize=8)
plt.legend(fontsize=14)
plt.tick_params(labelsize=14)

# ATTEMPT TO PLOT ML ESTIMATOR (NOT ALWAYS AVAILABLE FOR EVERY PARAMETER)
try:
    plt.loglog(N, ml, "o-", color="tab:green", label="ML", linewidth=2, markersize=8)
    plt.loglog(N1[9:], new1_ml[9:], "o--", color="tab:green", label="ML (extended)", linewidth=2, markersize=8)
    ml_plotting = 1
except Exception:
    ml_plotting = 0
    print("ML PLOT FAILED")

# ADD REFERENCE LINE (tau = N/50) WITHOUT LETTING IT RESCALE THE Y AXIS
ylim = plt.gca().get_ylim()
plt.plot(N1, N1 / 50.0, "--k", label=r"$\tau = N/50$", linewidth=2)
plt.ylim(ylim)

# LABEL X AXIS WITH PARAMETER NAME
plt.xlabel(r"number of samples, $N$: %s" % param, fontsize=16)

# ESTIMATE TAU FROM 'NEW' METHOD'S INTERSECTION WITH THE REFERENCE LINE
tau_new = int(cross_y_array[min_pos])
if (len(np.where(minus_array < 0)[0]) == 0):
    tau_new = "NO INT"

# ESTIMATE TAU FROM ML METHOD'S INTERSECTION (IF ML WAS SUCCESSFULLY PLOTTED)
# THEN LABEL Y AXIS
if (ml_plotting == 1):
    try:
        tau_ml = int(cross_y_ml_array[min_pos_ml])
    except Exception:
        tau_ml = "nan"
    if (len(np.where(minus_array_ml < 0)[0]) != 0):
        tau_new = "NO INT"

    plt.ylabel(r"$\tau$ estimates:" + "\n" + r"$\tau$ 'new' X: %s ; $\tau$ 'ml' X: %s" % (tau_new, tau_ml), fontsize=16,
               labelpad=15)
else:
    plt.ylabel(r"$\tau$ estimates:" + "\n" + r"$\tau$ 'new' X estimate: %.0f" % cross_y_array[min_pos], fontsize=16,
               labelpad=15)

# SAVE FIGURE
plt.savefig("%s.png" % save_filename, bbox_inches="tight")

print("--------------------------------------------")
print("PLOTTING DONE\nSECTION TIME: %s\nTOTAL RUN TIME: %s" % (
datetime.datetime.now() - section_time, datetime.datetime.now() - start_timer))

# -----------------------------------------
# END

print("--------------------------------------------")
print("TOTAL TIME TAKEN FOR %s: %s" % (param, datetime.datetime.now() - start_timer))