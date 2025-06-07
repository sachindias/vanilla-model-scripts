import matplotlib.pyplot as plt
import scipy as sci
import scipy.stats
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib as mpl
import datetime

mpl.rcParams['axes.linewidth'] = 2.5

yellow_colour = "#C1B732"
green_colour = "#67D8B9"
pink_colour = "#DC267F"
blue_colour = "#648FFF"


###################################
#EXTRACT DATA
###################################
line1 = "READ SERR 1 2"
filename = "GRO_RGS_02_08_24"

XRB = open("%s.qdp" %filename, "r")
lines = XRB.readlines()
XRB.close()

wave = []
model_flux = []
counts = []
wave_err = []
counts_err = []

count = 0
line_count = 1
for line in lines:
    if (line == '%s\n' %line1 or 
        line == '@%s.pco\n' %filename or
        line == '!\n'):
        print("     IGNORED:", line)
        line_count = line_count + 1
    else:
        try:
            if (line == "NO NO NO NO NO\n"):
                print(line, line_count)
            
            model_temp = np.float(line.split()[4])
            model_flux.append(model_temp)
            
            wave_temp = np.float(line.split()[0])
            wave.append(wave_temp)
            
            counts_temp = np.float(line.split()[2])
            counts.append(counts_temp)
            
            wave_err_temp = np.float(line.split()[1])
            wave_err.append(wave_err_temp)
            
            counts_err_temp = np.float(line.split()[3])
            counts_err.append(counts_err_temp)
            
            line_count = line_count + 1               
        except Exception:
            count = count + 1
            line_count = line_count + 1
 
print ("no. Exceptions: ", count)
wave = np.array(wave)   
model_flux = np.array(model_flux) 
counts = np.array(counts) 
wave_err = np.array(wave_err)
counts_err = np.array(counts_err)
ratio = (counts/model_flux)
counts_err_ratio = ratio*(counts_err/counts)

###################################
#PLOTTING
###################################

#SETS X WIDTH FOR SINGLE OF DOUBLE COLUMN SPREAD
double_column = False
date_now = str(datetime.datetime.now())[0:10]
minor_fontsize = 22
major_fontsize = 26
if (double_column == True):
    x_size = 24
else:
    x_size = 12
    
fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [5, 3]})
fig.set_figheight(12*(12/15))
fig.set_figwidth(x_size)

axs[0].hlines(y=model_flux[0], xmin = wave[0] - wave_err[0], xmax=wave[0] + wave_err[0], color = "k", linewidth = 1.4)

for n in range(129):
    if (n > 0):
        axs[0].hlines(y=model_flux[n], xmin = wave[n] - wave_err[n], xmax=wave[n] + wave_err[n], color = "k", linewidth = 1.4)
        axs[0].vlines(x=wave[n]- wave_err[n], ymin=model_flux[n], ymax=model_flux[n-1], color = "k", linewidth = 1.4)

for n in range(129,258):
    axs[0].hlines(y=model_flux[n], xmin = wave[n] - wave_err[n], xmax=wave[n] + wave_err[n], color = pink_colour, linewidth = 1.4)
    axs[0].vlines(x=wave[n]- wave_err[n], ymin=model_flux[n], ymax=model_flux[n-1], color = pink_colour, linewidth = 1.4)

for n in range(258,387):
    axs[0].hlines(y=model_flux[n], xmin = wave[n] - wave_err[n], xmax=wave[n] + wave_err[n], color = blue_colour, linewidth = 1.4)
    axs[0].vlines(x=wave[n]- wave_err[n], ymin=model_flux[n], ymax=model_flux[n-1], color = blue_colour, linewidth = 1.4)



axs[0].errorbar(wave[0:129], counts[0:129], yerr = counts_err[0:129], xerr = wave_err[0:129], color = "k", linestyle = "", marker = "", linewidth = 1.3)
axs[0].errorbar(wave[129:258], counts[129:258], yerr = counts_err[129:258], xerr = wave_err[129:258], color = pink_colour, linestyle = "", marker = "", linewidth = 1.3)
axs[0].errorbar(wave[258:387], counts[258:387], yerr = counts_err[258:387], xerr = wave_err[258:387], color = blue_colour, linestyle = "", marker = "", linewidth = 1.3)

axs[0].set_yscale("log")
axs[0].set_yticks([1e-3, 1e-2, 1e-1])
axs[0].set_yticklabels([r'1$\times10^{-3}$', r'1$\times10^{-2}$', r'1$\times10^{-1}$'])
axs[0].set_xlim(11.02,24)
axs[0].tick_params(axis="y", labelsize=18, direction="in", length = 5, width = 2)
axs[0].tick_params(which='minor', length=3, direction="in", width = 2)
axs[0].set_ylim(0.0005, 0.2)
axs[0].set_xticks ([], [])

axs[0].set_ylabel("Counts (s$^{-1}$ Å$^{-1}$ cm$^{-2}$)", fontsize = major_fontsize, labelpad=15)
axs[0].tick_params(axis="y", labelsize=minor_fontsize)


axs[1].errorbar(wave[0:129], ratio[0:129], yerr = counts_err_ratio[0:129], xerr = wave_err[0:129], color = "k", linestyle = "", marker = "", linewidth = 1.3)
axs[1].errorbar(wave[129:258], ratio[129:258], yerr = counts_err_ratio[129:258], xerr = wave_err[129:258], color = pink_colour, linestyle = "", marker = "", linewidth = 1.3)
axs[1].errorbar(wave[258:387], ratio[258:387], yerr = counts_err_ratio[258:387], xerr = wave_err[258:387], color = blue_colour, linestyle = "", marker = "", linewidth = 1.3)
axs[1].axhline(y = 1, color = "lime", linestyle = "-", alpha = 1)

axs[1].set_xlim(11,24)
axs[1].set_ylim(0.501, 1.49)
axs[1].set_xlabel(r"Wavelength (Å)", fontsize = major_fontsize, labelpad = 20)
axs[1].set_ylabel("Ratio", fontsize = major_fontsize, labelpad=66)

axs[1].set_yticks([0.6,1, 1.4])
axs[1].set_yticklabels([0.6,1, 1.4])
axs[1].tick_params(axis="y", labelsize=minor_fontsize, direction="in", length=5, width = 2)
axs[1].tick_params(axis="x", labelsize=minor_fontsize, direction="in", length=5, width = 2)
axs[1].tick_params(which='minor', length=3, direction="in", width = 2)
axs[1].xaxis.set_minor_locator(MultipleLocator(1))
axs[1].yaxis.set_minor_locator(MultipleLocator(0.2))

for n in range(2):
    axs[n].axvspan(22.7,23.2, color = "gray", alpha = 0.4)

plt.subplots_adjust(wspace=0, hspace=0)
for axis in ['top', 'bottom', 'left', 'right']:
    axs[0].spines[axis].set_linewidth(2)
    axs[1].spines[axis].set_linewidth(2)
plt.savefig("Plots/GRO_RGS_SPECTRA_%s.pdf" %date_now, bbox_inches = "tight")
plt.show()