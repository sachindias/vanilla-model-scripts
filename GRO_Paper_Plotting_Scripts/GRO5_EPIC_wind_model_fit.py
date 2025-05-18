import matplotlib.pyplot as plt
import scipy as sci
import scipy.stats as st
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import datetime

filenames = ["rev0966_MCMC_min_no_wind", "rev0966_MCMC_min"]

yellow_colour = "#C1B732"
green_colour = "#67D8B9"
pink_colour = "#DC267F"
blue_colour = "#648FFF"

colour = ["k", pink_colour]

no_spectra = 1

#SETS X WIDTH FOR SINGLE OF DOUBLE COLUMN SPREAD
double_column = False
date_now = str(datetime.datetime.now())[0:10]
minor_fontsize = 22
major_fontsize = 26
if (double_column == True):
    x_size = 24
else:
    x_size = 12
    

fig, axs = plt.subplots(3, 1, gridspec_kw={'height_ratios': [6,3,3]})
fig.set_figheight(11)
fig.set_figwidth(x_size)

for k in range(2):
    
    print("\n-----------------------------------------------")
    print("\nREADING FILE: %s\n" %filenames[k])
    
    line1 = "READ SERR 1 2"
    filename = filenames[k]
    
    XRB = open("%s.qdp" %filename, "r")
    lines = XRB.readlines()
    XRB.close()
    
    energy = []
    flux = []
    counts = []
    energy_err = []
    counts_err = []
    
    count = 0
    for line in lines:
        if (line == '%s\n' %line1 or 
            line == '@%s.pco\n' %filename or
            line == '!\n'):
            print("     IGNORED:", line)
        else:
            try:
                flux_temp = np.float(line.split()[4])
                flux.append(flux_temp)
                
                energy_temp = np.float(line.split()[0])
                energy.append(energy_temp)
                
                counts_temp = np.float(line.split()[2])
                counts.append(counts_temp)
                
                energy_err_temp = np.float(line.split()[1])
                energy_err.append(energy_err_temp)
                
                counts_err_temp = np.float(line.split()[3])
                counts_err.append(counts_err_temp)
                           
            except Exception:
                count = count + 1
     
    print ("no. Exceptions: ", count)
    energy = np.array(energy)   
    flux = np.array(flux) 
    counts = np.array(counts) 
    energy_err = np.array(energy_err)
    counts_err = np.array(counts_err)
    ratio = (counts/flux)
    
    
    E_split = energy >= 2.35
    energy_high = energy[E_split]
    flux_high = flux[E_split]
    ratio_high = ratio[E_split] 
    counts_high = counts[E_split] 
    energy_err_high = energy_err[E_split]
    counts_err_high = counts_err[E_split]
    counts_err_high_ratio = ratio_high*(counts_err_high/counts_high)
    
    E_split=np.logical_not(E_split)
    energy_low = energy[E_split]
    flux_low = flux[E_split]
    ratio_low = ratio[E_split]
    counts_low = counts[E_split] 
    energy_err_low = energy_err[E_split]
    counts_err_low = counts_err[E_split]
    counts_err_low_ratio = ratio_low*(counts_err_low/counts_low)
    

    #PLOTTING BOTH MODELS
    axs[0].hlines(y=flux_low[0], xmin = energy_low[0] - energy_err_low[0], xmax=energy_low[0] + energy_err_low[0], color = colour[k], linewidth = 1.2)
    
    for n in range(1,len(energy_low)):
        axs[0].hlines(y=flux_low[n], xmin = energy_low[n] - energy_err_low[n], xmax=energy_low[n] + energy_err_low[n], color = colour[k], linewidth = 1.2)
        axs[0].vlines(x=energy_low[n]- energy_err_low[n], ymin=flux_low[n], ymax=flux_low[n-1], color = colour[k], linewidth = 1.2)
                
    for n in range(1,len(energy_high)):
        axs[0].hlines(y=flux_high[n], xmin = energy_high[n] - energy_err_high[n], xmax=energy_high[n] + energy_err_high[n], color = colour[k], linewidth = 1.2)
        axs[0].vlines(x=energy_high[n]- energy_err_high[n], ymin=flux_high[n], ymax=flux_high[n-1], color = colour[k], linewidth = 1.2)
       
    axs[0].errorbar(energy_low, counts_low, yerr = counts_err_low, xerr = energy_err_low, color = "k", linestyle = "", marker = "", linewidth = 1)
    axs[0].errorbar(energy_high, counts_high, yerr = counts_err_high, xerr = energy_err_high, color = "k", linestyle = "", marker = "", linewidth = 1)
    
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_xlim(energy_low[0] - energy_err_low[0],9.7)
    loc = np.where(counts_low == np.max(counts_low))[0][0]
    print(0.5*(counts_high[-1] - counts_err_high[-1]))
    axs[0].set_ylim(0.5*(counts_high[-1] - counts_err_high[-1]), 1.1*counts_low[loc] + counts_err_low[loc])
    
    axs[0].set_ylabel("Counts (s$^{-1}$ keV$^{-1}$)", fontsize = major_fontsize, labelpad = 10)
    axs[0].tick_params(which='minor', length=3, direction="in", width = 1.5)
    axs[0].set_xticks ([], [])
    tick_labels = [10,100, 1000]
    axs[0].set_yticks(tick_labels)
    axs[0].set_yticklabels(['{:g}'.format(tick) for tick in tick_labels])
    axs[0].tick_params(axis="y", labelsize=minor_fontsize, direction="in", width = 2)
    
    
    #PLOTTING RATIOS FOR WITHOUT AND WITH WIND MODELS
    axs[1+k].errorbar(energy_low, ratio_low, yerr = counts_err_low_ratio, xerr = energy_err_low, color = colour[k], linestyle = "", marker = "", linewidth = 1)
    axs[1+k].errorbar(energy_high, ratio_high, yerr = counts_err_high_ratio, xerr = energy_err_high, color = colour[k], linestyle = "", marker = "", linewidth = 1)
    axs[1+k].axhline(y = 1, color = "lime", linestyle = "-", alpha = 1)
    
    axs[1+k].set_xscale("log")
    axs[1+k].set_xlim(energy_low[0] - energy_err_low[0],9.7)
    axs[1+k].set_ylim(0.76, 1.14)
    axs[1+k].set_xlabel("Energy (keV)", fontsize = major_fontsize, labelpad = 20)
    axs[1+k].set_ylabel("Ratio", fontsize = major_fontsize, labelpad = 34)
    
    axs[1+k].tick_params(which='minor', length=3, direction="in", width = 2)
    axs[1+k].tick_params(axis="x", labelsize=minor_fontsize, direction="in", width = 2)
    tick_labels = [1,2,5]
    axs[1+k].set_xticks(tick_labels)
    axs[1+k].set_xticklabels(['{:g}'.format(tick) for tick in tick_labels])
    axs[1+k].tick_params(axis="y", labelsize=minor_fontsize, direction="in", width = 2)
    axs[1+k].set_yticks([0.8, 0.9, 1.0, 1.1])
    
    plt.subplots_adjust(wspace=0, hspace=0)

#PLOTTING GABS AND SMEDGE ENERGIES
for i in range(3):

    axs[i].axvline(x = 6.68421, color = "k", linestyle = "--", alpha = 0.4)
    axs[i].axvline(x = 6.94937, color = "k", linestyle = "--", alpha = 0.4)
    axs[i].axvline(x = 7.86196, color = "k", linestyle = "--", alpha = 0.4)
    axs[i].axvline(x = 8.23026, color = "k", linestyle = "--", alpha = 0.4)
    axs[i].axvline(x = 8.44593, color = "k", linestyle = "-.", alpha = 0.4)
    
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[i].spines[axis].set_linewidth(2)

plt.savefig("GRO_J1655_40_GRO5_wind_model_%s.pdf" %date_now, bbox_inches = "tight")
plt.show()

