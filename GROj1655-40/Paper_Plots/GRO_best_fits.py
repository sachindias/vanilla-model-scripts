import matplotlib.pyplot as plt
import scipy as sci
import scipy.stats
import numpy as np
import matplotlib as mpl
import datetime
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

yellow_colour = "#C1B732"
green_colour = "#67D8B9"
pink_colour = "#DC267F"
blue_colour = "#648FFF"

colour = ["k", pink_colour]

###################################
#Set Up
###################################
file_chooser = 3 #starting at 0

filenames = [
            "rev0956_MCMC_min",
             "rev0964a_MCMC_min",
             "rev0964b_MCMC_min",
             "rev0965_MCMC_min",
             "rev0966_MCMC_min",
             "rev0970_MCMC_min"
             ]

labelpady_1_list = [30] + [10] * 5

labelpady_2_list = [13] + [20] * 5

ylim_low  = [0.61, 0.87, 0.87, 0.92, 0.89, 0.89]
ylim_high = [1.58, 1.17, 1.17, 1.12, 1.12, 1.21]

yellow_colour = "#C1B732"
green_colour = "#67D8B9"
pink_colour = "#DC267F"
blue_colour = "#648FFF"

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
    
fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3,2]})
fig.set_figheight(12*(12/15))
fig.set_figwidth(x_size)
mpl.rcParams['axes.linewidth'] = 2

###################################
#Extracting Data
###################################
for k in range(1):
    
    print("\n-----------------------------------------------")
    print("\nREADING FILE: %s\n" %filenames[k])
    
    line1 = "READ SERR 1 2"
    filename = filenames[file_chooser]
    
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
    
    if (no_spectra == 1):
    
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
        
        ###################################
        #PLOTTING
        ###################################
        
        if (k == 0):
            
            axs[0].hlines(y=flux[0], xmin = energy[0] - energy_err[0], xmax=energy[0] + energy_err[0], color = "k", linewidth = 1.3)
            
            for n in range(1,len(energy_low)):
                axs[0].hlines(y=flux_low[n], xmin = energy_low[n] - energy_err_low[n], xmax=energy_low[n] + energy_err_low[n], color = "k", linewidth = 1.2)
                axs[0].vlines(x=energy_low[n]- energy_err_low[n], ymin=flux_low[n], ymax=flux_low[n-1], color = "k", linewidth = 1.2)
                        
            for n in range(1,len(energy_high)):
                axs[0].hlines(y=flux_high[n], xmin = energy_high[n] - energy_err_high[n], xmax=energy_high[n] + energy_err_high[n], color = "k", linewidth = 1.2)
                axs[0].vlines(x=energy_high[n]- energy_err_high[n], ymin=flux_high[n], ymax=flux_high[n-1], color = "k", linewidth = 1.2)
              
            
            axs[0].errorbar(energy_low, counts_low, yerr = counts_err_low, xerr = energy_err_low, color = "k", linestyle = "", marker = "", linewidth = 1)
            axs[0].errorbar(energy_high, counts_high, yerr = counts_err_high, xerr = energy_err_high, color = "k", linestyle = "", marker = "", linewidth = 1)
            
            
            axs[0].set_xscale("log")
            axs[0].set_yscale("log")
            axs[0].set_yticks([0.1,0.2,0.5,2,5,10,20,50, 100], [0.1,0.2,0.5,2,5,10,20,50, 100])
            axs[0].set_xlim(energy_low[0] -energy_err_low[0] ,energy_high[-1] +energy_err_high[-1] )
            axs[0].set_xticks ([], [])
            
            labelpady_1 = labelpady_1_list[file_chooser]
            
            labelpady_2 = labelpady_2_list[file_chooser]
            
            labelpadx = 20
            
            tick_labels = [1,10,100, 1000]
            axs[0].set_yticks(tick_labels)
            axs[0].set_yticklabels(['{:g}'.format(tick) for tick in tick_labels])
            loc_high = np.where(np.array(counts_low) == max(counts_low))[0][0]
            axs[0].set_ylim(0.5*(counts_high[-1] - counts_err_high[-1]), counts_low[loc_high] + counts_err_low[loc_high])
            
            axs[0].set_ylabel("Counts (s$^{-1}$ keV$^{-1}$)", fontsize = major_fontsize, labelpad = labelpady_1)
            axs[0].tick_params(axis="y", labelsize=minor_fontsize)
        
        axs[k+1].errorbar(energy_low, ratio_low, yerr = counts_err_low_ratio, xerr = energy_err_low, color = "k", linestyle = "", marker = "", linewidth = 1)
        axs[k+1].errorbar(energy_high, ratio_high, yerr = counts_err_high_ratio, xerr = energy_err_high, color = "k", linestyle = "", marker = "", linewidth = 1)
        
        axs[k+1].set_xscale("log")
        axs[k+1].set_ylim(ylim_low[file_chooser], ylim_high[file_chooser]) #rev0956

        axs[k+1].set_xlim(energy_low[0] -energy_err_low[0] ,energy_high[-1] +energy_err_high[-1] )
        
        axs[k+1].set_xlim(energy_low[0] -energy_err_low[0] ,energy_high[-1] +energy_err_high[-1] )
        axs[k+1].set_xlabel("Energy (keV)", fontsize = major_fontsize, labelpad = labelpadx)
        axs[k+1].set_ylabel("Ratio", fontsize = major_fontsize, labelpad = labelpady_2)
        
        tick_labels = [1,2,5]
        axs[k+1].set_xticks(tick_labels)
        axs[k+1].set_xticklabels(['{:g}'.format(tick) for tick in tick_labels])
        
        if (file_chooser == 3):
            axs[k+1].set_yticks ([0.95,1, 1.05, 1.1], [0.95,1, 1.05, 1.1]) #for rev0965
        axs[k+1].tick_params(axis="y", labelsize=minor_fontsize)
        axs[k+1].tick_params(axis="x", labelsize=minor_fontsize)
        axs[k+1].axhline(y = 1, color = "lime", linestyle = "-", alpha = 1)
        
        plt.subplots_adjust(wspace=0, hspace=0)

axs[0].tick_params(axis = "y",which = "minor", direction = "in", width = 2, length = 3)
axs[0].tick_params(axis = "y",which = "major", direction = "in", width = 2, length = 5)

axs[1].tick_params(axis = "both",which = "minor", direction = "in", width = 2, length = 3)
axs[1].tick_params(axis = "both",which = "major", direction = "in", width = 2, length = 5)

if (file_chooser == 0):
    
    filename_inset = "rev0956_4_8_keV_FeK_line_EXT"
    
    print("\n-----------------------------------------------")
    print("\nREADING FILE: %s\n" %filename_inset)
    
    line1 = "READ SERR 1 2"
    
    XRB = open("%s.qdp" %filename_inset, "r")
    lines = XRB.readlines()
    XRB.close()
    
    energy_inset = []
    flux_inset = []
    counts_inset = []
    energy_err_inset = []
    counts_err_inset = []
    Gauss_inset = []
    
    count = 0
    for line in lines:
        if (line == '%s\n' %line1 or 
            line == '@%s.pco\n' %filename or
            line == '!\n'):
            print("     IGNORED:", line)
        else:
            try:
                flux_temp = float(line.split()[4])
                flux_inset.append(flux_temp)
                
                energy_temp = float(line.split()[0])
                energy_inset.append(energy_temp)
                
                counts_temp = float(line.split()[2])
                counts_inset.append(counts_temp)
                
                energy_err_temp = float(line.split()[1])
                energy_err_inset.append(energy_err_temp)
                
                counts_err_temp = float(line.split()[3])
                counts_err_inset.append(counts_err_temp)
                           
                Gauss_temp = float(line.split()[6])
                Gauss_inset.append(Gauss_temp)
                
            except Exception:
                count = count + 1
    
    print ("no. Exceptions: ", count)
    energy_inset = np.array(energy_inset)   
    flux_inset = np.array(flux_inset) 
    counts_inset = np.array(counts_inset) 
    energy_err_inset = np.array(energy_err_inset)
    counts_err_inset = np.array(counts_err_inset)
    ratio_inset = (counts_inset/flux_inset)
    counts_err_inset = ratio_inset*(counts_err_inset/counts_inset)
    Gauss_inset = np.array(Gauss_inset)
    
    # create inset plot
    sub_ax = inset_axes(
        parent_axes=axs[0],
        width="45%",
        height="50%",
        loc='lower left',
        bbox_to_anchor=(0.05, 0.11, 1, 1),  # (x0, y0, width, height)
        bbox_transform=axs[0].transAxes,
        borderpad=0 
    )    
    
    sub_ax.hlines(y=flux_inset[0], xmin = energy_inset[0] - energy_err_inset[0], xmax=energy_inset[0] + energy_err_inset[0], color = colour[k], linewidth = 1.3)
    
    for n in range(1,len(energy_inset)):
        plt.hlines(y=flux_inset[n], xmin = energy_inset[n] - energy_err_inset[n], xmax=energy_inset[n] + energy_err_inset[n], color = colour[k], linewidth = 1.3)
        plt.vlines(x=energy_inset[n] - energy_err_inset[n], ymin=flux_inset[n], ymax=flux_inset[n-1], color = colour[k], linewidth = 1.3)
     
    sub_ax.errorbar(energy_inset, counts_inset, yerr = counts_err_inset, xerr = energy_err_inset, color = "k", linestyle = "", marker = "", linewidth = 1.3)

    plt.xlim(4,8)
    plt.ylim(1.1, 6)
    
    sub_ax.xaxis.set_major_locator(MultipleLocator(1))
    sub_ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    sub_ax.yaxis.set_minor_locator(MultipleLocator(0.2))
  
    tick_labels = [4,5,6,7,8]
    sub_ax.set_xticks(tick_labels)
    sub_ax.set_xticklabels(['{:g}'.format(tick) for tick in tick_labels])
    tick_labels = [2,4,6]
    sub_ax.set_yticks(tick_labels)
    sub_ax.set_yticklabels(['{:g}'.format(tick) for tick in tick_labels])
    
    sub_ax.tick_params(axis = "both", which='major', length=4, direction="in", width = 2, labelsize=minor_fontsize*0.8)
    sub_ax.tick_params(axis = "both", which='minor', length=2, direction="in", width = 2, labelsize=minor_fontsize*0.8)
    

if (file_chooser == 0):
    plt.savefig("GRO_%s_MCMC_%s_INSET.pdf" %(filenames[0].split("_")[0],date_now), bbox_inches = "tight")
else:
    plt.savefig("GRO_%s_MCMC_%s.pdf" %(filenames[0].split("_")[0],date_now), bbox_inches = "tight")
 
plt.show()
