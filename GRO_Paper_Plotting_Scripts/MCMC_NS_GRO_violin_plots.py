import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import seaborn as sns
from os.path import exists
import matplotlib.patches as mpatches

start = datetime.datetime.now()
print(start)
print('\n')

###############
#SEABORN VIOLIN PLOTS
###############


#Extracting Data
REVS = [
"rev0964a",
"rev0964b",
"rev0965",
"rev0966"
]

params = [  "H__2", 
            "O_I__10", "O_II__11", 
            "Ne_I__13","Ne_II__14",
            "Fe__31",
            "Gamma__33","FracSctr__34",
            "a__37", 
            "i__38","Mbh__39","Mdd__40","Dbh__41", "hd__42",
            "Index1__51", "Index2__52",
            "gamma__54", "logxi__55", "Afe__57",
            "norm__60", 
            ]

big_array = np.array([None]*len(params))
big_array_temp = np.array([None]*len(params))
no_walkers = 200
scale_down_factor = 600
skip_factor = no_walkers*scale_down_factor

exclude_range = 1

print("\nMCMC")
for looper in range(len(REVS)):
    testtime1 = datetime.datetime.now()
    
    for i in range(3):
        
        hdul = fits.open('MCMC/%s_afree_2M_SAMPLE%s.fits' %(REVS[looper],i+1))
        cols = hdul[1].columns
        data = hdul[1].data
        
        for n in range(len(params)):
            temp_param_array = data[params[n]]     
            
            if (i == 0):
                big_array_temp[n] = temp_param_array
            else:
                big_array_temp[n] = np.concatenate((big_array_temp[n], temp_param_array), axis=0)

    multiplier = int(len(big_array_temp[0])/scale_down_factor)
    if (looper == 0):
        rev_array = np.array([REVS[looper]]*multiplier)
        NS_MCMC_array = np.array(["MCMC"]*multiplier)
    else:
        rev_array = np.concatenate((rev_array, np.array([REVS[looper]]*multiplier)), axis = 0)
        NS_MCMC_array = np.concatenate((NS_MCMC_array, np.array(["MCMC"]*multiplier)), axis = 0)

    for n in range(len(params)):
        
        temp_param_array = big_array_temp[n]
        temp_array = np.array([])
        for k in range(no_walkers):
            temp_array = np.concatenate((temp_array, temp_param_array[k::skip_factor]))
        if (n == 11 or n == 19):
            temp_array = np.log10(temp_array)         
        if (looper == 0):
            big_array[n] = temp_array

        else:
            big_array[n] = np.concatenate((big_array[n], temp_array), axis=0)        

    print(looper, "time:", datetime.datetime.now() - testtime1)

 
    
print("\nNS")
for looper in range(len(REVS)):
    testtime1 = datetime.datetime.now()    
    
    Nested_Sampling = open("NS/equal_weighted_post_%s.txt" %REVS[looper].split("_")[0], "r")
    
    lines = Nested_Sampling.readlines()
    Nested_Sampling.close()       
    
    multiplier = len(lines) - 1
    rev_array = np.concatenate((rev_array, np.array([REVS[looper]]*multiplier)), axis = 0)
    NS_MCMC_array = np.concatenate((NS_MCMC_array, np.array(["NS"]*multiplier)), axis = 0)
    
    for n in range(len(params)):
        NS_temp = []
        for line in lines:
            
            if (line != lines[0]):
                NS_temp.append(float(line.split(" ")[n]))

        big_array[n] =  np.concatenate((big_array[n], np.array(NS_temp)), axis = 0)
    print(looper, multiplier, "time:", datetime.datetime.now() - testtime1)


print("\nDATAFRAME")
d = {     
     "spin": big_array[8], 
     "incl": big_array[9], 
     "Mass BH": big_array[10], 
     "Mdd": big_array[11], 

     "REV": rev_array,
     "NS or MCMC": NS_MCMC_array,
     }

df = pd.DataFrame(data=d, index=np.linspace(0, len(big_array[0]), len(big_array[0])))

print("\nPLOTTING")

yellow_colour = "#C1B732"
green_colour = "#67D8B9"
pink_colour = "#DC267F"
blue_colour = "#648FFF"

y_labels = [
            "spin",
            r"Inclination ($^{\circ}$)", r"Mass (M$_{\odot}$)", 
            r"$\log{(\dot{M}})$",#" ($\log{(10^{18}}$ g s$^{-1}$))",
            ]

print("\nTOTAL TIME: ", datetime.datetime.now() - start)


no_plots = len(y_labels)
fig = plt.figure(figsize=(12, 5*2))
grid = fig.add_gridspec(2,2) #rows then columns

label_array = []

df_params = [
             "spin",
             "incl", "Mass BH", 
             "Mdd"
             ]

for n in range(len(df_params)):
    start_plot_time = datetime.datetime.now()
    k = int(n/2)
           
    ax = fig.add_subplot(grid[k, n%2])

    temp = sns.violinplot(data=df, x="REV", y=df_params[n], 
                   hue='NS or MCMC', 
                   split=True,
                   palette=[blue_colour,pink_colour], 
                   cut = 0, bw=.2, inner = None
                   )
    
    label_array.append(temp)
    ax.get_legend().remove()
    
    mcmc_patch = mpatches.Patch(color=blue_colour, label="MCMC")
    ns_patch = mpatches.Patch(color=pink_colour, label="Nested Sampling")
    
    # Add legend below the bottom plots, centered
    fig.legend(handles=[mcmc_patch, ns_patch], loc="lower center", fontsize=30, 
               ncol=2, bbox_to_anchor=(0.5, -0.1), frameon=False)
    
    plt.xlabel('GRO', size=34, labelpad = 10)
    plt.ylabel(y_labels[n], size=34, labelpad = 15)
    plt.yticks(fontsize=26)
    plt.xticks(fontsize=26)
    
    plt.xticks([
                0,1,2,3,
                ], 
                [
                2,3,4,5
                ])
               
    print("PLOTTING TIME: ", df_params[n], ":", datetime.datetime.now() - start_plot_time)
    
    plt.tick_params(axis="both", direction = "in", width = 1.5, length = 4)
    
    sns.despine()
    sns.set_style("ticks")
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
     
    if (n == 0):
        loc = [-1.0, -0.5, 0.0, 0.5]
        plt.ylabel("Spin", size=32, labelpad = 10)
        plt.xlabel("")
        ax.set_xticklabels([])    
    if (n == 1):
        loc = [84.2, 84.4, 84.6, 84.8, 85.0]
        plt.ylabel(y_labels[n], labelpad=7)
        plt.xlabel("")
        ax.set_xticklabels([])
    if (n == 2):
        loc = [4,5,6,7,8,9]
        plt.ylabel(y_labels[n], labelpad=45)
    if (n == 3):
        loc = [0.0, 0.4, 0.8, 1.2, 1.6]
        ax.set_ylim(0, 1.7)
    ax.set_yticks(loc)
    
    
labels = ["MCMC", "Nested\nSampling"]   
 
plt.subplots_adjust(wspace = 0.5, hspace = 0.03)
plt.text(4, 0.8, " ", color = "white")

date_now = str(datetime.datetime.now())[0:10]
plt.savefig("GRO_Violin_Plots_%s.pdf" %date_now, bbox_inches="tight")
plt.show()

print("\nTOTAL TIME: ", datetime.datetime.now() - start)


