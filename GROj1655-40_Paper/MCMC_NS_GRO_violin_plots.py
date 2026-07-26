import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches

#THE BELOW FUNCTIONS ARE IN A SEPARATE FILE TO MAKE THIS SCRIPT EASIER TO READ
#THIS FILE SHOULD BE IN THE SAME DIRECTORY AS THIS SCRIPT BEFORE RUNNING
from NS_MCMC_comparison_violin_plots import get_params, extract_MCMC_data, extract_NS_data, combine_NS_MCMC_arrays, combined_NS_MCMC_array_concatenated

#TIMER
start_timer = datetime.datetime.now()
print(start_timer)
print('\n')

pink_colour = "#DC267F"
blue_colour = "#648FFF"

REVS = [
    "rev0964a",
    "rev0964b",
    "rev0965",
    "rev0966"
    ]

all_MCMC_base_filenames = []
all_NS_filenames = []
all_params = []
all_combined_values = []

MCMC_file_location = "MCMC"
NS_file_location = "NS"

#GET MCMC & NS FILENAMES & PARAMS
for n in range(len(REVS)):
    MCMC_base_filename = ("%s_afree_2M_SAMPLE" %REVS[n])
    NS_filename = ("equal_weighted_post_%s" %REVS[n])
    
    all_MCMC_base_filenames.append(MCMC_base_filename)
    all_NS_filenames.append(NS_filename)
    
    params = get_params(MCMC_file_location, MCMC_base_filename, start_timer)
    all_params.append(params)

#ENUSRES THAT ALL THE SAME PARAMETERS ARE USED
if (len(set(len(x) for x in all_params)) != 1):
    min_length = min(len(x) for x in all_params)
    for n in range(len(all_params)):
        all_params[n] = all_params[n][0:min_length]  

#DEAL WITH EACH REVOLUTION IN TURN & SAVE IN LARGE OVERALL ARRAY
for n in range(len(all_MCMC_base_filenames)):
    print("COMBINING VALUES FOR: GRO%s\n" %(n+1))
    #SELECT FILE
    MCMC_base_filenames = all_MCMC_base_filenames[n]
    NS_filename = all_NS_filenames[n]
    params = all_params[n]
    
    #MCMC
    MCMC_values = extract_MCMC_data(MCMC_file_location, MCMC_base_filenames, params, start_timer,
                                no_walkers = 200, scale_down_factor = 600)
        
    #NESTED SAMPLING
    NS_values = extract_NS_data(NS_file_location, NS_filename, params)

    #COMBINING
    combined_values = combine_NS_MCMC_arrays(MCMC_values, NS_values, params, REVS[n])
    all_combined_values.append(combined_values)
    
    #CONCATENATE INTO ONE ARRAY
    concat_values = combined_NS_MCMC_array_concatenated(all_combined_values)

#SORT INTO PANDAS DATAFRAME
print("\nDATAFRAME")
d = {     
     "spin": concat_values[8], 
     "incl": concat_values[9], 
     "Mass BH": concat_values[10], 
     "Mdd": np.log10(concat_values[11]), 

     "REV": concat_values[-2],
     "NS or MCMC": concat_values[-1],
     }

df = pd.DataFrame(data=d, index=range(len(concat_values[0])))


print("\nPLOTTING, TOTAL TIME: ", datetime.datetime.now() - start_timer)

y_labels = [
            "spin",
            r"Inclination ($^{\circ}$)", 
            r"Mass (M$_{\odot}$)", 
            r"$\log{(\dot{M}})$",#" ($\log{(10^{18}}$ g s$^{-1}$))",
            ]

fig = plt.figure(figsize=(12, 5*2))
grid = fig.add_gridspec(2,2) #ROWS THEN COLUMNS

#PARAMETERS TO PLOT
df_params = [
             "spin",
             "incl", "Mass BH", 
             "Mdd"
             ]

#PLOTS PARAMETERS ONE BY ONE
for n in range(len(df_params)):
    k = int(n/2)    
    ax = fig.add_subplot(grid[k, n%2])

    per_rev_violin_plot = sns.violinplot(data=df, x="REV", y=df_params[n], 
                   hue='NS or MCMC', 
                   split=True,
                   palette=[blue_colour,pink_colour], 
                   cut = 0, bw=.2, inner = None
                   )
    ax.get_legend().remove()
    
    #ADD LEGEND BELOW THE FINAL PLOT
    mcmc_patch = mpatches.Patch(color=blue_colour, label="MCMC")
    ns_patch = mpatches.Patch(color=pink_colour, label="Nested Sampling")
    fig.legend(handles=[mcmc_patch, ns_patch], loc="lower center", fontsize=30, 
               ncol=2, bbox_to_anchor=(0.5, -0.1), frameon=False)
    
    #PLOT SETTINGS
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

    print("PLOTTING TIME: ", df_params[n], ":", datetime.datetime.now() - start_timer)
 
#FINAL ADJUSTMENTS
plt.subplots_adjust(wspace = 0.5, hspace = 0.03)
plt.text(4, 0.8, " ", color = "white")

#SAVING & PLOTTING
date_now = str(datetime.datetime.now())[0:10]
plt.savefig("GRO_Violin_Plots_%s.pdf" %date_now, bbox_inches="tight")
plt.show()

print("\nTOTAL TIME: ", datetime.datetime.now() - start_timer)
