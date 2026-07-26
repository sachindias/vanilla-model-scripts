import datetime
from NS_MCMC_comparison_violin_plots import get_params, extract_MCMC_data, extract_NS_data, combine_NS_MCMC_arrays, violin_plotter
    
#TIMER
start_timer = datetime.datetime.now()
print(start_timer)
print('\n')

#ALL MCMC FILES
all_MCMC_base_filenames = [
    "rev0956_afree_2M_SAMPLE",
    "rev0964a_afree_2M_SAMPLE",
    "rev0964b_afree_2M_SAMPLE",
    "rev0965_afree_2M_SAMPLE",
    "rev0966_afree_2M_SAMPLE",
    "rev0970_afree_2M_SAMPLE",
    ]
MCMC_file_location = "MCMC"

#ALL NS FILES
all_NS_filenames = [
    "equal_weighted_post_rev0956",
    "equal_weighted_post_rev0964a",
    "equal_weighted_post_rev0964b",
    "equal_weighted_post_rev0965",
    "equal_weighted_post_rev0966",
    "equal_weighted_post_rev0970"
    ]
NS_file_location = "NS"

#REV NAMES
GRO = [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        ]

all_combined_values = []
all_params = []
 
for n in range(len(all_MCMC_base_filenames)):
    print("ENSURING SAME PARAMS ON: GRO%s\n" %GRO[n])
    MCMC_base_filenames = all_MCMC_base_filenames[n]
    params = get_params(MCMC_file_location, MCMC_base_filenames, start_timer)
    all_params.append(params)
    
#ENUSRES THAT ALL THE SAME PARAMETERS ARE USED
#THIS IS JUST A SIMPLE PLOTTER AND NEEDS ALL PARAMETERS TO BE THE SAME
#SOMETHING MORE BESPOKE COULD BE WRITTEN TO DEAL WITH THOSE CASES
if (len(set(len(x) for x in all_params)) != 1):
    check_spin_present = 1 if 'a__37' in all_params[0] else 0
    for n in range(len(all_params)):
        if (('a__37' in all_params[n]) != check_spin_present):
            check_spin_present = 2
            
    if (check_spin_present == 2):
        for params in all_params:
            if 'a__37' in params:
                params.remove('a__37')
                
if (len(set(len(x) for x in all_params)) != 1):
    min_length = min(len(x) for x in all_params)
    for n in range(len(all_params)):
        all_params[n] = all_params[n][0:min_length]      

#DEAL WITH EACH REVOLUTION IN TURN & SAVE IN LARGE OVERALL ARRAY
for n in range(len(all_MCMC_base_filenames)):
    print("COMBINING VALUES FOR: GRO%s\n" %GRO[n])
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
    combined_values = combine_NS_MCMC_arrays(MCMC_values, NS_values, params, GRO[n])
    all_combined_values.append(combined_values)

#MAKE THE VIOLIN PLOTS
violin_plotter(all_combined_values, all_params[0], "REV", "GRO", save_filename = "GRO_MCMC_NS_Violin_Comparison", plot_in_console = True)

print("\nTOTAL TIME: ", datetime.datetime.now() - start_timer)

