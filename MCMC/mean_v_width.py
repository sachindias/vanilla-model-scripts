import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import datetime
import os
import math

#-----------------------------------------
#FUNCTIONS

#CHECK WHICH PARAMETERS ARE NEEDED FOR PLOTTING
def check_params_present(hdus_suffix, base_filename):
    
    #ALL FITS FILES SHOULD HAVE THE SAME FIELDS
    hdul = fits.open('%s%s.fits' %(base_filename, hdus_suffix))
    cols = hdul[1].columns
    data = hdul[1].data
    
    checked_params = []
    
    #TEST FOR N_O_III AND N_NE_III
    try:
        OIII_array = data['O_III__12']
        NeIII_array = data['Ne_III__15']
        print("N_O_III AND N_NE_III INCL.")
        checked_params.extend([['O_III__12'], ['Ne_III__15']])
    except Exception:
        checked_params.extend(['NO_DATA','NO_DATA'])
        
    #TEST IF SPIN FREE
    try:
        a_array = data["a__37"]
        print("SPIN FREE")
        checked_params.extend([["a__37"]])
    except Exception:
        checked_params.extend(['NO_DATA'])
        
    #CHECK FOR GABS
    try:
        LineE_array = data['LineE__61']
        Strength_array1 = data['Strength__63']
        Strength_array2 = data['Strength__66']
        Strength_array3 = data['Strength__69']
        checked_params.extend([['LineE__61', 'Strength__63', 'Strength__66', 'Strength__69']])
        print("GABS INCL.")
    except Exception:
        checked_params.extend(['NO_DATA'])
        
    #CHECK FOR SMEDGE
    try:
        Smedge_E_array = data['edgeE__73']
        Smedge_Tau_array = data['MaxTau__74']
        Smedge_w_array = data['width__76']
        print("SMEDGE INCL.")
        checked_params.extend([['edgeE__73', 'MaxTau__74', 'width__76']])
    except Exception:
        checked_params.extend(['NO_DATA'])
        
    #CHECK FOR NUSTAR DATA
    try:
        factor_1_array = data['factor__61']
        factor_2_array = data['factor__121']
        checked_params.extend([['factor__61', 'factor__121']])
        print("NuSTAR DATA INCL.")
    except Exception:
        print("ONLY XMM DATA")
        checked_params.extend(['NO_DATA'])
        
    return checked_params

#COLLATES THE LIST OF ALL PARAMETERS FOR CHECK FITS FILES
def get_params(hdu_suffix, base_filename, start_timer):
    
    #PARAMETERS THAT MIGHT BE MISSING
    checked_params = check_params_present(hdu_suffix, base_filename)
    
    #PARAMETERS THAT ARE ALWAYS PRESENT
    standard_params = [
        ["H__2", "O_I__10", "O_II__11"],
        ["Ne_I__13","Ne_II__14"],
        [ "Fe__31","Gamma__33","FracSctr__34"],
        ["i__38","Mbh__39","Mdd__40","Dbh__41", "hd__42", 
         "Index1__51", "Index2__52", "gamma__54",
         "logxi__55", "Afe__57", "norm__60"],
         ]
    
    #SORT PARAMETERS INTO THE CORRECT ORDER
    params = ["FIT_STATISTIC"]
    for n in range(len(checked_params)):
        try:
            params.extend(standard_params[n])
        except Exception:
            print("RAN OUT OF STANDARD PARAMS, CONTINUING WITH CHECKED ONLY")
        if (checked_params[n] != 'NO_DATA'):
            params.extend(checked_params[n])

    print("\nPARAM LIST MADE (LENGTH: %s PARAMS), TIME:" %len(params), datetime.datetime.now() - start_timer)  
    print("--------------------------------------------\n") 
            
    return params

#CHECKS THE LENGTH OF BATCHES AND FINDS WHICH BATCHES HAVE WHICH ITERATIONS
def length_checker(base_filename, hdu_list, no_walkers):
    
    iteration_counter = 0
    for n in range(len(hdu_list)):
        #LOADS DATA
        hdul = fits.open('%s%s.fits' %(base_filename, hdu_list[n]))
        cols = hdul[1].columns
        data = hdul[1].data

        #START ITERATION OF BATCH (OR GROUP OF BATCHES) AND
        #TOTAL COUNT OF ITERATIONS AT THE END OF THE BATCH
        iteration_starting_point = iteration_counter + 1
        iteration_counter = iteration_counter + len(data)/no_walkers
        
        #MACTHES ITERATIONS WITH BATCHES (OR GROUPS OF BATCHES)
        if (n == 0):
            all_iterations = np.linspace(iteration_starting_point, int(iteration_counter), int(len(data)/no_walkers))
            batch_locatons = np.repeat(hdu_list[n], len(data)/no_walkers)
            start_iterations = np.repeat(iteration_starting_point, len(data)/no_walkers)
        else:
            all_iterations = np.concatenate((all_iterations, np.linspace(iteration_starting_point, int(iteration_counter), int(len(data)/no_walkers))), axis = 0)
            batch_locatons = np.concatenate((batch_locatons, np.repeat(hdu_list[n], int(len(data)/no_walkers))), axis = 0)
            start_iterations = np.concatenate((start_iterations, np.repeat(iteration_starting_point, int(len(data)/no_walkers))), axis = 0)

        
        print("%s%s.fits : %s : %s" %(base_filename,hdu_list[n], len(data), iteration_counter))
        hdul.close()

    print("\nLENGTH CHECKED, TIME:", datetime.datetime.now() - start_timer) 
    print("--------------------------------------------") 
    
    #all_iterations: EVERY ITERATION ACROSS ALL BATCHES
    #batch_locatons: WHICH BATCH EACH ITERATION IS IN
    #start_iterations: THE BATCH STARTING ITERATION FOR EACH ITERATION
    return all_iterations, batch_locatons, start_iterations

#CACULATES RELATIVE CHANGE IN MEAN RELATIVE TO THE WIDTH OF THE DISTRIBUTION
def change_mean_v_width(Sample1,Sample2):
    
    #CALCULATE MEAN & ERROR
    mean1 = np.mean(Sample1)
    mean2 = np.mean(Sample2)
    
    q_16_1 = np.quantile(Sample1, 0.16)
    q_84_1 =np.quantile(Sample1, 0.84)
    width_1 = q_84_1 - q_16_1 
    
    q_16_2 = np.quantile(Sample2, 0.16)
    q_84_2 = np.quantile(Sample2, 0.84)
    width_2 = q_84_2 - q_16_2
    
    rel_change_mean = np.sqrt(2)*(mean1 - mean2)/np.sqrt(width_1**2 + width_2**2)
    
    return rel_change_mean, width_1, width_2

#PLOTS THE DIFFERENT DISTRIBUTIONS FOR THE VARIOUS BATCH
def plotter(MAXI_or_GRO, params, base_filename, selected_batch_labels,
            start_timer, plot_in_console, save_filename = False, comparison_batch_locations = False,
            comparison_locs = False, start_position_in_batch = False, selected_batches = False):
    
    print("\nPLOTTING\n")
    
    #SET UP FIGURE SIZES AND THE NUMBER OF PLOTS NEEDED
    plt.figure(figsize=(15, 10*len(params)))
    number_plots = len(params)

    colours = [
        "#E74C3C",  # RED
        "#2ECC71",  # GREEN
        "#3498DB",  # SKY BLUE
        "#F39C12",  # GOLDEN ORANGE
        "#F1C40F",  # YELLOW
        "#9B59B6",  # PURPLE
        "#E67E22",  # ORANGE
        "#1ABC9C",  # TEAL
        "#FF5733",  # BRIGHT RED-ORANGE
        "#33FF57",  # BRIGHT GREEN
        "#3357FF",  # BRIGHT BLUE
        "#8E44AD",  # DEEP PURPLE
        "#16A085",  # SEA GREEN
        "#27AE60",  # EMERALD
        "#2980B9",  # DEEP BLUE
        "#D35400",  # BURNT ORANGE
        "#C0392B",  # DARK RED
        "#7F8C8D",  # GRAY
        "#34495E"   # NAVY GRAY
        ]

    #PLOT DISTRIBUTIONS & MEASURE OF WIDTH
    for n in range(len(params)):
        plt.subplot(number_plots*2,2,n+1)
        plot_timer = datetime.datetime.now() #ADDITIONAL TIMER FOR PLOT TIMINGS
        
        for k in range(len(selected_batch_labels)):
            if (MAXI_or_GRO == "GRO"):
                #LOAD GRO DATA
                hdul = fits.open('%s%s.fits' %(base_filename, selected_batches[k]))
                cols = hdul[1].columns
                data = hdul[1].data
                
            elif (MAXI_or_GRO == "MAXI"):
                #LOAD MAXI DATA
                hdul = fits.open('%s%s.fits' %(base_filename, comparison_batch_locations[k]))
                data = hdul[1].data[int((comparison_locs[k]-start_position_in_batch[k]+1)*200):int((comparison_locs[k]-start_position_in_batch[k]+1+1000)*200)]
            
            #PLOT DISTRIBUTION, MEAN & WIDTH
            plt.hist(data[params[n]], bins = 20, histtype = "step", linewidth = 2, color = colours[k], label = selected_batch_labels[k], alpha = 0.5)
            plt.axvline(x = np.mean(data[params[n]]), linestyle = "--", color = colours[k])
            plt.axvline(x = np.quantile(data[params[n]], 0.16), linestyle = "-", color = colours[k])
            plt.axvline(x = np.quantile(data[params[n]], 0.84), linestyle = "-", color = colours[k])
           
        #PLOTTING SETTINGS
        plt.legend() 
        plt.yticks([])
        plt.xticks(fontsize = 14)
        plt.xlabel(params[n], fontsize = 16)
        
        print('%s PARAM: %s, PLOT TIME: %s' %(n, params[n], datetime.datetime.now() - plot_timer))

    #SAVE PLOT TO A PNG FILE
    if (save_filename != False):
        save_timer = datetime.datetime.now() #ADDITIONAL TIMER FOR PLOT SAVING
        plt.savefig("%s.png" %save_filename, bbox_inches = "tight")
        print("\nSAVE TIME: %s" %(datetime.datetime.now() - save_timer))

    #PLOT IN CONSOLE
    if (plot_in_console):
        plt.show()

    print("--------------------------------------------")
    print("DATA PLOTTED, TIME: %s\n" %(datetime.datetime.now() - start_timer))  

#CALCUALTE CHANGE IN MEAN RELATIVE TO WIDTH
def print_rel_change(MAXI_or_GRO, params, base_filename, selected_batch_labels,
                     comparison_batch_locations = False,
                     comparison_locs = False,
                     start_position_in_batch = False,
                     selected_batches = False):
    print("\nCHANGE IN MEAN RELATIVE TO WIDTH \n")

    final_change_in_mean = [] #KEEPS TRACK OF THE FINAL RELATIVE CHANGE IN THE MEAN
    comparison_strings = [] #KEEPS TRACK OF STRING TO PRINT ALL RELATIVE CHANGES IN THE MEAN

    for n in range(len(params) - 1):
        param_timer = datetime.datetime.now() #ADDITIONAL TIMER FOR EACH PARAMETER

        param_string = params[n] + " "*(12-len(params[n])) + ";"
        
        for k in range(len(selected_batch_labels)-1):
            if (MAXI_or_GRO == "GRO"):
                #LOADS UP TWO BATCHES TO COMPARE
                hdul = fits.open('%s%s.fits' %(base_filename, selected_batches[k]))
                cols = hdul[1].columns
                data = hdul[1].data
                
                hdul2 = fits.open('%s%s.fits' %(base_filename, selected_batches[k+1]))
                cols2 = hdul2[1].columns
                data2 = hdul2[1].data
            
            elif (MAXI_or_GRO == "MAXI"):
                #LOADS UP TWO LOCATIONS TO COMPARE (AT COMPARISON LOCATIONS + 1000 ITERATIONS)
                hdul = fits.open('%s%s.fits' %(base_filename, comparison_batch_locations[k]))
                data = hdul[1].data[int((comparison_locs[k]-start_position_in_batch[k]+1)*200):int((comparison_locs[k]-start_position_in_batch[k]+1+1000)*200)]
                
                hdul2 = fits.open('%s%s.fits' %(base_filename, comparison_batch_locations[k+1]))
                data2 = hdul2[1].data[int((comparison_locs[k+1]-start_position_in_batch[k+1]+1)*200):int((comparison_locs[k+1]-start_position_in_batch[k+1]+1+1000)*200)]
                
            
            #COMPARES BATCHES
            rel_change_mean, width_1, width_2 =  change_mean_v_width(data[params[n]], data2[params[n]]) 
            
            #TRACKS IF SPACES NEEDED FOR STRING
            if (rel_change_mean < 0):
                space_values1 = 0
            else:
                space_values1 = 1
                
            #TURNS WIDTH OF DISTRIBUTION INTO POWER OF 10
            scaling_value = max([width_1, width_2])
            scaling_magnitude = math.floor(np.log10(scaling_value))
            new_width_1 = width_1/10**(math.floor(np.log10(scaling_value)))
            new_width_2 = width_2/10**(math.floor(np.log10(scaling_value)))

            #TRACKS IF SPACES NEEDED FOR STRING
            if (scaling_magnitude < 0):
                space_values2 = 0
            else:
                space_values2 = 1

            #CREATES STRING
            param_string = param_string + " " + selected_batch_labels[k] + " V " + selected_batch_labels[k+1] + " " + " "*space_values1 + ("%.4f" %rel_change_mean) + " " + ("(%.4f, %.4f)e%s" %(new_width_1, new_width_2, scaling_magnitude)) + " "*space_values2 + " ;"
        
        #APPEND TO ARRAYS
        comparison_strings.append(param_string)   
        final_change_in_mean.append(rel_change_mean)    
        print('%s PARAM: %s, INDIVIDUAL TIME: %s' %(n, params[n], datetime.datetime.now() - param_timer))

    print("\n\n")
    print("COMPARISONS :")
    for n in range(len(comparison_strings)):
        print(comparison_strings[n])

    print("\nMEAN OF FINAL COMPARISSONS: %.4f PERCENT" %(100*np.mean(abs(np.array(final_change_in_mean)))))

def main_function(MAXI_or_GRO, hdu_list, base_filename, save_filename, start_timer, plot_in_console = False):
    params = get_params(hdu_list[0], base_filename, start_timer)


    if (MAXI_or_GRO == "GRO"):
        #SELECTS THE BATCHES TO USE (& ADDS LABELS TO THEM)
        #AC STANDS FOR 'AFTER CLUSTERING'
        #REFERING TO THE BATCH AFTER THE LAST CLUSTERING FIX
        #IN THIS EXAMPLE THERE ISN'T ANY CLUSTING FIX
        selected_batches = [""]
        selected_batch_labels = ["AC"] 
    
        for n in range(len(hdu_list)):
            if (n > 0):
                if (n%20 == 0):
                    selected_batches.append(hdu_list[n])
                    selected_batch_labels.append("AC +%s" %n)
    
        ####################################
        #PLOTTING
        ####################################
        
        plotter(MAXI_or_GRO, params, base_filename, selected_batch_labels,
                start_timer, plot_in_console, save_filename, 
                comparison_batch_locations = False,
                comparison_locs = False,
                start_position_in_batch = False, 
                selected_batches = selected_batches)
    
        ###########################################
        #CALCUALTE CHANGE IN MEAN RELATIVE TO WIDTH 
        ###########################################
    
        print_rel_change(MAXI_or_GRO, params, base_filename, selected_batch_labels,
                             comparison_batch_locations = False,
                             comparison_locs = False,
                             start_position_in_batch = False,
                             selected_batches = selected_batches)
    
    elif (MAXI_or_GRO == "MAXI"):
        #GETS WHICH BATCHES HAVE WHICH ITERATIONS & STARTING ITERATIONS FOR THOSE BATCHES
        all_iterations, batch_locatons, start_iterations = length_checker(base_filename, hdu_list, no_walkers)

        comparison_locs = [9e3, 4.9e4, 9.9e4, 1.99e5]
        start_position_in_batch = []

        comparison_batch_locations = []
        for n in range(len(comparison_locs)):
            which_batches = np.where(all_iterations == comparison_locs[n])
            comparison_batch_locations.append(batch_locatons[which_batches][0])
            start_position_in_batch.append(start_iterations[which_batches][0])

        params = get_params(hdu_list[0], base_filename, start_timer)

        selected_batch_labels  = ["AC", "AC+5", "AC+10", "AC+20"]

        ####################################
        #PLOTTING
        ####################################

        plotter(MAXI_or_GRO, params, base_filename, selected_batch_labels,
                start_timer, plot_in_console, save_filename, 
                comparison_batch_locations = comparison_batch_locations,
                comparison_locs = comparison_locs,
                start_position_in_batch = start_position_in_batch, 
                selected_batches = False)

        ###########################################
        #CALCUALTE CHANGE IN MEAN RELATIVE TO WIDTH 
        ###########################################

        print_rel_change(MAXI_or_GRO, params, base_filename, selected_batch_labels,
                             comparison_batch_locations = comparison_batch_locations,
                             comparison_locs = comparison_locs,
                             start_position_in_batch = start_position_in_batch,
                             selected_batches = False)

#-----------------------------------------
#MAIN SCRIPT 1

print("--------------------------------------------") 

#TIMER
start_timer = datetime.datetime.now()
print(start_timer)
print('\n')

base_filename = "rev0956_27_9_23_afree_2M"
save_filename = 'GRO1_mean_v_error'

no_walkers = 200
iterations_per_batch = 10000
plot_in_console = True

#THIS WOULD HAVE TO BE MANUALLY SET PER CHAIN TO CAPTURE ALL THE FILENAMES
#THE BELOW EXAMPLE IS FOR GRO1
hdu_list = [   
        "",
        "_ext",
        ]

ext_condition = 0
file_counter = 1

while(ext_condition < 1):
    try:
        file_counter = file_counter + 1
        if (os.path.getsize("%s_ext%s.fits" %(base_filename, file_counter)) > 0):
            hdu_list.append("_ext%s" %file_counter)
            print("%s_ext%s.fits" %(base_filename, file_counter))
    except Exception:
        ext_condition = 2 
   
print("DATA SORTED, TIME:", datetime.datetime.now() - start_timer, "\n")         
print("--------------------------------------------") 
    
main_function("GRO", hdu_list, base_filename, save_filename, start_timer, plot_in_console = True)

#WRAPPING UP
print("TOTAL TIME: %s" %(datetime.datetime.now() - start_timer))

#-----------------------------------------
#MAIN SCRIPT 2

print("--------------------------------------------") 

#TIMER
start_timer = datetime.datetime.now()
print(start_timer)
print('\n')

base_filename = "rev3533_1_9_22_amax_2M"
save_filename = 'Obs2_amax_mean_v_error'

no_walkers = 200
length_each_run = 10000

#THIS WOULD HAVE TO BE MANUALLY SET PER CHAIN TO CAPTURE ALL THE FILENAMES
#THE BELOW EXAMPLE IS FOR MAXI J1820+070, Obs 2, a = 0.998
hdu_list = []

#THE BATCHING IS QUITE COMPLICATED AS SOME BATCHES CONTAIN THE LAST BATCH TOO
#THIS CHECKS THROUGH AND DEALS WITH THAT
try:
    if (os.path.getsize("%s_P2_ext.fits" %base_filename)/1000/1000 < 500):
        hdu_list.append("_P2")
    else:
        print("EXTENDED FILE: P2")
except Exception:
    print("EXTENDED FILE P2_ext DOES NOT EXIST")

n = 2
while (n < 20): 
    try:
        if (os.path.getsize("%s_P2_ext%s.fits" %(base_filename, n+1))/1000/1000 < 500):
            hdu_list.append("_P2_ext%s" %n)
            n = n + 1
        else:
            print("EXTENDED FILE: P2_ext%s" %n)
            n = n + 1
    except Exception:
        print("FILE DOES NOT EXIST: P2_ext%s" %(n+1))
        hdu_list.append("_P2_ext%s" %n)
        break
    
print("\nDATA SORTED", datetime.datetime.now() - start_timer) 
print("--------------------------------------------") 
    
main_function("MAXI", hdu_list, base_filename, save_filename, start_timer, plot_in_console = True)

#WRAPPING UP
print("TOTAL TIME: %s" %(datetime.datetime.now() - start_timer))
