import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import datetime
import os

#TIMER
start_timer = datetime.datetime.now()
print(start_timer)
print('\n')

#-----------------------------------------
#FUNCTIONS

#A VISUAL CHECK TO ENSURE ALL BATCHES ARE OF THE CORRECT LENGTH
#AND RETURNS THE NUMBER OF BATCHES NOT BEING PLOTTED
def length_checking(no_batches_to_plot, hdu_suffixes, base_filename, start_timer):
    
    skipped_batches = 0
    
    if (no_batches_to_plot != 'All'):
        #PLOTS THE LAST no_batches_to_plot BATCHES ONLY
        for n in range(len(hdu_suffixes) - no_batches_to_plot):
            hdu_suffixes.pop(0)
            skipped_batches = skipped_batches + 1
        
    for h in range(len(hdu_suffixes)):
        hdul = fits.open('%s%s.fits' %(base_filename, hdu_suffixes[h]))
        cols = hdul[1].columns
        data = hdul[1].data
        
        print("%s%s.fits; LENGTH: %s" %(base_filename, hdu_suffixes[h], len(data)))
        hdul.close()
    
    print("\nLENGTH CHECKED, TIME:", datetime.datetime.now() - start_timer)  
    print("--------------------------------------------\n") 
    
    return skipped_batches, hdu_suffixes
 
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

#EXTRACTS THE DATA FROM EACH OF THE BATCHES AND PUTS INTO AN ARRAY PER PARAMETER
def extract_trace_plots_data(hdu_suffixes, base_filename, params, start_timer):

    i = 1
    params_array_walker = [None]*len(params) #EMPTY ARRAY FOR ALL WALKER LOCATIONS
    
    #EXTRACTS DATA FOR EACH FIELD AND EACH BATCH AND ADDS TO THE ARRAY
    while (i < len(params) + 1):
        for n in range(len(hdu_suffixes)):   
            hdul = fits.open('%s%s.fits' %(base_filename, hdu_suffixes[n]))
            cols = hdul[1].columns
            data = hdul[1].data
            
            param_array = data[params[i-1]]
            
            if (n == 0):
                params_array_walker[i-1] = param_array
            else:
                params_array_walker[i-1] = np.hstack((params_array_walker[i-1], param_array))
                
        i = i + 1
        print('PARAM: %s, TIME: %s' %(params[i - 2], datetime.datetime.now() - start_timer))

    print("\nDATA EXTRACTED\n")  
    print("--------------------------------------------")
    
    return params_array_walker

#PLOT NORMAL TRACE PLOTS
def plot_trace_plots(params_array_walker, params, start_timer, no_walkers = 200, 
                     thin_factor = False, iterations_per_batch = 10000,
                     save_filename = False, plot_in_console = False,
                     no_individual_walkers_to_plot = False):
    
    #GET THE NUMBER OF ITERATIONS PER WALKER
    no_iterations = int(len(params_array_walker[0])/no_walkers)
    iterations = np.linspace(1,no_iterations,no_iterations)
    
    #SET UP FIGURE SIZES AND THE NUMBER OF PLOTS NEEDED
    plt.figure(figsize=(30, 10*len(params)))
    number_plots = len(params)
    
    #PLOT THE CHAIN FOR EACH WALKER
    for j in range(len(params_array_walker)):
        plt.subplot(number_plots*2,2,j+1)
        plot_timer = datetime.datetime.now() #ADDITIONAL TIMER FOR PLOT TIMINGS
        
        for i in range(no_walkers):
            walker_chain = params_array_walker[j][i::no_walkers]
            
            #THINNING IF THERE IS A LOT OF DATA PRESENT
            if (thin_factor == False):
                plt.plot(iterations, walker_chain, 'gray', alpha = 0.2)
            else:
                thinned_iterations = iterations[0::thin_factor]
                thinned_walker_chain = walker_chain[0::thin_factor]
                
                plt.plot(thinned_iterations, thinned_walker_chain, 'gray', alpha = 0.2)
        
        #PLOT INDIVIDUAL WALKERS FOR CLARITY (1 - 20)
        if (no_individual_walkers_to_plot != False):
            plot_individual_walkers(j, iterations, params_array_walker, no_individual_walkers_to_plot = 4, no_walkers = 200)
        
        #PLOT LINES TO INDICATE BATCHES
        for n in range(int(len(iterations)/iterations_per_batch)):
            plt.axvline(x = iterations_per_batch*(n), color = "blue", linestyle = '--', alpha = 0.3)
        
        #PLOT SETTINGS
        plt.xlim(1, no_iterations)   
        plt.ylabel(params[j], fontsize = 14)
        plt.xlabel('iterations', fontsize = 14)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
    
        print('PARAM: %s, PLOT TIME: %s' %(params[j], datetime.datetime.now() - plot_timer))
    
    #INCREASE SIZE BETWEEN PLOTS
    plt.tight_layout(pad=3.0)
    
    #SAVE PLOT TO A PNG FILE
    if (save_filename != False):
        save_timer = datetime.datetime.now() #ADDITIONAL TIMER FOR PLOT SAVING
        plt.savefig("%s_trace_plots.png" %save_filename, bbox_inches = "tight")
        print("\nSAVE TIME: %s" %(datetime.datetime.now() - save_timer))
    
    #PLOT IN CONSOLE
    if (plot_in_console):
        plt.show()
        
    print("\nTRACE PLOTS CREATED, TIME:", datetime.datetime.now() - start_timer)  
    print("--------------------------------------------\n") 

#PLOT INDIVIDUAL WALKER CHAINS TO SEE HOW THEY MOVE AS OPPOSED TO THE BULK
def plot_individual_walkers(param_counter, iterations, params_array_walker, no_individual_walkers_to_plot, no_walkers = 200):
    
    colours = [
    "#000000",  # BLACK
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
    
    for n in range(no_individual_walkers_to_plot):
        plt.plot(iterations, params_array_walker[param_counter][10*(n)::no_walkers], color = colours[n], label = "walker %s" %(10*(n)) , alpha = 1)

#PLOT MEAN & MEDIAN TRACE PLOTS
def plot_mean_median_trace_plots(base_filename, hdu_suffixes, params, start_timer, save_filename = False, plot_in_console = False):
    
    #SET UP FIGURE SIZES AND THE NUMBER OF PLOTS NEEDED
    plt.figure(figsize=(30, 10*len(params)))
    number_plots = len(params)

    #PLOT THE WALKER BULK MEAN & MEDIAN FOR EACH PARAMETER
    for j in range(len(params)):
        plt.subplot(number_plots*2,2,j+1)
        plot_timer = datetime.datetime.now() #ADDITIONAL TIMER FOR PLOT TIMINGS

        #PLOT FOR EACH BATCH
        for i in range(len(hdu_suffixes)):
            hdul = fits.open('%s%s.fits' %(base_filename, hdu_suffixes[i]))
            cols = hdul[1].columns
            data = hdul[1].data
            
            batch_iterations = [1 + (i*1e4), 1e4 + (i*1e4)]
            
            #CALCUALTE CENTILES
            q16 = np.quantile(data[params[j]], 0.16)
            q50 = np.quantile(data[params[j]], 0.50)
            q84 = np.quantile(data[params[j]], 0.84)
            
            #CALCUALTE MEAN & STANDARD DEVIATION
            mean = np.mean(data[params[j]])
            stan_dev = np.std(data[params[j]])
            
            #PLOT CENTILES
            plt.plot(batch_iterations, [q50,q50], "g-")
            plt.plot(batch_iterations, [q16,q16], "g-")
            plt.plot(batch_iterations, [q84,q84], "g-")   
            
            #PLOT MEAN & STANDARD DEVIATION
            plt.plot(batch_iterations, [mean,mean], "k-")
            plt.plot(batch_iterations, [mean + stan_dev,mean + stan_dev], "r-")
            plt.plot(batch_iterations, [mean - stan_dev,mean - stan_dev], "r-")
            
        #PLOTTING PARAMS
        plt.xlim(1, len(hdu_suffixes)*1e4)
        plt.xlabel("iterations", fontsize = 14)
        plt.ylabel(params[j], fontsize = 14)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)

        print('PARAM: %s, PLOT TIME: %s' %(params[j], datetime.datetime.now() - plot_timer))

    #INCREASE SIZE BETWEEN PLOTS
    plt.tight_layout(pad=3.0)

    #SAVE PLOT TO A PNG FILE
    if (save_filename != False):
        save_timer = datetime.datetime.now() #ADDITIONAL TIMER FOR PLOT SAVING
        plt.savefig("%s_mean_stan_dev_trace_plots.png" %save_filename, bbox_inches = "tight")
        print("\nSAVE TIME: %s" %(datetime.datetime.now() - save_timer))

    #PLOT IN CONSOLE
    if (plot_in_console):
        plt.show()
        
    print("\nMEAN & STANDARD DEVIATION TRACE PLOTS CREATED, TIME:", datetime.datetime.now() - start_timer)  
    print("--------------------------------------------\n") 

#BRINGS TOGETHER OTHER FUNCTIONS TO PLOT TRACE PLOTS
def plotter(no_batches_to_plot, hdu_list, base_filename, start_timer, no_walkers, 
            thin_factor, iterations_per_batch, save_filename, plot_in_console, 
            no_individual_walkers_to_plot, plot_trace = False, 
            plot_mean_median_trace = False):
    
    skipped_batches, hdu_suffixes = length_checking(no_batches_to_plot, hdu_list, base_filename, start_timer)
    params = get_params(hdu_list[0], base_filename, start_timer)

    ####################################
    #PLOT TRACE PLOTS FOR EACH PARAMETER
    ####################################
    if (plot_trace):   
        params_array_walker = extract_trace_plots_data(hdu_suffixes, base_filename, params, start_timer)
        plot_trace_plots(params_array_walker, params, start_timer, no_walkers, 
                         thin_factor, iterations_per_batch,
                         save_filename, plot_in_console,
                         no_individual_walkers_to_plot)

    ###################################################
    # PLOT MEAN & MEDIAN TRACE PLOTS FOR EACH PARAMETER
    ###################################################
    if (plot_mean_median_trace):
        plot_mean_median_trace_plots(base_filename, hdu_suffixes, params, start_timer, save_filename, plot_in_console)


#-----------------------------------------
#MAIN SCRIPT

print("--------------------------------------------") 

#SET MAIN PARAMETERS
base_filename = "rev0956_27_9_23_afree_2M"
save_filename = 'GRO1_afree_final_ten'

no_walkers = 200
iterations_per_batch = 10000
thin_factor = 20
no_batches_to_plot = 10 #CHANGE TO 'All' TO PLOT EVERYTHING (MIGHT BE VERY SLOW OR CRASH)
no_individual_walkers_to_plot = 4
plot_in_console = True

#THIS WOULD HAVE TO BE MANUALLY SET PER CHAIN TO CAPTURE ALL THE FILENAMES
#THE BELOW EXAMPLE IS FOR GRO1
hdu_list = [
        ""    
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

#PLOT TRACE PLOTS
plotter(no_batches_to_plot, hdu_list, base_filename, start_timer, no_walkers, 
            thin_factor, iterations_per_batch, save_filename, plot_in_console, 
            no_individual_walkers_to_plot, 
            plot_trace = True, 
            plot_mean_median_trace = True)

#WRAPPING UP
print("TOTAL TIME: %s" %(datetime.datetime.now() - start_timer))