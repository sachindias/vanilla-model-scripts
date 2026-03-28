import matplotlib.pyplot as plt
import datetime
import numpy as np
from astropy.io import fits
import datetime

#-----------------------------------------
#FUNCTIONS

#CHECK WHICH PARAMETERS ARE NEEDED FOR PLOTTING
def check_params_present(MCMC_file_location, MCMC_base_filenames):
    
    #ALL FITS FILES SHOULD HAVE THE SAME FIELDS
    hdul = fits.open('%s/%s1.fits' %(MCMC_file_location, MCMC_base_filenames))
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
def get_params(MCMC_file_location, MCMC_base_filenames, start_timer):
    
    #PARAMETERS THAT MIGHT BE MISSING
    checked_params = check_params_present(MCMC_file_location, MCMC_base_filenames)
    
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
    params = []
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

#EXTRACTS THE DATA FROM EACH OF THE MCMC SAMPLES AND PUTS INTO AN ARRAY PER PARAMETER
def extract_MCMC_data(MCMC_file_location, MCMC_base_filenames, params, start_timer, no_walkers, scale_down_factor):

    i = 0
    MCMC_values = [[]] * len(params) #EMPTY ARRAY FOR ALL MCMC PARAMETER VALUES
    
    #EXTRACTS DATA FOR EACH FIELD AND EACH SAMPLE AND ADDS TO THE ARRAY
    while (i < len(params)):
        for n in range(3):   
            hdul = fits.open('%s/%s%s.fits' %(MCMC_file_location, MCMC_base_filenames, n+1))
            cols = hdul[1].columns
            data = hdul[1].data
            
            param_array = data[params[i]]
            
            if (n == 0):
                MCMC_values[i] = param_array
            else:
                MCMC_values[i] = np.concatenate((MCMC_values[i], param_array))
                
        print('PARAM %s: %s, TIME: %s' %(i+1, params[i], datetime.datetime.now() - start_timer))
        i = i + 1

    print("\nMCMC DATA EXTRACTED")  
    print("--------------------------------------------")
    
    thinned_MCMC_values = MCMC_thinning(params, MCMC_values, no_walkers, scale_down_factor, start_timer)
    
    return thinned_MCMC_values

#THINS THE MCMC CHAIN PER WALKER
def MCMC_thinning(params, MCMC_values, no_walkers, scale_down_factor, start_timer):
    thinned_MCMC_values = [[]] * len(params)
    
    for n in range(len(params)):
        for k in range(no_walkers):
            param_values_thinned = MCMC_values[n][k::scale_down_factor]
            
            if (n == 0) & (k == 0):
                thinned_MCMC_values[n] = param_values_thinned
            else:
                thinned_MCMC_values[n] = np.concatenate((thinned_MCMC_values[n], param_values_thinned))
                
        print('PARAM %s: %s, TIME: %s' %(n + 1, params[n], datetime.datetime.now() - start_timer))
    print("\nMCMC DATA THINNED")  
    print("--------------------------------------------")
    
    return thinned_MCMC_values

#EXTRACTS THE DATA FROM THE NS FILE AND PUTS INTO AN ARRAY PER PARAMETER
def extract_NS_data(NS_file_location, NS_filename, params):
    Nested_Sampling = open("%s/%s.txt" %(NS_file_location, NS_filename), "r")
    lines = Nested_Sampling.readlines()
    Nested_Sampling.close()
    
    #IGNORE FIRST LINE AND SORT VALUES INTO THE ARRAY
    print("IGNORED FIRST LINE:", lines[0])
    NS_values = [[] for _ in range(len(params))] #EMPTY ARRAY FOR ALL NS PARAMETER VALUES
    for line in lines[1:]:  # skip first line
        values = line.split()
        for n in range(len(params)):
            NS_values[n].append(float(values[n]))

    NS_values = unlog_NS_Mdd_and_norm_R(params, NS_values, lines)

    print("\nNS DATA EXTRACTED")  
    print("--------------------------------------------\n")

    return NS_values

#UNLOGS NS PARAMETERS TO BRING IN LINE WITH MCMC VALUES
def unlog_NS_Mdd_and_norm_R(params, NS_values, lines):
    #CHECK WHICH POSITIONS Mdd & norm_R ARE
    Mdd_pos = np.where(np.array(params) == "Mdd__40")[0][0] 
    norm_R_pos = np.where(np.array(params) == "norm__60")[0][0]
    
    #UNLOG PARAMETERS
    for n in range(len(NS_values[0])):
        NS_values[Mdd_pos][n] = 10**float(NS_values[Mdd_pos][n])
        NS_values[norm_R_pos][n] = 10**float(NS_values[norm_R_pos][n])

    return NS_values


#-----------------------------------------
#MAIN EXAMPLE SCRIPT

#TIMER
start_timer = datetime.datetime.now()
print(start_timer)
print('\n')

rev = "0956"
spin = "afree"

step_sizes = [
            200,
            250,
            300, 
            ]

#MCMC
MCMC_base_filenames = "rev0956_afree_2M_SAMPLE"
MCMC_file_location = "MCMC"

params = get_params(MCMC_file_location, MCMC_base_filenames, start_timer)
no_plots = len(params)
plt.figure(figsize=(15, 12*no_plots))

MCMC_values = extract_MCMC_data(MCMC_file_location, MCMC_base_filenames, params, start_timer, 200, 600)
        
for n in range(len(MCMC_values)):
    plt.subplot(no_plots*2,2,n+1)
    plt.hist(MCMC_values[n], bins = 30, histtype="step", color = "k", linestyle = "--", linewidth = 2, density = True, label = "MCMC", alpha = 0.7) 
    
    #PLOT SETTINGS
    plt.title(params[n], fontsize = 16)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)

#NESTED SAMPLING
colours = [
    "tab:blue", "tab:red","tab:green", 
    "tab:orange", "tab:purple", "tab:brown", 
    "tab:pink", "tab:gray", "tab:olive", "tab:cyan"
    ]

NS_filename = "equal_weighted_post"

for i in range(len(step_sizes)):
    print("\n%s STEPS" %step_sizes[i])
    NS_file_location = ("GRO_J1655_40_EPIC_rev%s_%ssteps_%s_fit/chains" %(rev, step_sizes[i], spin))    
    NS_values = extract_NS_data(NS_file_location, NS_filename, params)
    
    for n in range(len(MCMC_values)):
        plt.subplot(no_plots*2,2,n+1)
        plt.hist(NS_values[n], bins = 30, histtype="step", color = colours[i], linestyle = "-", linewidth = 2, density = True, label = "NS %s steps" %step_sizes[i], alpha = 0.9)
        plt.legend()  

plt.savefig("GRO_J1655_40_rev%s_%s_StepSize_Param_Comparison.png" %(rev, spin), bbox_inches = "tight")
plt.show()

print("TOTAL TIME: %s" %(datetime.datetime.now() - start_timer))