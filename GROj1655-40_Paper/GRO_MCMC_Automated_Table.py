from astropy.io import fits
import numpy as np
import datetime

start = datetime.datetime.now()
print(start)

#CALCULATES 16th, 50th and 84th CENTILES
def quantiles(data):     
        
        med = np.median(data)
        err_up = np.quantile(data, 0.84) - med
        err_down = med - np.quantile(data, 0.16)
        
        #MULTIPLIER FOR norm_R
        if (j == 19):
            med = med*10**4
            err_up = err_up*10**4
            err_down = err_down*10**4
        
        return med, err_up, err_down
   
#ROUNDS CENTILES TO 3 DECIMAL PLACES
def rounder(med, err_up, err_low):
    
    str_med = str(round(med, 3))
    str_eu = str(round(err_up, 3))
    str_ed = str(round(err_low, 3))
    
    length_m = len(str_med) - str_med.find(".")
    if (length_m < 4):
        str_med = str_med + "0"
        
    length_eu = len(str_eu) - str_eu.find(".")
    if (length_eu < 4):
        str_eu = str_eu + "0"
        
    length_ed =len(str_ed) - str_ed.find(".")
    if (length_ed < 4):
        str_ed = str_ed + "0"
        
    return str_med, str_eu, str_ed

#########################
#MAIN CODE###############
#########################

no_walkers = 200
scale_down_factor = 600
skip_factor = no_walkers*scale_down_factor

REVS = [
        "rev0956",
        "rev0964a",
        "rev0964b",
        "rev0965",
        "rev0966",
        "rev0970",
        ]

Overleaf_params = [
        "$N_H$",
	    "$N_{O_I}$",
	    "$N_{O_{II}}$",
	    "$N_{Ne_{I}}$",
	    "$N_{Ne_{II}}$",
	    "$N_{Fe}$",
        
        "$\Gamma_{S}$",
        "$f_{scat}$",
        
        "$a$",
    	"$i$",
    	"$M$",
   		r"$\dot{M}$",
    	"$D$",
    	"$f_{col}$",
        
        "$I_1$",
        "$I_2$",
	    "$\Gamma_{R}$",
	    r"$\log(\xi)$",
	    "$A_{Fe}$",
	    r"$norm_R (\times 10^{-4})$",
        
        r"$E_1$",
        r"$\tau_1$",
        r"$\tau_2$",
        r"$\tau_3$",
        
        r"$E_c$",
        r"$\tau_{max}$",
        r"$W$",
        ]

big_array = np.array([None]*len(Overleaf_params))

big_array[0] = r"\textsc{ISMabs}"
big_array[6] = r"\textsc{SIMPL}"
big_array[8] = r"\textsc{KERRBB}"
big_array[14] = r"\textsc{relxillCp}"
big_array[20] = r"\textsc{gabs$_1$}"
big_array[22] = r"\textsc{gabs$_2$}"
big_array[23] = r"\textsc{gabs$_3$}"
big_array[24] = r"\textsc{smedge}"

for n in range(len(big_array)):
    try:
        if (len(big_array[n]) > 1):
            big_array[n] = big_array[n] + " & "
    except Exception:
        big_array[n] = "& "
        
    big_array[n] = big_array[n] + Overleaf_params[n] + " & "

#########################
#PARAMS##################
#########################

for n in range(len(REVS)):
    print("ON FILE:", n+1)

    params = ["H__2", 
              "O_I__10", "O_II__11",
              "Ne_I__13", "Ne_II__14",
              "Fe__31",
              "Gamma__33","FracSctr__34",
              "a__37",
              "i__38","Mbh__39","Mdd__40","Dbh__41", "hd__42",
              "Index1__51", "Index2__52",
              "gamma__54", "logxi__55", "Afe__57",
              "norm__60",
              ]

    #EXTRACTING DATA AND THINNING
    for j in range(20):
        for k in range(3):
    
            hdul = fits.open('MCMC/%s_afree_2M_SAMPLE%s.fits' %(REVS[n],k+1))
            cols = hdul[1].columns
            data = hdul[1].data
            
            if (k == 0):
                temp_param_array = data[params[j]]   
            else:
                temp_param_array = np.concatenate((temp_param_array, data[params[j]]), axis=0)
    
        temp_array = np.array([])
        for k in range(no_walkers):
            temp_array = np.concatenate((temp_array, temp_param_array[k::skip_factor]))

        med, err_up, err_low = quantiles(temp_array)
        str_med, str_eu, str_ed = rounder(med, err_up, err_low)
        big_array[j] = big_array[j] + str_med + "$^{+" + str_eu + "}_{-" + str_ed + "}$" + " & "  

#########################
#GABS####################
#########################

    #FOR MODELS WITHOUT A GABS COMPONENT
    if (n < 1):
        big_array[20] = big_array[20] + "-" + " & "
        big_array[21] = big_array[21] + "-" + " & "
        big_array[22] = big_array[22] + "-" + " & "
        big_array[23] = big_array[23] + "-" + " & "
        
    #FOR MODELS WITH A GABS COMPONENT
    if (n > 0):
        params = params + ['LineE__61', 'Strength__63', 'Strength__66', 'Strength__69']
    
        for j in range(20,24):
            for k in range(3):
        
                hdul = fits.open('MCMC/%s_afree_2M_SAMPLE%s.fits' %(REVS[n],k+1))
                cols = hdul[1].columns
                data = hdul[1].data
                
                if (k == 0):
                    temp_param_array = data[params[j]]   
                else:
                    temp_param_array = np.concatenate((temp_param_array, data[params[j]]), axis=0)
        
            temp_array = np.array([])
            for k in range(no_walkers):
                temp_array = np.concatenate((temp_array, temp_param_array[k::skip_factor]))

            med, err_up, err_low = quantiles(temp_array)
            str_med, str_eu, str_ed = rounder(med, err_up, err_low)
            big_array[j] = big_array[j] + str_med + "$^{+" + str_eu + "}_{-" + str_ed + "}$" + " & "  
 

#########################
#SMEDGE##################
#########################
   
    #FOR MODELS WITHOUT A SMEDGE COMPONENT
    if (n < 4):
        big_array[24] = big_array[24] + "-" + " & "
        big_array[25] = big_array[25] + "-" + " & "
        big_array[26] = big_array[26] + "-" + " & "

    #FOR MODELS WITH A SMEDGE COMPONENT
    if (n > 3):
        params = params + ['edgeE__73', 'MaxTau__74', 'width__76']
           
        for j in range(24,27):
            for k in range(3):
        
                hdul = fits.open('MCMC/%s_afree_2M_SAMPLE%s.fits' %(REVS[n],k+1))
                cols = hdul[1].columns
                data = hdul[1].data
                
                if (k == 0):
                    temp_param_array = data[params[j]]   
                else:
                    temp_param_array = np.concatenate((temp_param_array, data[params[j]]), axis=0)
        
            temp_array = np.array([])
            for k in range(no_walkers):
                temp_array = np.concatenate((temp_array, temp_param_array[k::skip_factor]))

            med, err_up, err_low = quantiles(temp_array)
            str_med, str_eu, str_ed = rounder(med, err_up, err_low)
            big_array[j] = big_array[j] + str_med + "$^{+" + str_eu + "}_{-" + str_ed + "}$" + " & "  
         
#PRINTS OUT TABLE FOR LaTeX USE
print("\n")
print(r"\hline")
print(r"\hline")    
hlines = np.array([6, 8, 14, 20, 22, 23, 24])  
for n in range(len(big_array)):
    if (len(np.where(hlines == n)[0]) > 0):
        print(r"\hline")
    print(big_array[n][:-3] + r"\\")
print(r"\hline")
print(r"\hline")

print("\nTOTAL TIME: ", datetime.datetime.now() - start)    