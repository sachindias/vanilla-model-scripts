from astropy.io import fits
import numpy as np
import datetime

#TIMER
start = datetime.datetime.now()
print(start)

#CALCULATES 16th, 50th and 84th CENTILES
def quantiles(data):     
        
    if (j == 11 or j == 19):
        data = 10**data

    med = np.median(data)
    err_up = np.quantile(data, 0.84) - med
    err_down = med - np.quantile(data, 0.16)
    
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
    
    #ENSURES CONSITENT SIGNIFICANT FIGURES
    length_m = len(str_med) - str_med.find(".")
    if (length_m < 4):
        str_med = str_med + "0"
        
    length_eu = len(str_eu) - str_eu.find(".")
    if (length_eu < 4):
        str_eu = str_eu + "0"
        
    length_ed =len(str_ed) - str_ed.find(".")
    if (length_ed < 4):
        str_ed = str_ed + "0"
    
    #FIXES WHEN THE MEDIAN OR ERRORS SIT AT THE MODEL BOUNDARY 
    if (med == 85.00):
        str_med = "85.000"
    if (err_up == 0.0):
        str_eu = "X"
    if (err_low == 0.0):       
        str_ed = "X"
        
    #FIX WHEN MEDIAN IS VERY SMALL
    if (med < 0.001):
        str_med = str(round(med, 4))
        str_eu = str(round(err_up, 4))
        str_ed = str(round(err_low, 4))
        
        length_eu = len(str_eu) - str_eu.find(".")
        if (length_eu < 5):
            str_eu = str_eu + "0"
    
    return str_med, str_eu, str_ed

##########
#MAIN CODE
##########

FILES = [
        "NS/equal_weighted_post_rev0956",
        "NS/equal_weighted_post_rev0964a",
        "NS/equal_weighted_post_rev0964b",
        "NS/equal_weighted_post_rev0965",
        "NS/equal_weighted_post_rev0966",
        "NS/equal_weighted_post_rev0970",
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

table_rows = np.array([None]*len(Overleaf_params))

#ADDS MODEL COMPONENTS & '&'S
table_rows[0] = r"\textsc{ISMabs}"
table_rows[6] = r"\textsc{SIMPL}"
table_rows[8] = r"\textsc{KERRBB}"
table_rows[14] = r"\textsc{relxillCp}"
table_rows[20] = r"\textsc{gabs$_1$}"
table_rows[22] = r"\textsc{gabs$_2$}"
table_rows[23] = r"\textsc{gabs$_3$}"
table_rows[24] = r"\textsc{smedge}"

for n in range(len(table_rows)):
    try:
        if (len(table_rows[n]) > 1):
            table_rows[n] = table_rows[n] + " & "
    except Exception:
        table_rows[n] = "& "
        
    table_rows[n] = table_rows[n] + Overleaf_params[n] + " & "
    
#######
#PARAMS
#######   

for n in range(len(FILES)):   
    print("ON FILE:", n+1, ":", FILES[n])

    #OPEN FILES
    Nested_Sampling = open("%s.txt" %FILES[n], "r")
    lines = Nested_Sampling.readlines()
    Nested_Sampling.close()
    
    NS_array = np.array([None]*len(Overleaf_params))
    
    #CONVERTS DATA INTO USABLE ARRAY
    for h in range(len(lines[0].split(" "))):
        temp = []
        for k in range(len(lines)):
            if (k > 0):
                temp.append(float(lines[k].split(" ")[h]))
        temp = np.array(temp)
        NS_array[h] = temp

    #APPEND VALUES TO LIST OF ALL TABLE ROWS
    for j in range(20):
        
        med, err_up, err_low = quantiles(NS_array[j])
        str_med, str_eu, str_ed = rounder(med, err_up, err_low)
        part_1 = str_med
        if (str_eu == "X"):
            part_2 = ""
        else:
            part_2 = "$^{+" + str_eu + "}"
        if (str_ed == "X"):
             part_3 = ""
        else:
             part_3 = "_{-" + str_ed + "}$" 
             
        if (part_2 == "" and part_3 != ""):
            table_rows[j] = table_rows[j] + part_1 + "$" + part_3 + " & "  
                
        else:
            table_rows[j] = table_rows[j] + part_1 + part_2 + part_3 + " & "  
    
    #####
    #GABS
    #####


    if (n < 1):
        table_rows[20] = table_rows[20] + "-" + " & "
        table_rows[21] = table_rows[21] + "-" + " & "
        table_rows[22] = table_rows[22] + "-" + " & "
        table_rows[23] = table_rows[23] + "-" + " & "
        
    if (n > 0):
        for j in range(20,24):
        
            #APPEND VALUES TO LIST OF ALL TABLE ROWS
            med, err_up, err_low = quantiles(NS_array[j])
            str_med, str_eu, str_ed = rounder(med, err_up, err_low)
            table_rows[j] = table_rows[j] + str_med + "$^{+" + str_eu + "}_{-" + str_ed + "}$" + " & "  

    #######
    #SMEDGE
    #######
   
    if (n < 4):
        table_rows[24] = table_rows[24] + "-" + " & "
        table_rows[25] = table_rows[25] + "-" + " & "
        table_rows[26] = table_rows[26] + "-" + " & "

    
    if (n > 3):           
        for j in range(24,27):
            
            #APPEND VALUES TO LIST OF ALL TABLE ROWS
            med, err_up, err_low = quantiles(NS_array[j])
            str_med, str_eu, str_ed = rounder(med, err_up, err_low)
            table_rows[j] = table_rows[j] + str_med + "$^{+" + str_eu + "}_{-" + str_ed + "}$" + " & "   
    
#PRINT THE TOTAL LIST OF STRINGS TO BUILD THE TABLE
print("\n")
print(r"\hline")
print(r"\hline")    
hlines = np.array([6, 8, 14, 20, 22, 23, 24])  
for n in range(len(table_rows)):
    if (len(np.where(hlines == n)[0]) > 0):
        print(r"\hline")
    print(table_rows[n][:-3] + r"\\")
print(r"\hline")
print(r"\hline")

#PRINT RUNTIME
print("\nTOTAL TIME: ", datetime.datetime.now() - start) 