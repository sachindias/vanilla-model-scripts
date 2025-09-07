import numpy as np
from astropy.io import fits
import datetime

start = datetime.datetime.now()
print(start)
print('\n')

#LOAD INITIAL CHAIN POINTS FITS FILE
hdul = fits.open('ICP_rev3531_1_9_22_a0.fits') 
cols = hdul[1].columns
data = hdul[1].data

#LOAD FIRST RUN FITS FILE TO BE OVERWRITTEN
hdulI = fits.open('rev3531_1_9_22_a0_2M.fits', mode = 'update')
colsI = hdulI[1].columns
dataI = hdulI[1].data

#SET THESE NUMBERS!!!!!!
M_dd_Edd_max = 23.2539
M_dd_Edd_min = 18.5592
no_walkers = 200
MAXI_GRO = "MAXI"

#MAXI or GRO
if (MAXI_GRO == "MAXI"):
    R_norm = 1.67457e-3
    R_norm_max = 100 * R_norm
    R_norm_min = 1e-6 * R_norm_max
elif (MAXI_GRO == "GRO"):
    R_norm_max = 2.0
    R_norm_min = 1e-9

#TEST FOR CERTAIN FREE PARAMETERS IN MODEL

#TEST FOR N_O_III AND N_NE_III
try:
    OIII_array = data['O_III__12']
    NeIII_array = data['Ne_III__15']
    N_III = 1
    print("N_O_III AND N_NE_III INCL.")
except Exception:
    N_III = 0
    
#TEST IF SPIN FREE
try:
    a_array = data["a__37"]
    spin = 1
    print("SPIN FREE")
except Exception:
    spin = 0
    
#CHECK FOR GABS
try:
    LineE_array = data['LineE__61']
    Strength_array1 = data['Strength__63']
    Strength_array2 = data['Strength__66']
    Strength_array3 = data['Strength__69']
    gabs = 1
    print("GABS INCL.")
except Exception:
    gabs = 0
    
#CHECK FOR SMEDGE
try:
    Smedge_E_array = data['edgeE__73']
    Smedge_Tau_array = data['MaxTau__74']
    Smedge_w_array = data['width__76']
    smedge = 1
    print("SMEDGE INCL.")
except Exception:
    smedge = 0
    
#CHECK FOR NUSTAR DATA
try:
    factor_1_array = data['factor__61']
    factor_2_array = data['factor__121']
    NuSTAR_data = 1
    print("NuSTAR DATA INCL.")
except Exception:
    NuSTAR_data = 0
    print("ONLY XMM DATA")


#EXTRACT DATA FROM ICP FILE (IF NOT ALREADY)
H_array = data['H__2']
OI_array = data['O_I__10']
OII_array = data['O_II__11']
OIII_array = data['O_III__12']
NeI_array = data['Ne_I__13']
NeII_array = data['Ne_II__14']
NeIII_array = data['Ne_III__15']
Fe_array = data['Fe__31']       

Gamma_K_array = data['Gamma__33']
FrSc_array = data['FracSctr__34']

i_array = data['i__38']
Mbh_array = data['Mbh__39']
Mdd_array = data['Mdd__40']
Dbh_array = data['Dbh__41']
hd_array = data['hd__42']

I1_array = data['Index1__51']
I2_array = data['Index2__52']
Gamma_R_array = data['gamma__54']
logxi_array = data['logxi__55']
Afe_array = data['Afe__57']
norm_array = data['norm__60']

FS = data['FIT_STATISTIC']

#EMPTY ARRAYS TO POPULATE WITH 
#ACCEPTED STARTING POINTS
H_a = []
OI_a = []
OII_a = []
OIII_a = []
NeI_a = []
NeII_a = []
NeIII_a = []
Fe_a = []

Gamma_K_a = []
FrSc_a = []

a_a = []
i_a = []
Mbh_a = []
Mdd_a = []
Dbh_a = []
hd_a = []

I1_a = []
I2_a = []
Gamma_a = []
logxi_a = []
Afe_a = []
norm_a = []

LineE_a = []
Strength_1_a = []
Strength_2_a = []
Strength_3_a = []

Smedge_E_a = []
Smedge_Tau_a = []
Smedge_w_a = []

factor_1_a = []
factor_2_a = []

FS_a = []

#UPPER & LOWER LIMITS, AS WELL AS ARRAY OF ALL PARAMETERS
arrays = [Gamma_K_array,FrSc_array,
          i_array,Mbh_array,Mdd_array,Dbh_array,hd_array,
          I1_array,I2_array,
          Gamma_R_array,logxi_array,Afe_array,norm_array]

lower_limits = [
                1, #Gamma_K
                0, #FrSc
                0, #i
                0, #M
                1e-6 * M_dd_Edd_min, #Mdd
                0, #D
                1, #hd
                0, #I1
                0, #I2
                1, #Gamma_R
                0, #log(xi)
                0.5, #A_Fe
                R_norm_min #norm_R
                ]

higher_limits = [
                5, #Gamma_K
                1, #FrSc
                85, #i
                50, #M
                M_dd_Edd_max, #Mdd
                30, #D 
                10, #hd
                10, #I1
                10, #I2
                3.6, #Gamma_R
                4.7, #log(xi)
                10, #A_Fe
                R_norm_max  #norm_R
                ]

#INSERT SPIN PARAMETER INFO
if (spin == 1):
    arrays = arrays[:2] + [a_array] + arrays[2:] 
    higher_limits = higher_limits[:2] + [0.998] + higher_limits[2:] 
    lower_limits = lower_limits[:2] + [-0.998] + lower_limits[2:] 

#INSERT COLUMN DENSITY PARAMETER INFO
N_H_insert_array = [H_array, OI_array, OII_array,
                    NeI_array, NeII_array, Fe_array
                    ]

if (N_III == 1):
    N_H_insert_array = N_H_insert_array[:3] + [OIII_array] + N_H_insert_array[3:5] + [NeIII_array] + N_H_insert_array[5:]

arrays = N_H_insert_array + arrays
n = 0
for n in range(len(N_H_insert_array)):
    higher_limits = np.insert(higher_limits, n, 1e6)
    lower_limits = np.insert(lower_limits, n, 0)
    
#INSERT GABS & SMEDGE PARAMETER INFO    
if (gabs == 1):
    arrays = arrays + [LineE_array, Strength_array1, Strength_array2, Strength_array3]
    higher_limits = higher_limits + [6.6, 0, 0, 0]
    lower_limits = lower_limits + [6.9, 1, 1, 1]

if (smedge == 1):
    arrays = arrays + [Smedge_E_array, Smedge_Tau_array, Smedge_w_array]
    higher_limits = higher_limits + [8, 0, 0.01]
    lower_limits = lower_limits + [9, 1, 1]
 
#INSERT CONSTANTS PARAMETER INFO 
if (NuSTAR_data == 1):
    arrays = arrays + [factor_1_array, factor_2_array]
    higher_limits = higher_limits + [0.5, 0.5]
    lower_limits = lower_limits + [3, 3]
    
     
#LOOKS THROUGH ALL POINTS FROM INITIAL 1000
#PICKS 200 WHICH SATISFY THE LIMITS SET ABOVE
for n in range(len(FS)):
    if (len(Gamma_K_a) < no_walkers):
        True_False = []
        for i in range(len(higher_limits)):
            if (arrays[i][n] > lower_limits[i] and arrays[i][n] < higher_limits[i]):
                True_False.append(True)
            else:
                True_False.append(False)
        if all(True_False):

                #ASSIGNS EACH VALUE OF THE 200 SELECTED 
                #POINTS INTO THE CORRECT ARRAY
                H_a.append(H_array[n])
                OI_a.append(OI_array[n])
                OII_a.append(OII_array[n])
                NeI_a.append(NeI_array[n])
                NeII_a.append(NeII_array[n])
                Fe_a.append(Fe_array[n])
                
                if (N_III == 1):
                    OIII_a.append(OIII_array[n])
                    NeIII_a.append(NeIII_array[n])
                
                Gamma_K_a.append(Gamma_K_array[n])
                FrSc_a.append(FrSc_array[n])
                
                if (spin == 1):
                    a_a.append(a_array[n])
                
                i_a.append(i_array[n])
                Mbh_a.append(Mbh_array[n])
                Mdd_a.append(Mdd_array[n])
                Dbh_a.append(Dbh_array[n])
                hd_a.append(hd_array[n])
                
                I1_a.append(I1_array[n])
                I2_a.append(I2_array[n])
                Gamma_a.append(Gamma_R_array[n])
                logxi_a.append(logxi_array[n])
                Afe_a.append(Afe_array[n])
                norm_a.append(norm_array[n])
                
                if (gabs == 1):
                    LineE_a.append(LineE_array[n])
                    Strength_1_a.append(Strength_array1[n])
                    Strength_2_a.append(Strength_array2[n])
                    Strength_3_a.append(Strength_array3[n])
                    
                if (smedge == 1):
                    Smedge_E_a.append(Smedge_E_array[n])
                    Smedge_Tau_a.append(Smedge_Tau_array[n])
                    Smedge_w_a.append(Smedge_w_array[n])
                
                if (NuSTAR_data == 1):
                    factor_1_a.append(factor_1_array[n])
                    factor_2_a.append(factor_2_array[n])
                
                FS_a.append(FS[n])


#REPACKAGING DATA
dataI['H__2'] = H_a
dataI['O_I__10'] = OI_a
dataI['O_II__11'] = OII_a
dataI['O_III__12'] = OIII_a
dataI['Ne_I__13'] = NeI_a
dataI['Ne_II__14'] = NeII_a
dataI['Ne_III__15'] = NeIII_a
dataI['Fe__31'] = Fe_a   

if (N_III == 1): 
    dataI['O_III__12'] = OIII_a
    dataI['Ne_III__15'] = NeIII_a

dataI['Gamma__33'] = Gamma_K_a 
dataI['FracSctr__34'] = FrSc_a

if (spin == 1):
    dataI['a__37'] = a_a

dataI['i__38'] = i_a
dataI['Mbh__39'] = Mbh_a
dataI['Mdd__40'] = Mdd_a
dataI['Dbh__41'] = Dbh_a
dataI['hd__42'] = hd_a

dataI['Index1__51'] = I1_a
dataI['Index2__52'] = I2_a
dataI['gamma__54'] = Gamma_a
dataI['logxi__55'] = logxi_a
dataI['Afe__57'] = Afe_a
dataI['norm__60'] = norm_a  

if (gabs == 1):
    dataI['LineE__61'] = LineE_a
    dataI['Strength__63'] = Strength_1_a
    dataI['Strength__66'] = Strength_2_a
    dataI['Strength__69'] = Strength_3_a
                    
if (smedge == 1):
    dataI['edgeE__73'] = Smedge_E_a
    dataI['MaxTau__74'] = Smedge_Tau_a
    dataI['width__76'] = Smedge_w_a
                
if (NuSTAR_data == 1):   
    dataI['factor__61'] = factor_1_a
    dataI['factor__121'] = factor_2_a

dataI['FIT_STATISTIC'] = FS_a

hdulI.flush()
hdulI.close()

end = datetime.datetime.now()
print(end - start)
