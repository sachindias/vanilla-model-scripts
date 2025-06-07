import numpy as np
from astropy.io import fits
import datetime

start = datetime.datetime.now()
print(start)
print('\n')

#OPEN FILES TO USE IN CODE
hdul = fits.open('<FILENAME>.fits')
cols = hdul[1].columns
data = hdul[1].data

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
    a = 1
    print("SPIN FREE")
except Exception:
    a = 0
    
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


#LOAD IN NOT ALREADY LOADED DATA
data1 = []

H_array = data['H__2']
OI_array = data['O_I__10']
OII_array = data['O_II__11']
NeI_array = data['Ne_I__13']
NeII_array = data['Ne_II__14']
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

FS_array = data['FIT_STATISTIC']
length = len(FS_array)

data1.append([H_array, OI_array, OII_array, 
              NeI_array, NeII_array, Fe_array,
              Gamma_K_array, FrSc_array, 
              i_array, Mbh_array, Mdd_array, Dbh_array, hd_array, 
              I1_array, I2_array, 
              Gamma_R_array, logxi_array, Afe_array, norm_array])

param_index = [2,10,11,13,14,31,
               33,34,38,39,40,41,42,
               51,52,54,55,57,60]

data1 = data1[0]

if (N_III == 1):
    data1 = data1[:3] + [OIII_array] + data1[3:5] + [NeIII_array] + data1[5:]
    param_index = param_index[:3] + [12] + param_index[3:5] + [15] + param_index[5:]

if (a == 1):
    loc = np.where(np.array(param_index) == 34)[0][0] + 1
    data1 = data1[:loc] + [a_array] + data1[loc:] 
    param_index = param_index[:loc] + [37] + param_index[loc:] 

if (gabs == 1):
    data1 = data1 + [LineE_array, Strength_array1, Strength_array2, Strength_array3]
    param_index = param_index + [61, 63, 66, 69]

if (smedge == 1):
    data1 = data1 + [Smedge_E_array, Smedge_Tau_array, Smedge_w_array]
    param_index = param_index + [73,74,76]

#PRINT VALUE OF MIN FIT STAT. VALUE
min_FS = np.min(FS_array)
min_FS_pos = np.where(FS_array == min_FS)
min_FS_pos = min_FS_pos[0][0]
print("FS: %.2f @ n = %s \n" %(min_FS,min_FS_pos))

#PRINT newpar FOR XSPEC
#THESE ARE ALSO THE PARAMETERS OF THE BEST FIT FROM MCMC
for n in range(len(data1)):
    print("newpar %s %s" %(param_index[n],data1[n][min_FS_pos]))

hdul.close()

end = datetime.datetime.now()
print('\n')
print(end - start)
