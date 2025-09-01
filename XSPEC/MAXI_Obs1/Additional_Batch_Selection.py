import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import datetime
import corner
from matplotlib import rcParams
import shutil

start = datetime.datetime.now()
print(start)
print('\n')

#SET THESE NUMBERS!!!!!!
no_walkers = 200
multi = 1 #1 FOR MULTIPLE FILES, 0 FOR ONE FILE

def new_filename(BASE_FILE):
    try:
        if ("ext" in BASE_FILE):
            NEW_FILE = BASE_FILE.split("ext")[0] + "ext" + str(int(BASE_FILE.split("ext")[1]) + 1)
        else:
            NEW_FILE = BASE_FILE + "ext1"
    
        print(BASE_FILE)
        print(NEW_FILE)
        print(COPY_FILE, "\n")
        
    except Exception:
        print("ERROR 1: NAMING\n")    

    return NEW_FILE

def copy_file(COPY_FILE, NEW_FILE):
    try:
        shutil.copy("%s.fits" %COPY_FILE, "%s.fits" %NEW_FILE)
    except Exception:
        print("ERROR 2: COPYING\n")    
        
    return NEW_FILE

def process(BASE_FILE,COPY_FILE, no_walkers):
    #LOAD BASE FITS FILES
    hdul = fits.open('%s.fits' %BASE_FILE)
    cols = hdul[1].columns
    data = hdul[1].data
    
    NEW_FILENAME = new_filename(BASE_FILE)
    NEW_FILE = copy_file(COPY_FILE, NEW_FILENAME)
    
    #TAKE FITS FILE OF 200 AND REWRITE VALUES WITHIN CORRECT REGIONS
    hdulI = fits.open("%s.fits" %NEW_FILE, mode = 'update')
    colsI = hdulI[1].columns
    dataI = hdulI[1].data
    
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
    
    #TEST FOR CERTAIN FREE PARAMETERS IN MODEL
    #AND EXTRACT THE LAST POINT FOR EACH WALKER
    
    #TEST FOR N_O_III AND N_NE_III
    try:
        OIII_array = data['O_III__12'][-200:]
        NeIII_array = data['Ne_III__15'][-200:]
        N_III = 1
        print("N_O_III AND N_NE_III INCL.")
    except Exception:
        N_III = 0
        
    #TEST IF SPIN FREE
    try:
        a_array = data["a__37"][-200:]
        spin = 1
        print("SPIN FREE")
    except Exception:
        spin = 0
        
    #CHECK FOR GABS
    try:
        LineE_array = data['LineE__61'][-200:]
        Strength_array1 = data['Strength__63'][-200:]
        Strength_array2 = data['Strength__66'][-200:]
        Strength_array3 = data['Strength__69'][-200:]
        gabs = 1
        print("GABS INCL.")
    except Exception:
        gabs = 0
        
    #CHECK FOR SMEDGE
    try:
        Smedge_E_array = data['edgeE__73'][-200:]
        Smedge_Tau_array = data['MaxTau__74'][-200:]
        Smedge_w_array = data['width__76'][-200:]
        smedge = 1
        print("SMEDGE INCL.")
    except Exception:
        smedge = 0
        
    #CHECK FOR NUSTAR DATA
    try:
        factor_1_array = data['factor__61'][-200:]
        factor_2_array = data['factor__121'][-200:]
        NuSTAR_data = 1
        print("NuSTAR DATA INCL.")
    except Exception:
        NuSTAR_data = 0
        print("ONLY XMM DATA")
    
    #EXTRACT DATA FROM FIRST FILE (IF NOT ALREADY)
    H_array = data['H__2'][-200:]
    OI_array = data['O_I__10'][-200:]
    OII_array = data['O_II__11'][-200:]
    NeI_array = data['Ne_I__13'][-200:]
    NeII_array = data['Ne_II__14'][-200:]
    Fe_array = data['Fe__31'][-200:]       
    
    Gamma_K_array = data['Gamma__33'][-200:]
    FrSc_array = data['FracSctr__34'][-200:]
    
    i_array = data['i__38'][-200:]
    Mbh_array = data['Mbh__39'][-200:]
    Mdd_array = data['Mdd__40'][-200:]
    Dbh_array = data['Dbh__41'][-200:]
    hd_array = data['hd__42'][-200:]
    
    I1_array = data['Index1__51'][-200:]
    I2_array = data['Index2__52'][-200:]
    Gamma_R_array = data['gamma__54'][-200:]
    logxi_array = data['logxi__55'][-200:]
    Afe_array = data['Afe__57'][-200:]
    norm_array = data['norm__60'][-200:]
    
    FS = data['FIT_STATISTIC'][-200:]
    
    #REPACKAGING DATA
    dataI['H__2'] = H_array
    dataI['O_I__10'] = OI_array
    dataI['O_II__11'] = OII_array
    dataI['Ne_I__13'] = NeI_array
    dataI['Ne_II__14'] = NeII_array
    dataI['Fe__31'] = Fe_array   
    
    if (N_III == 1): 
        dataI['O_III__12'] = OIII_array
        dataI['Ne_III__15'] = NeIII_array
    
    dataI['Gamma__33'] = Gamma_K_array
    dataI['FracSctr__34'] = FrSc_array
    
    if (spin == 1):
        dataI['a__37'] = a_array
    
    dataI['i__38'] = i_array
    dataI['Mbh__39'] = Mbh_array
    dataI['Mdd__40'] = Mdd_array
    dataI['Dbh__41'] = Dbh_array
    dataI['hd__42'] = hd_array
    
    dataI['Index1__51'] = I1_array
    dataI['Index2__52'] = I2_array
    dataI['gamma__54'] = Gamma_R_array
    dataI['logxi__55'] = logxi_array
    dataI['Afe__57'] = Afe_array
    dataI['norm__60'] = norm_array  
    
    if (gabs == 1):
        dataI['LineE__61'] = LineE_array
        dataI['Strength__63'] = Strength_array1
        dataI['Strength__66'] = Strength_array2
        dataI['Strength__69'] = Strength_array3
                        
    if (smedge == 1):
        dataI['edgeE__73'] = Smedge_E_array
        dataI['MaxTau__74'] = Smedge_Tau_array
        dataI['width__76'] = Smedge_w_array
                    
    if (NuSTAR_data == 1):   
        dataI['factor__61'] = factor_1_array
        dataI['factor__121'] = factor_2_array
    
    dataI['FIT_STATISTIC'] = FS
    
    hdulI.flush()
    hdulI.close()    

if (multi == 1):
    
    #MAKE A COPY OF THE REWRITABLE FILE
    BASE_FILE_ALL = ["rev3531_1_9_22_a0_2M",
                     "rev3531_1_9_22_amax_2M_P2",
                     "rev3531_1_9_22_afree_2M_P2_ext1"]

    COPY_FILE_ALL = ["rev3531_1_9_22_a0_rewritable_200",
                    "rev3531_1_9_22_amax_rewritable_200",
                    "rev3531_1_9_22_afree_rewritable_200"]

    for CP_COUNT in range(len(BASE_FILE_ALL)):

        BASE_FILE = BASE_FILE_ALL[CP_COUNT]
        COPY_FILE = COPY_FILE_ALL[CP_COUNT]

        process(BASE_FILE,COPY_FILE, no_walkers)

elif (multi == 0):
    
    BASE_FILE = "rev3531_1_9_22_a0_2M"
    COPY_FILE = "rev3531_1_9_22_a0_rewritable_200"
    
    process(BASE_FILE,COPY_FILE, no_walkers)

end = datetime.datetime.now()
print(end - start)