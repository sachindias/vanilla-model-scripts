import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import datetime
import corner
from matplotlib import rcParams

start = datetime.datetime.now()
print(start)
print('\n')

#SET UP
revs = [
        "0956",
        "0964a",
        "0964b",
        "0965",
        "0966",
        "0970"
        ]

no_walkers = 200
scale_down_factor = 600
skip_factor = no_walkers*scale_down_factor

#GO THROUGH EACH OBS
for n in range(len(revs)):
    big_array = np.array([])
    for i in range(3):
        
        #OPEN FILES TO USE IN THE CODE
        hdul = fits.open('rev%s_afree_2M_SAMPLE%s.fits' %(revs[n], i+1))
        cols = hdul[1].columns
        data = hdul[1].data
        
        FS_array = data['FIT_STATISTIC']
        big_array = np.concatenate((big_array, FS_array))
    
    #SORT INTO WALKERS & THIN
    thin_array = np.array([])
    for k in range(no_walkers):
        thin_array = np.concatenate((thin_array, big_array[k::skip_factor]))
    
    print("rev%s: Avg Chi-sq = %.2f" %(revs[n], np.mean(thin_array)))

hdul.close()
 
end = datetime.datetime.now()
print('\n')
print(end - start)
