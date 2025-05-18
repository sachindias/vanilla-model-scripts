########################################
#############IMPORT LIBRARIES###########
########################################
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
from astropy.io import fits
import my_corner as mc

start_timer = datetime.datetime.now()

########################################
############LOADING IN DATA#############
########################################

custom_values = [r"f$_{col}$"]

REVS = [
        "rev0964a",
        "rev0964b",
        "rev0965",
        "rev0966",
        ]

params = [  
    "H__2", 
    "O_I__10", "O_II__11", 
    "Ne_I__13","Ne_II__14",
    "Fe__31",
    "Gamma__33","FracSctr__34",
    "a__37", 
    "i__38","Mbh__39","Mdd__40","Dbh__41", "hd__42",
    "Index1__51", "Index2__52",
    "gamma__54", "logxi__55", "Afe__57",
    "norm__60", 
            ]

total_array = np.array([None]*len(params))
total_array_temp = np.array([None]*len(params))
no_walkers = 200
scale_down_factor = 600
skip_factor = no_walkers*scale_down_factor #scales down the data to a managable amount

print("\nMCMC")
for looper in range(len(REVS)):
    j = looper
    testtime1 = datetime.datetime.now()
    
    for i in range(3):
        
        hdul = fits.open('MCMC/%s_afree_2M_SAMPLE%s.fits' %(REVS[looper],i+1))
        cols = hdul[1].columns
        data = hdul[1].data
        
        for n in range(len(params)):
            temp_param_array = data[params[n]]
            
            if (i == 0):
                total_array_temp[n] = temp_param_array
            else:
                total_array_temp[n] = np.concatenate((total_array_temp[n], temp_param_array), axis=0)
    
    for n in range(len(params)):
        
        temp_param_array = total_array_temp[n]
        temp_array = np.array([])
        for k in range(no_walkers):
            temp_array = np.concatenate((temp_array, temp_param_array[k::skip_factor]))
        
        if (n == 8):
            a_array = temp_array
        if (n == 13):
            hd_array = temp_array

    print(looper, "time:", datetime.datetime.now() - testtime1)


    df = pd.DataFrame({
                       'a': a_array,
                       'hd':hd_array,
                       })


    #PLOTTING THE CONTOURS AND HISTOGRAMS
    if (j == 0):
        fig, ax, limits, hist_2d_info, contour_2d_info, hist_1d_info = mc.total_plotter(df, 
                                     top_bottom = "bottom", colour = "tab:red", lw = 2, 
                                     cmap = plt.cm.Greys, 
                                     bin_no = 25,
                                     hist_box = "off",
                                     contours = True,
                                     hist2d = False, mode = "off", lw_mm = 2,
                                     colour_c = "RRR", values = "off",
                                     custom_values = custom_values,
                                     value_size = 16,
                                     sigma_type = "actual",
                                     sigma_max = 3,
                                     contour_type = "single",
                                     tick_size = 30
                                     )
        
        base_a = np.median(df["a"])
        base_hd = np.median(df["hd"])
                
    if (j == 1):
    #plottting a 2nd set of data on the same corner plot as the 1st
        fig, ax, limits, hist_2d_info, contour_2d_info, hist_1d_info = mc.total_plotter(df, 
                                     fig, ax, limits, multi_fig = True,
                                     top_bottom = "bottom", colour = "tab:blue", lw = 2, 
                                     cmap = plt.cm.Greys, 
                                     bin_no = 25,
                                     hist_box = "off",
                                     contours = True,
                                     hist2d = False, mode = "off", lw_mm = 2,
                                     colour_c = "BBB", values = "off",
                                     sigma_type = "actual",
                                     sigma_max = 3,
                                     contour_type = "single",
                                     tick_size = 30)

    if (j == 2):
    #plottting a 3rd set of data on the same corner plot as the 1st
        fig, ax, limits, hist_2d_info, contour_2d_info, hist_1d_info = mc.total_plotter(df, 
                                     fig, ax, limits, multi_fig = True,
                                     top_bottom = "bottom", colour = "tab:green", lw = 2, 
                                     cmap = plt.cm.Greys, 
                                     bin_no = 25,
                                     hist_box = "off",
                                     contours = True,
                                     hist2d = False, mode = "off", lw_mm = 2,
                                     colour_c = "GGG", values = "off",
                                     sigma_type = "actual",
                                     sigma_max = 3,
                                     contour_type = "single",
                                     tick_size = 30)

    if (j == 3):
    #plottting a 4th set of data on the same corner plot as the 1st
        fig, ax, limits, hist_2d_info, contour_2d_info, hist_1d_info = mc.total_plotter(df, 
                                    fig, ax, limits, multi_fig = True,
                                    top_bottom = "bottom", colour = "m", lw = 2, 
                                    cmap = plt.cm.Greys, 
                                    bin_no = 25,
                                    hist_box = "off",
                                    contours = True,
                                    hist2d = False, mode = "off", lw_mm = 2,
                                    colour_c = "MMM", values = "off",
                                    sigma_type = "actual",
                                    sigma_max = 3,
                                    contour_type = "single",
                                    tick_size = 14)
    

    ax[1,0].plot(np.median(df["a"]), np.median(df["hd"]), color = "w", marker = "o")

  
########################################
###################KEY##################
########################################

labels = ["GRO2", "GRO3", "GRO4", "GRO5"]
colours = ["tab:red", "tab:blue", "tab:green", "m"]
x_label_pos = 0.5*(limits[0][0] + limits[1][0])
y_label_pos = 0.95*limits[2][0]

x_scale_for_key = -0.4 #arbitary scale for the key,change to move it around
ax[0,0].text(x_scale_for_key*x_label_pos, y_label_pos, "KEY", color = "k", fontsize = 18)

for n in range(len(labels)):
    ax[0,0].text(x_scale_for_key*x_label_pos, y_label_pos - (0.1 * y_label_pos * (n + 1)),
                 labels[n], color = colours[n], 
                 fontsize = 18)

########################################
####PLOTTING SPIN FCOL RELATIONSHIP#####
########################################

ax[1,0].set_xlabel(r"$a$", fontsize = 20, labelpad = 10)
ax[1,1].set_xlabel(r"$f_{col}$", fontsize = 20, labelpad = 10)
ax[1,0].set_ylabel(r"$f_{col}$", fontsize = 20, labelpad = 10)

#determining the spin - colour correction relationship
def r_isco(a):
    
    Z1_A = (1 - a**2)**(1/3)
    Z1_B = (1 + a)**(1/3)
    Z1_C = (1 - a)**(1/3)
    
    Z1 = 1 + Z1_A*(Z1_B + Z1_C)
    
    Z2 = (3*a**2 + Z1**2)**(1/2)
    
    r_isco_A = 3 - Z1
    r_isco_B = 3 + Z1 + 2*Z2
    r_isco_C = (r_isco_A*r_isco_B)**(1/2)
    
    if (a >= 0):
        r_isco_final = 3 + Z2 - r_isco_C
    else:
        r_isco_final = 3 + Z2 + r_isco_C
    
    return r_isco_final

spins = np.linspace(-0.998, 0.998, 1000) #creates an array of 1000 spin values
r_iscos = [] #empty array to fill using ".append"
f_colz = [] #empty array to fill using ".append"
for n in range(len(spins)): #cycle through the array of 1000 spin values
    r_iscos.append(r_isco(spins[n])) #calculate r_isco for each spin and adds
                                     #them to the r_iscos array
    f_colz.append(np.sqrt(r_iscos[n])*0.78) #calculate f_col for each r_isco
                                            #multipled by a constant 0.78
                                            #and adds them to the f_colz array

ax[1,0].plot(spins, f_colz, "w--", linewidth = 3)
ax[0,0].set_xlim(left = -0.998)


########################################
#######TIMERS, PLOTTERS & SAVING########
########################################

save_timer = datetime.datetime.now()
date_now = str(datetime.datetime.now())[0:10]
plt.savefig("GRO_spin_fcol_corner_plot_%s.pdf" %date_now, bbox_inches = "tight")
print("Save-time: ", datetime.datetime.now() - save_timer)

plot_timer = datetime.datetime.now()
plt.show() #use only once everything is plotted
 
print("Plot-time: ", datetime.datetime.now() - plot_timer)
print("Total Runtime: ", datetime.datetime.now() - start_timer)
