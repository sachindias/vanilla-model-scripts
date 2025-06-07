#Courtesy of Sachin Dias O.o
#Do enjoy!! *0*

########################################
#############IMPORT LIBRARIES###########
########################################
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import datetime
import sys
import scipy.stats as stats
#import seaborn as sns #technically not needed unless this option is used

border_width = 2 #set the border width
mpl.rcParams['axes.linewidth'] = border_width #widens all the borders for a nicer look

########################################
##############FIGURE PLOTTER############
########################################

###################################################################################
#Below: All the options for plotting the figures, defaults are also given in 
#function

#df = the Dataframe, this is essential, it can't plot anything without it

#you can ignore fig, ax and limits unless you are plotting a 2nd data set, 
#in which case you MUST include "fig, ax, limits" as outputs of the first 
#use of total_plotter, and "fig, ax, limits, multi_fig = True" as inputs 
#to the 2nd, 3rd, 4th, etc total_plotter used. BEWARE these plots can be 
#hectric to read...

#top_bottom = sets whether the plot is in the top diagonal or bottom diagonal
    #OPTIONS: "top" or "bottom"
    
#hist_box = sets whether the histograms are fully contained in boxes or only on
#the x axis
    #OPTIONS: "on" or "off"
    
#contours = sets whether contours are plotted or not
#hist_2d = sets whether a 2d histogram based upon matplotlib.pyplot.hist2d 
#is plotted or not    
#points = sets whether individual points are plotted or not
    #OPTIONS for all 3: True or False (self explanitory which is which)
    
#mean = sets whether the location of the mean (from 1D hist) is plotted or not
#median = sets whether the location of the median (from 1D hist) is plotted or not
#mode = sets whether the location of the mode (from 1D hist) is plotted or not
    #OPTIONS: "off" >>> don't plot
    #         "on" >>> plot on everything
    #         "on_hist" >>> plot on 1D histograms only
    #         "on_plots" >>> plot on the countour/2d hist/point plots only
    
#lw = sets the linewidth of the 1D histogram
#lw_mm = sets the linewidth of the mean, median or mode lines
#bin_no = number of bins for both the 1D histogram and countour/2d hist/point plots
#colour = sets colour of the points, 1D histograms, or Seaborn plots if incliuded
#alpha = sets alpha (opacity) of the points, 2D histograms or Seaborn plots 
#        if incliuded
#alpha_hist = sets alpha (opacity) of the 1D histograms
#markersize = sets the marker size of the points plots

#cmap = sets the colour map for the 2D histograms, for reasons unknown to me
#       it has a very specific format
    #OPTIONS: plt.cm.<colour scheme> <- a list of colour maps can be found online
    #                                   just take the name of this e.g "viridis" or 
    #                                   "Greys" and use it, see the function itself
    #                                   for an example.

#colour_c = sets the colour for the contour plots.
    #OPTIONS: "RBG" >>> from outside - in, red, blue, green
    #         "transparent" >>> no colouring on the contour plots
    #         PLEASE ADD MORE AS YOU SEE FIT, LOOK FOR "COLOUR CONTOUR OPTIONS"
    
#values = sets whether the median +/- error (the range between the 16th and 84th)
#         quantiles are shown
    #OPTIONS: "on" or "off"
    
#custom_values = put in a list custom names for each value
#value_size = set the size of the text for the values
#label_size = set the size of the text for the axes labels
#tick_size = sets the size of the text for the axes ticks

#error_lines = sets whether the 16th and 84th quantiles are plotted as lines
    #OPTIONS: "on" or "off"
    
#error_colour = sets the colour of the 16th and 84th quantile lines 
#error_linestyle = sets the linestyle of the 16th and 84th quantile lines  

#sigma_type = sets whether contours are plotted using either the actual
#             1, 2 and 3 sigma values in 2D or rather percentages of the data, based
#             on what 1,2 and 3 sigma percenatages are in 1 dimension (i.e ~ 68%,
#             95%, 99%)
    #OPTIONS: "actual" >>> plot actual sigma values in 2D (~ 39%, 86%, 99% of data)
    #         "percentage" >>> plot percentage of data based upon sigma 
    #         percenatages in 1 dimension (~ 68, 95, 99% of data)
#sigma_max = sets whether to use up up 1, 2 or 3 sigma 
    #OPTIONS: 1, 2, 3
#contour_type = sets whether to plot multiple contours or just the maximum sigma
    #OPTIONS: "single" or "multi"
    
#priors = sets whether to plot priors or not and if so what extent to plot
    #OPTIONS: "off" 
    #         "on"  
    #         "on_stretch" >>> will plot priors such that the peak is visible
    #                          NOTICE: this will effect the look of the contour plots

#prior_no = the array with which histograms (from 0 - X) need priors plotter
#prior_type = what types of priors need plotted for each value in the prior_no list

#PLEASE NOTE: THE PRIOR PARAMETERS HAVE BEEN SPECIALISED TO MY DOCUMENTATION
#             AND MY SPECIFIC PRIORS, PLEASE CHECK IF THESE ARE APPLICABLE,
#             IF MORE PRIORS NEED ADDING OR THE CURRENT PRIORS NEED CHANING
#             FOR MORE INFO CHECK scipy.stats DOCUMENTATION       
#h1 = prior parameter 1
#h2 = prior_parameter 2
#h3 = prior parameter 3

#prior_col = sets the colour for priors
#prior_style = sets the line style for the priors

#label = sets the name and displays the label for the set of data
###################################################################################
def total_plotter(df, fig = 1, ax = 1, limits = [[],[]], multi_fig = False, 
               top_bottom = "bottom", hist_box = "on",
               contours = True, hist2d = False, points = False,
               mean = "off", median = "off", mode = "off",
               lw = 1, lw_mm = 1, bin_no = 20, colour = "k",
               alpha = 1, alpha_hist = 1, markersize = 1,
               cmap = plt.cm.viridis, colour_c = "RBG",
               values = "off", custom_values = [None], value_size = 14, label_size = 14,
               tick_size = 12, error_lines = "off", error_colour = "k",
               error_linestyle = "--", sigma_type = "percentage",
               sigma_max = 3, contour_type = "multi",
               priors = "off", prior_no = [], prior_type = [],
               h1 = [], h2 = [], h3 = [], prior_col = "k",
               prior_style = "--", label = " "):
    
    run_timer = datetime.datetime.now()
    
    no_cols = len(df.columns) #number of parameters to plot
    no_points = len(df.iloc[:, 0]) #total number of points
    #no_points = 20000  #restricts the number of points only included
                        #for testing, otherwise leave commented
    
    plot_c = colour #colour for the points, 2D histograms or Seaborn plots
    alpha_p = alpha #alpha for the points, 2D histograms or Seaborn plots
    markersize_p = markersize #marker size for the point plots
    cmap_p = cmap #colour map colours for the 2D histograms
    
    hist_2d_info = []
    contour_2d_info = []
    modes = []
    
    #creates the figure
    if (multi_fig == False):
        print("FIRST DATA SET")
        
        fig, ax = plt.subplots(no_cols, no_cols, figsize = (no_cols*5/1.5,no_cols*5/1.5))  
        
        #This is the space between plots, I personally like a bit of space,
        #but if you do not, then reduce these to 0, or increase if you want more space,
        #however note this will decrease each plot's size, so you may need to incease
        #the figure size above in line with this. The "top" plots have more space
        #just to fit the values on, if they are selected to be on.
        space = 0.05
        plt.subplots_adjust(wspace=space, hspace=space)
        
        if (top_bottom == "top" and values == "on"):
            plt.subplots_adjust(wspace=space*5, hspace=space*5)
        
        min_lim = []
        max_lim = []
        max_y = []
        
    #needed as True if multiple sets of data are plotted on the same graph
    elif (multi_fig == True):
        print("NEXT DATA SET")
        
        min_lim = limits[0]
        max_lim = limits[1]
        max_y = limits[2]

    #plotting the histograms
    for n in range(no_cols):
        col = df.iloc[:, n][0:no_points] #extracts data from Dataframe
        
        hist_1d_info = [[],[]]
        
        height, bins1d, patches = ax[n,n].hist(col, bins = bin_no, 
                                density = True, color = "gray", 
                                edgecolor = plot_c, linewidth = lw,
                                histtype = "step", alpha = alpha_hist)
        
        hist_1d_info[0].append(height)
        hist_1d_info[1].append(bins1d)
     
        #calculates the value of the median and range between 16th and 84 quantile
        q50 = np.percentile(col, 50)
        q16 = np.percentile(col, 16)
        q84 = np.percentile(col, 84)
            
        dx_down = q50 - q16
        dx_up = q84 - q50
        
        #plots the value of the median and range between 16th and 84 quantile
        if (values == "on"):
            if (custom_values == [None]):
                ax[n,n].set_title(r"%s = %.2f$^{+%.2f}_{-%.2f}$" %(df.columns[n], q50, dx_up, dx_down),
                              size = value_size)
            else:
                ax[n,n].set_title(r"%s = %.2f$^{+%.2f}_{-%.2f}$" %(custom_values[n], q50, dx_up, dx_down),
                              size = value_size)
        
        #plots the 16th and 84th quantiles as vertical lines
        if (error_lines == "on"):
            ax[n,n].axvline(q16, color = error_colour, linestyle = error_linestyle, linewidth = lw_mm)
            ax[n,n].axvline(q84, color = error_colour, linestyle = error_linestyle, linewidth = lw_mm)
        
        #saves the min and max bin values
        if (multi_fig == False):
            min_lim.append(np.min(bins1d))
            max_lim.append(np.max(bins1d))
            max_y.append(np.max(height))
        
        #sets limits based upon max and min for both data sets
        elif (multi_fig == True):
            if (np.min(bins1d) < min_lim[n]):
                min_lim[n] = np.min(bins1d)
            if ((np.max(bins1d) > max_lim[n])):
                max_lim[n] = np.max(bins1d)
            if ((np.max(height) > max_y[n])):
                max_y[n] = np.max(height)
        
        #removes those pesky y values, who needs them XD
        ax[n,n].set_yticks([])
        
        #plotting mean
        if (mean == "on" or mean == "on_hist"):
            ax[n,n].axvline(np.mean(col), color = "tab:red", linewidth = lw_mm)
        
        #plotting median, as with mean, feel free to change colours
        if (median == "on" or median == "on_hist"):
            ax[n,n].axvline(np.median(col), color = "tab:green", linewidth = lw_mm)

        #calculates the mode as the middle of the highest bin in 1-dimension
        #saves into the list "modes"
        mode_pos = np.where(height == np.max(height))[0][0]
        mode_no = (bins1d[mode_pos] + bins1d[mode_pos + 1])/2
        modes.append(mode_no)
        
        #plotting mode, again change colour if you want
        if (mode == "on" or mode == "on_hist"):
            ax[n,n].axvline(mode_no , color = "tab:blue", linewidth = lw_mm)
        
        #This sets which x labels to include, depending on if you want the contour 
        #plots in the top diagonal or bottom diagonal, the labels are determined by 
        #the name given in the Dataframe, so name them appropriately :)
        if (top_bottom == "bottom"):
            if (n != no_cols - 1):
                ax[n,n].set_xticks([])
            else:
                ax[n,n].set_xlabel(df.columns[no_cols - 1], size = label_size)
                ax[n,n].tick_params(axis='x', which='major', labelsize = tick_size, width = border_width)
    
        elif (top_bottom == "top"):
            ax[n,n].set_xlabel(df.columns[n], size = label_size)
            ax[n,n].tick_params(axis='x', which='major', labelsize = tick_size, width = border_width)
        
        #sets the limits of the histograms, to the same as what the contours will
        #be set to
        ax[n,n].set_xlim(min_lim[n], max_lim[n])
    
        #this removes the box around the histograms, personal preference whether,
        #you prefer this or not. 
        if (hist_box == "off"):
           ax[n,n].spines['top'].set_visible(False)
           ax[n,n].spines['left'].set_visible(False)
           ax[n,n].spines['right'].set_visible(False)
           
        #plots the specified priors
        if (priors == "on" or priors == "on_stretch"):
            if (n in prior_no):
                pn = np.where(np.array(prior_no) == n)[0][0] #pos in prior_no array
                
                #x and y values, straight forward enough yeh?
                x_prior = np.linspace(0.01*min_lim[n], 5*max_lim[n], 300)
                y_prior = prior_plotter(x_prior, pn, prior_type, h1, h2, h3)
           
                ax[n,n].plot(x_prior, y_prior, color = prior_col, linestyle = prior_style,
                             linewidth = lw)
                
                #if "on_stretch" is selected, will make sure the peak of the 
                #prior is visible
                if (priors == "on_stretch"):
                    max_y_pos = np.where(y_prior == np.max(y_prior))[0][0]
                    if (1.05 * x_prior[max_y_pos] > max_lim[n]):
                        ax[n,n].set_xlim(right = 1.05*x_prior[max_y_pos])
                        max_lim[n] = 1.05 * x_prior[max_y_pos]
                        
                    if (0.95 * x_prior[max_y_pos] < min_lim[n]):
                        ax[n,n].set_xlim(left = 0.95*x_prior[max_y_pos])
                        max_lim[n] = 0.95 * x_prior[max_y_pos]          
    
    #progress bar (this calculates the total no. of contour plots)
    total_prog = 0
    for pc in range(no_cols):
        total_prog = total_prog + pc
        
    prog_counter = 1  
       
    #plots contours in the bottom diagonal       
    if (top_bottom == "bottom"):    
        for o in range(no_cols - 1):
            for n in range(o + 1, no_cols):           
                hist_info, contour_info = contour_plots(df, contours, hist2d, points, o, n, ax, no_points, no_cols,
                              bin_no, max_lim, min_lim, median, mean, mode, lw_mm, modes,
                              plot_c, alpha_p, markersize_p, cmap_p, colour_c, label_size,
                              tick_size, sigma_type, sigma_max, contour_type)                    
 
                hist_2d_info.append(hist_info)
                contour_2d_info.append(contour_info)
                
                #progress bar
                sys.stdout.write('\r')
                if (prog_counter < total_prog):
                    sys.stdout.write("[%-40s] %d%%" 
                                     %('='*(int(np.round((40/total_prog)*prog_counter))), 
                                     (100/total_prog)*prog_counter)) 
                else:
                    sys.stdout.write("[%-20s] %d%%\n" %('='*40, 
                                     100)) 
                sys.stdout.flush()           
                
                prog_counter = prog_counter + 1
                
        for o in range(no_cols):
            for n in range(0, o):   
                ax[n,o].set_visible(False)
    
    #plots contours in the top diagonal             
    if (top_bottom == "top"):
        for o in range(no_cols):
            for n in range(0, o):           
                hist_info, contour_info = contour_plots(df, contours, hist2d, points, o, n, ax, no_points, no_cols,
                              bin_no, max_lim, min_lim, median, mean, mode, lw_mm, modes, 
                              plot_c, alpha_p, markersize_p, cmap_p, colour_c, label_size,
                              tick_size, sigma_type)                  
    
                hist_2d_info.append(hist_info)
                contour_2d_info.append(contour_info)
                
                #progress bar
                sys.stdout.write('\r')
                if (prog_counter < total_prog):
                    sys.stdout.write("[%-40s] %d%%" 
                                     %('='*(int(np.round((40/total_prog)*prog_counter))), 
                                     (100/total_prog)*prog_counter)) 
                else:
                    sys.stdout.write("[%-20s] %d%%\n" %('='*40, 
                                     100)) 
                sys.stdout.flush()           
                
                prog_counter = prog_counter + 1
  
        for o in range(no_cols -1):
            for n in range(o + 1, no_cols):   
                ax[n,o].set_visible(False)

    limits = [min_lim, max_lim, max_y]
    print("\nRuntime: ", datetime.datetime.now() - run_timer, "\n")  

    return fig, ax, limits, hist_2d_info, contour_2d_info, hist_1d_info

########################################
#######CONTOUR/HIST/POINT PLOTTER#######
########################################

#THIS PLOTS THE CONTOURS, IGNORE INPUTS, IMPORTANT ONES HAVE ALREADY BEEN EXPLAINED
def contour_plots(df, contours, hist2d, points, o, n, ax, no_points, no_cols, bin_no,
                  max_lim, min_lim, median, mean, mode, lw_mm, modes, plot_c, alpha_p,
                  markersize_p, cmap_p, colour_c, label_size, tick_size, sigma_type,
                  sigma_max, contour_type):
    
    mulitple_counter = 0 #a counter to deal with duplicates
    
    if (contours == True):
        
        ########################################
        #########SEABORN COUNTOUR PLOTS#########
        ########################################
        #This will use Seaborn to plot contours, they look far nicer, but 
        #idk what the contour values are? maybe 10th percentiles? either way,
        #the option to use them is here, but they take a very long time if
        #the data set is large. ALSO if using, please comment out: MY OWN
        #CONTOUR PLOTS or things may look weird
        
        #sns.kdeplot(df.iloc[:, o][0:no_points], df.iloc[:, n][0:no_points],
        #           color = plot_c, ax = ax[n,o], alpha = alpha_p)
        
        
        ########################################
        ##########MY OWN COUNTOUR PLOTS#########
        ########################################
        
        #plots an invisible histogram, basically just to get the data,
        #I'm sure there is a better way out there, maybe np.hist???
        hist_info = ax[n,o].hist2d(df.iloc[:, o][0:no_points], 
                                      df.iloc[:, n][0:no_points], 
                                      bins = bin_no, 
                                      alpha = 0)
        
        #widths of bins in each direction (bin widths are more or less equal)
        x_dist = hist_info[1][1] - hist_info[1][0]
        y_dist = hist_info[2][1] - hist_info[2][0]
        
        #max and min x-axis and y-axis values
        min_x = hist_info[1][0]
        min_y = hist_info[2][0]
        max_x = hist_info[1][bin_no]
        max_y = hist_info[2][bin_no]
        
        bin_area = (x_dist * y_dist)
        total_mass_density = 0
        contour_info = [[],[],[], []]
        
        #calculates the total mass density across all bins
        for p in range(bin_no):
            for q in range(bin_no):
                mass_density = hist_info[0][p][q] * bin_area
                total_mass_density = total_mass_density + mass_density
    
        #calculates the mass density for each bin, p for x-axis, q for y-axis,
        #and saves these densities and their positions
        for p in range(bin_no):
            for q in range(bin_no):
                per_mass_density = (hist_info[0][p][q] * bin_area)/total_mass_density
                contour_info[0].append(per_mass_density)
                contour_info[1].append(p)
                contour_info[2].append(q)

        #orders a column of decending mass densities, UNRELATED to p and q positions
        #in contour_info
        sorted_array = np.sort(contour_info[0])
        contour_info[3] = sorted_array[::-1]
        
        if (sigma_type == "actual"):
            sigma_array = [0.3935, 0.8647, 0.9889] # 1, 2 and 3 sigmas worth of data, with
                                                   # each sigma being the value of sigma in      
                                                   # 2 dimensions based upon the generalised
                                                   # regularised incomplete gamma function
                                                   # type Q(n/2,0,m^2/2) into Wolfram|Alpha
                                                   # where n is no. dimensions, m is no sigma
                                                       
        elif (sigma_type == "percentage"):
            sigma_array = [0.6827, 0.9545, 0.9973] # 1, 2 and 3 "sigmas" worth of data,
                                                   # Not actually sigmas, but 1 ~ 68.27% 
                                                   # of the data, 2 ~ 95.45% of the data,
                                                   # and 3 ~ 99.73% of the data, to see 
                                                   # the exact percentages uncomment
                                                   # print statement below (see ***)
                                                    
        ########################################
        #########COLOUR CONTOUR OPTIONS#########
        ########################################
        #please add more if you want
        #also I am aware RBG may be confusing with RGB...
        yellow_colour = "#C1B732"
        green_colour = "#67D8B9"
        pink_colour = "#DC267F"
        blue_colour = "#648FFF"
        
        if (colour_c == "RBG"):
            colours = ["r", "b", "g"]
            alpha_cont = 0.5
        if (colour_c == "RRR"):
            colours = ["r", "r", "r"]
            alpha_cont = 0.4
        if (colour_c == "BBB"):
            colours = ["b", "b", "b"]
            alpha_cont = 0.4
        if (colour_c == "GGG"):
            colours = ["g", "g", "g"]
            alpha_cont = 0.4
        if (colour_c == "MMM"):
            colours = ["m", "m", "m"]
            alpha_cont = 0.4
        if (colour_c == "YYY"):
            colours = [yellow_colour,yellow_colour,yellow_colour]
            alpha_cont = 0.4
        if (colour_c == "G2G2G2"):
            colours = [green_colour,green_colour,green_colour]
            alpha_cont = 0.4
        if (colour_c == "PPP"):
            colours = [pink_colour,pink_colour,pink_colour]
            alpha_cont = 0.4
        if (colour_c == "B2B2B2"):
            colours = [blue_colour,blue_colour,blue_colour]
            alpha_cont = 0.4
            

        
 
        #for each of the "sigmas", working backwards, so 3 sigma, 2 sigma, 1 sigma,
        #starts at the bin with largest mass density and works down until the sigma
        #limit is reached. At this point it saves their positions into the positional
        #array, which is either 0 for not included, and 1 for include (i.e in sigma
        #range). If you want to look at the positional array, NOTE the x-axis increases
        #vertically downwards and the y-axis increases horizontally rightwards
        
        if (contour_type == "multi"):
            sigma_array = sigma_array[:sigma_max]
        elif (contour_type == "single"):
            sigma_array = [sigma_array[sigma_max-1]]  
        
        for sigma_count in range(len(sigma_array)):
            positional_array = np.zeros((bin_no, bin_no))
            sigma_total = 0 #the total actual sigma value for counting
            loc_counter = 0 #the number of bins from 0 being the highest mass density
            sigma_aim = sigma_array[-(sigma_count+1)]
            while (sigma_total < sigma_aim):
                sigma_total = sigma_total + contour_info[3][loc_counter]
                
                #this deals with places where the percentage mass density is the same for
                #multiple locations. It does this by sequentially moving through the positions
                #in an array if there are multiple or accepting the 0th position if there are
                #none
                location_array = np.where(contour_info[0] == contour_info[3][loc_counter])[0]
                if (len(location_array) > 1):
                    location = location_array[mulitple_counter]
                    mulitple_counter = mulitple_counter + 1
                    
                    if (mulitple_counter == len(location_array)):
                        mulitple_counter = 0  
                            
                else:
                    location = location_array[mulitple_counter]
                    
                location_x = contour_info[1][location]
                location_y = contour_info[2][location]
                
                positional_array[location_x][location_y] = 1
                
                loc_counter = loc_counter + 1
            
            mulitple_counter = 0 #in case it hasn't updated due to hitting sigma limit
            
            #***
            #print("Attempted sigma = %s, Actual sigma = %.4f" %(sigma_aim, sigma_total)) 
            #***

            #cycles through all bin positions which need to be plotted
            #i.e where the positional_array = 1
            for loc_x_count in range(bin_no):
                for loc_y_count in range(bin_no):
                    if (positional_array[loc_x_count][loc_y_count] == 1):    
                        
                        #Plots the filled in contour plots, if the colour (colour_c) is 
                        #set to "transparent", then no filled in contour plots will be
                        #plotted, just their outlines (see below)
                        #for more colour options, please add them above in
                        #COLOUR CONTOUR OPTIONS, default is "RBG" = red, blue, green
                        if (colour_c != "transparent"):
                            ax[n,o].fill_between([min_x + loc_x_count*x_dist,
                                                  min_x + (loc_x_count + 1)*x_dist],
                                                 y1 = min_y + loc_y_count*y_dist,
                                                 y2 = min_y + (loc_y_count + 1)*y_dist,
                                                 color = colours[sigma_count],
                                                 alpha = alpha_cont)
                        
                        
                        #Plots the outline of the contour plots, based
                        #on if there are other points next to it. The
                        #colour is set to a default of black ("k")
                        if (loc_y_count == 0):
                            ax[n,o].plot([min_x + loc_x_count*x_dist, 
                                          min_x + (loc_x_count + 1)*x_dist], 
                                         [min_y,min_y],
                                         color = "k")
                            
                        else:
                            if (positional_array[loc_x_count][loc_y_count - 1] == 0):
                                ax[n,o].plot([min_x + loc_x_count*x_dist, 
                                              min_x + (loc_x_count + 1)*x_dist], 
                                             [min_y + loc_y_count*y_dist,
                                              min_y + loc_y_count*y_dist],
                                             color = "k")           
                        
                        if (loc_y_count == bin_no - 1):
                            ax[n,o].plot([min_x + loc_x_count*x_dist, 
                                          min_x + (loc_x_count + 1)*x_dist], 
                                         [max_y, max_y],
                                         color = "k")
                            
                        else:
                            if (positional_array[loc_x_count][loc_y_count + 1] == 0):
                                ax[n,o].plot([min_x + loc_x_count*x_dist, 
                                              min_x + (loc_x_count + 1)*x_dist],
                                              [min_y + (loc_y_count + 1)*y_dist,
                                              min_y + (loc_y_count + 1)*y_dist],
                                              color = "k")
                            
                        if (loc_x_count == 0):
                            ax[n,o].plot([min_x,min_x], 
                                         [min_y + loc_y_count*y_dist,
                                          min_y + (loc_y_count + 1)*y_dist],
                                         color = "k")
                                        
                        else:       
                            if (positional_array[loc_x_count - 1][loc_y_count] == 0):
                                ax[n,o].plot([min_x + loc_x_count*x_dist,
                                              min_x + loc_x_count*x_dist],
                                             [min_y + (loc_y_count)*y_dist,
                                              min_y + (loc_y_count + 1)*y_dist],
                                             color = "k")
 
                        if (loc_x_count == bin_no - 1):
                            ax[n,o].plot([max_x, max_x], 
                                         [min_y + loc_y_count*y_dist,
                                          min_y + (loc_y_count + 1)*y_dist],
                                         color = "k")

                        else:                                 
                            if (positional_array[loc_x_count + 1][loc_y_count] == 0):
                                ax[n,o].plot([min_x + (loc_x_count + 1)*x_dist,
                                              min_x + (loc_x_count + 1)*x_dist],
                                             [min_y + (loc_y_count)*y_dist,
                                              min_y + (loc_y_count + 1)*y_dist],
                                             color = "k")
     
    #This will actually plot the 2d histograms based on the function hist2d   
    if (hist2d == True):
        hist_info = ax[n,o].hist2d(df.iloc[:, o][0:no_points], 
                                   df.iloc[:, n][0:no_points], 
                                   bins = bin_no, cmap = cmap_p, 
                                   alpha = alpha_p)  
        if (contours == False):
            contour_info = None
    
    #This will plot every single point in the data, beware for a lot of points,
    #this will make the histograms or the contour plots non-visible
    if (points == True):
        print("plotting")
        ax[n,o].plot(df.iloc[:, o][0:no_points], df.iloc[:, n][0:no_points], 
                     linestyle = "", marker = "o", color = plot_c, 
                     alpha = alpha_p, markersize = markersize_p)
        
        if (contours == False):
            contour_info = None
            if (hist2d == False):
                hist_info = None
    
    #if you ask it to plot nothing... why?
    if (contours == False and hist2d == False and points == False):
        contour_info = None       
        hist_info = None
                
    #sets contour plot limits to that of the histograms
    ax[n,o].set_xlim(min_lim[o], max_lim[o])
    ax[n,o].set_ylim(min_lim[n], max_lim[n])
    
    #removes x- & y-axis ticks and labels unless they're the bottom/left most
    if (o != 0):
        ax[n,o].set_yticks([])
        ax[n,o].set_ylabel("")
    else:
        ax[n,o].set_ylabel(df.columns[n], size = label_size)
        ax[n,o].tick_params(axis='y', which='major', labelsize = tick_size, width = border_width)
    if (n != no_cols - 1):
        ax[n,o].set_xticks([])
        ax[n,o].set_xlabel("")
    else:
        ax[n,o].set_xlabel(df.columns[o], size = label_size)
        ax[n,o].tick_params(axis='x', which='major', labelsize = tick_size, width = border_width)
       
    #plots vertical median, mode and mean lines, BASED UPON THE HISTOGRAMS, not the
    #contour plots, as these are not necessarily the same.
    #Colours are chosen based on what I thought was good, but feel free to change
    if (median == "on" or median == "on_plots"):
        ax[n,o].axhline(np.median(df.iloc[:, n][0:no_points]), color = "tab:green", linewidth = lw_mm)
        ax[n,o].axvline(np.median(df.iloc[:, o][0:no_points]), color = "tab:green", linewidth = lw_mm)
    if (mean == "on" or mean == "on_plots"):
        ax[n,o].axhline(np.mean(df.iloc[:, n][0:no_points]), color = "tab:red", linewidth = lw_mm)
        ax[n,o].axvline(np.mean(df.iloc[:, o][0:no_points]), color = "tab:red", linewidth = lw_mm)
    if (mode == "on" or mode == "on_plots"):
        ax[n,o].axhline(modes[n], color = "tab:blue", linewidth = lw_mm)
        ax[n,o].axvline(modes[o], color = "tab:blue", linewidth = lw_mm)

    return hist_info, contour_info

########################################
##############PRIOR PLOTTER#############
########################################

#this generates the y_axis for the priors for the posteriors specified
def prior_plotter(x_prior, pn, prior_type, h1, h2, h3):
    
    #Gaussian prior
    if (prior_type[pn] == "gaussian"):
        y_prior = stats.norm.pdf(x_prior, h1[pn], h2[pn])
    
    #gamma prior
    elif (prior_type[pn] == "gamma"):
        y_prior = stats.gamma.pdf(x_prior, h1[pn], scale = 1/h2[pn])  
        
    #shifted log normal prior
    elif (prior_type[pn] == "shiftedlognormal"):
        h2 = np.exp(h2[pn]) #has to convert from ln to exp for the scipy function
        y_prior = stats.lognorm.pdf(x_prior, h3[pn], h1[pn], h2) 
        
    return y_prior
    
#################################################################################
#When running the actual contours, none of the outputs are required unless you
#want to plot multiple figures on the same plots or want to extract some 
#information for checking purposes

#fig, ax = the plots themselves
#limits = [[min histogram values],[maximum histogram values]]
#hist_2d_info = some information about the contour plots themselves
    #>>>these will be set to "None" if just points are plotted
    #>>>otherwise will provide the heights and bin information for every plot
    #>>>Format: hist_2d_info[which plot][height = 0, x_axis = 1, y_axis = 2]
#contour_2d_info = some more information for plotting the contours
    #>>>these will be set to "None" if just points or just histograms are plotted
    #>>>otherwise will provide the mass density and bin positions for every plot
    #>>>Format: contour_2d_info[which plot][mass density = 0
                                           #x bin position = 1,
                                           #y bin position = 2, 
                                           #decending mass density = 3]
                #>>> NOTE decending mass densities are not at all correlated 
                #    with the x and y bin positions!!
                
#NOTE: Using a Pandas dataframe is not strictly necessary, 
    #but as I built this for using those data frames originally, 
    #if you don't want to use them you will need to change the 
    #functions yourself! 
#################################################################################