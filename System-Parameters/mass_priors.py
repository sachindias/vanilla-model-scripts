import matplotlib.pyplot as plt
import scipy as sci
import scipy.stats as st
import numpy as np
import datetime
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

#TIMER
start = datetime.datetime.now()

mpl.rcParams['axes.linewidth'] = 2.5

#-----------------------------------------
#FUNCTIONS

#CALCULATES THE ISOTROPIC DISTRIBUTION ON INCLINATION
def iso_prior(n, i_start = 0, i_end = 90):
    I1 = np.pi*i_start/180
    I2 = np.pi*i_end/180
    
    low = np.cos(I2)
    upp = np.cos(I1)

    x = np.random.uniform(low, upp, n)
    y = (np.arccos(x) * 180)/np.pi
    
    return (y)

#CALCULATES MINIMUM INCLINATION IN RADIANS
def min_incl(f_M, q, M_max):
    incl = np.arcsin((((1+q)**2)*f_M/M_max)**(1/3))
    incl = incl * 180/np.pi
    return incl

#CALCULATES MAXIMUM INCLINATION IN RADIANS
def max_incl(q):
    incl = np.arccos(0.462*(q/(1+q))**(1/3))
    incl = incl * 180/np.pi
    return incl

#CALCULATES THE MASS DISTRIBUTION
#N = NORMAL DISTRIBUTION
#U = UNIFORM DISTRIBUTION
#I = ISOTROPIC DISTRIBUTION
def M_dist_func(n, no_points, f_1, f_2, q_type, q_1, q_2, i_type, i_1, i_2):
    f_dist = np.random.normal(f_1[n], f_2[n], no_points)
    
    if (q_type[n] == "N"):
        q_dist = np.random.normal(q_1[n], q_2[n], no_points)
    if (q_type[n] == "U"):
        q_dist = np.random.uniform(q_1[n], q_2[n], no_points)

    if (i_type[n] == "N"): 
                
        i_1_rad = np.pi*i_1[n]/180
        i_2_rad = np.pi*i_2[n]/180
        i_dist = np.random.normal(i_1_rad, i_2_rad, no_points)
        
    if (i_type[n] == "I"):    
        i_dist = iso_prior(no_points, i_1[n], i_2[n])
        i_dist = np.pi*i_dist/180
    
    M_dist = (f_dist * (1 + q_dist)**2)/(np.sin(i_dist))**3
    M_dist = M_dist[ (M_dist <= 50) & (M_dist >= 2)] #KEEP MASS DISTRIBUTION WITHIN THE RANGE
    
    return M_dist, i_dist, q_dist, f_dist

#-----------------------------------------

M_total = []

M_max = 50.0
q = 0

XRB = [
       "GS 1354-64",
       "GRO J1655-40", 
       "GX 339-4",
       "V4641 Sgr", 
       "MAXI J1820+070"
       ]

f_1 = [5.73, 2.73, 1.91, 2.74, 5.18]
f_2 = [0.29, 0.09, 0.08, 0.04, 0.15]

q_1 = [0.12, 0.21, 0.0, 0.45, 0.072]
q_2 = [0.04, 0.45, 0.23, 0.04, 0.012]
q_type = ["N", "U", "U", "N", "N"] #WHICH TYPE OF DISTRIBUTION

i_1 = [29, 69, 19, 72.3, 63]
i_2 = [80, 3, 80, 4.1, 3]
i_type = ["I", "N", "I", "N", "N"] #WHICH TYPE OF DISTRIBUTION

#NUMBER OF ESTIMATES
no_points = 20000

#TO CHECK ONE DIST
Position = 3 #Out of 5 (INCLUDING 2 FOR GX 339-4)
Pos = Position - 1
XRB = [XRB[Pos]]
f_1 = [f_1[Pos]]
f_2 = [f_2[Pos]]

q_1 = [q_1[Pos]]
q_2 = [q_2[Pos]]
q_type = [q_type[Pos]]

i_1 = [i_1[Pos]]
i_2 = [i_2[Pos]]
i_type = [i_type[Pos]]


#CALCULATE MASS POSTERIOR
for n in range(len(XRB)):
    
    #PRINT MAX AND MIN INCLINATION VALUES BASED ON MASS RATIO
    print("Working on: %s" %XRB[n])
    print("min inclination:   %.2f deg" %min_incl(f_1[n], q, M_max))
    if (q_type[n] == "U"):
        print("max inclination:   %.2f deg" %max_incl(q_1[n]))
    if (q_type[n] == "N"):
        print("max inclination:   %.2f deg" %max_incl(q_1[n] - q_2[n]))        
    
    M_dist, i_dist, q_dist, f_dist = M_dist_func(n, no_points, f_1, f_2, q_type, q_1, q_2, i_type, i_1, i_2)
    
    #TOTAL   
    M_total = np.concatenate([M_total, M_dist])

    #ADD MORE TILL 20K POINTS (AS VALUES OUTSIDE 2 - 50 SOLAR MASSES ARE IGNORED)
    print("\n")
    while (len(M_total) < no_points):
        print("TOO SHORT: %s" %len(M_total))
        
        #ADDS 1000 POINTS PER TIME TILL THE POINTS ARE 20K
        M_dist_N, i_dist_N, q_dist_N, f_dist_N = M_dist_func(n, 1000, f_1, f_2, q_type, q_1, q_2, i_type, i_1, i_2)
        M_total = np.concatenate([M_total, M_dist_N])       

    print("REMOVING EXCESS %s POINTS" %(len(M_total) - no_points))
    M_total = M_total[0:no_points]
    print("TOTAL POINTS FINE: %s \n" %len(M_total))
    
    #HOW TO DEAL WITH THE SPECIAL CASE OF GX 339-4 WHERE WE USE 2 MASS FUNCTIONS
    if (XRB[n] == "GX 339-4"):
        M_total_2 = []
        M_dist_2, i_dist_2, q_dist_2, f_dist_2 = M_dist_func(n, no_points, [5.8], [0.5], q_type, q_1, q_2, i_type, i_1, i_2)

        M_total_2 = np.concatenate([M_total_2, M_dist_2])

        #ADDS 1000 POINTS PER TIME TILL THE POINTS ARE 20K
        print("\n")
        if (len(M_total_2) < no_points):
            while (len(M_total_2) < no_points):
                print("TOO SHORT: %s" %len(M_total_2))
                
                M_dist_N, i_dist_N, q_dist_N, f_dist_N = M_dist_func(n, 1000, [5.8], [0.5], q_type, q_1, q_2, i_type, i_1, i_2)
                M_total_2 = np.concatenate([M_total_2, M_dist_N])       
        if (len(M_total_2) >= no_points):
            print("REMOVING EXCESS %s POINTS" %(len(M_total_2) - no_points))
            M_total_2 = M_total_2[0:no_points]
            print("\nTOTAL POINTS FINE: %s \n" %len(M_total_2))

        M_total = np.concatenate([M_total, M_total_2], axis=None)
    
    #-------------------------------
    #ESTIMATE LOG NORMAL DISTRIBUTIONS
    
    dist = getattr(st, "lognorm")
    param_ln = dist.fit(M_total)
    
    arg = param_ln[:-2][0]
    loc = param_ln[-2]
    scale = param_ln[-1]
    
    #PRINT SHIFTED LOGNORMAL VARIABLES FOR PYTHON [AND XSPEC]
    print("VARIABLES")
    print("arg [h3]: %.4f" %arg)
    print("loc [h1]: %.4f" %loc)
    print("scale [h2]: %.4f || mu [ln(scale)]: %.4f\n" %(scale, (np.log(scale))))
    
    #PRINT PRIOR
    print("PRIOR")
    print("bayes 39 shiftedlognormal %.4f %.4f %.4f\n" %(loc, scale, arg))
   
    #MIN AND MAX MASSES IN DISTRIBUTION
    print("MIN MASS: %.4f" %(np.min(M_dist)))
    print("MAX MASS: %.4f" %(np.max(M_dist)))
   
    #GET VALUES TO PLOT
    x_fit_ln = np.linspace(0, 50, 1000)
    y_fit_ln = st.lognorm.pdf(x_fit_ln, arg, loc, scale)
    
    #-------------------------------
    
    #PLOTTING
    plt.figure(figsize = (10 ,7))
    
    #HISTOGRAM OF POSTERIOR
    y, x, _ = plt.hist(M_total, bins = 50, linewidth = 3, color = "k", alpha = 1, histtype = "step", density = True)
    
    #PRIOR LINE
    plt.plot(x_fit_ln, y_fit_ln, color = "tab:red", linewidth = 3)
    
    plt.xlabel(r'Mass$_{BH}$ ($M_{\odot}$)', fontsize = 30)
    plt.ylabel('posterior p(M|data)', fontsize = 30)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlim(0,50)