import numpy  as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import scipy as sci

XRBs = [
        "GS 1354-64",
        "GRO J1655-40", 
        "V4641 Sgr", 
        "MAXI J1820+070", 
        "GX 339-4"
        ]

#PARALLAX ESTIMATES & ERRORS (UNITS: mas)
Parx = [1.2762, 0.3284, 0.2135, 0.348, 0.2086]
err_Parx = [0.534, 0.048, 0.026, 0.033, 0.120]

#TO CHECK ONE DIST
Position = 4 #OUT OF 5
Pos = Position - 1
XRBs = [XRBs[Pos]]
Parx = [Parx[Pos]]
err_Parx = [err_Parx[Pos]]

#GENERATE POSTERIOR
for n in range(len(XRBs)):
    print("%s:" %XRBs[n], "\n")
    parx = Parx[n]
    err_parx = err_Parx[n]
    
    # distance r, in steps of dr (units: kpc)
    dr = 0.01
    r = np.arange(0, 20, dr)
    
    # likelihood (normal for parx data; as a function of 1/r)
    r_like = stat.norm.pdf(parx, 1/r, err_parx)
    
    # prior on distance (exponential from Gandhi et al. 2019; MNRAS)
    L = 2.17
    r_prior = r**2 * np.exp(-r/L)
    
    # normalise the prior
    A = sum(r_prior) * dr
    r_prior = r_prior / A
    
    # posterior = likelihood * prior / constant
    r_posterior = r_prior * r_like
    
    # normalise posterior
    A = sum(r_posterior) * dr
    r_posterior = r_posterior / A

    # normalise the likelihood
    A = sum(r_like) * dr
    r_like = r_like / A
    
    sum_e = 0
    mean = 0
    dist = []
    prob = []
    dx = []
    
    r_posterior_new = r_posterior[1:]
    r_new = r[1:]
    
    for count in range (len(r_posterior_new)):
    
        dist_temp = r_new[count]
    
        if (count != 0):
            dx_temp = r_new[count] - r_new[count - 1]
            if (count < 7):
                dx.append(dx_temp)    
            dx.append(dx_temp)
        elif (count == 0):
            dx_temp = r_new[0]
            if (count < 7):
                dx.append(r_new[0])
            
        sum_e = sum_e + (r_posterior_new[count] * dx_temp)
        mean = mean + (dist_temp * r_posterior[count] * dx_temp)
    
    mode_pos = np.where(r_posterior == np.max(r_posterior))[0][0]
    mode = r[mode_pos]
    

    #CALCULATE AND PRINT KEY VALUES
    alpha = mean/(mean - mode)
    theta = mode/(alpha - 1)
    
    print("alpha: %.4f" %alpha)
    print("theta: %.4f" %(theta))
    print("beta: %.4f" %(1/theta))
    print("Modal Distance:", r[np.where(r_posterior == max(r_posterior))[0][0]], "kpc\n\n")
    
    #GENERATE GAMMA PRIOR
    x_min = 0
    x_max = 20
    x = np.linspace(x_min, x_max, 1000)
    y = sci.stats.gamma.pdf(x, alpha, scale = theta)
    
    #PLOTTING COMMANDS
    plt.figure(figsize = (10 ,7))
    plt.xlim(0, 20)
    plt.xlabel("Distance (kpc)", fontsize = 30)
    plt.ylabel("posterior p(r|data)", fontsize = 30)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    
    plt.plot(r, r_posterior, "k", linewidth = 3)
    plt.plot(x, y, color = "tab:red", linewidth = 3)
    plt.axvline(x = r[np.where(r_posterior == max(r_posterior))[0][0]], color = 'tab:green', linewidth = 2, linestyle = "--")
    plt.show()