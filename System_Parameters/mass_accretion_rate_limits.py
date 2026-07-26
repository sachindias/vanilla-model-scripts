import numpy as np
import datetime

#TIMER
start = datetime.datetime.now()

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

#CALCULATES THE EDDINGTON LUMINOSITY (IN WATTS)
def L_Edd(M):
    G = 6.67408e-11         #Gravitational Constant ( m^3 Kg^-1 s^-2 )
    m_p = 1.67262e-27       #Mass of a Proton ( Kg )
    c = 299792458           #Speed of Light ( m s^-1 )
    sigma_T = 6.65245e-29   #Thomson Scattering Cross Section for Electron ( m^2 )
    M_sun = 1.989e30        #Mass of the Sun ( Kg )
    
    L_Edd = (4*np.pi*G*M*M_sun*m_p*c)/sigma_T #WATTS
    
    return L_Edd

#CALCULATES THE EDDINGTON MASS ACCRETION RATE (IN 10^18 GRAMMS PER SECOND)
#FOR NON-SPINNING BLACK HOLES
def Mdd_Edd(L_Edd):
    eta = 0.057
    c = 299792458
    
    Mdd_Edd = L_Edd/(eta*c**2)  #Kg s^-1
    Mdd_Edd = Mdd_Edd * 1000    #g s^-1
    Mdd_Edd = Mdd_Edd * 1e-18   #In KERRBB Units
    
    return Mdd_Edd

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

#SHIFTED LOG NORMAL PDF PARAMETERS (h# CORRESPONDS TO PARAMETER IN XSPEC)
#h1 = loc = shift
#h2 = ln(scale) = ln(exp(mu)) = mu
#h3 = sigma
h1 = [5.736, 1.256, 1.743, 3.475, 4.494]

h2 = [1.829, 1.542, 1.593, 1.166, 1.367]

h3 = [0.932, 0.159, 1.165, 0.188, 0.192]

#EMPTY LISTS FOR VALUES BASED ON 10TH & 90TH CENTILE MASSES
L_Edd_array = []
Mdd_Edd_array = []
Mass_90_array = []

L_Edd_array_10 = []
Mdd_Edd_array_10 = []
Mass_10_array = []

#FOR EACH XRB CALCULATE THE EDDINGTON MASS ACCRETION RATE
#FOR THE 10TH AND 90TH CENTILE MASSES
for n in range(len(h1)):    
    ln_hist = np.random.lognormal(h2[n], h3[n], 1000000)
    ln_hist = ln_hist + h1[n]
    
    percentile_90 = np.percentile(ln_hist, 90)
    Mass_90 = percentile_90
    L_Edd_var = L_Edd(Mass_90)
    Mdd_Edd_var = Mdd_Edd(L_Edd_var)
    
    percentile_10 = np.percentile(ln_hist, 10)
    Mass_10 = percentile_10
    L_Edd_var_10 = L_Edd(Mass_10)
    Mdd_Edd_var_10 = Mdd_Edd(L_Edd_var_10)
    
    L_Edd_array.append(L_Edd_var)
    Mdd_Edd_array.append(Mdd_Edd_var)
    Mass_90_array.append(Mass_90)
    
    L_Edd_array_10.append(L_Edd_var_10)
    Mdd_Edd_array_10.append(Mdd_Edd_var_10)
    Mass_10_array.append(Mass_10)
    
    #PRINT VALUES TO THE CONSOLE
    print(XRB[n],
          "\n", "90TH PERCENTILE MASS (%.4f SOLAR MASSES)" %percentile_90,
          "\n", "%s" %L_Edd_var, "W",
          "\n", "%.4f" %Mdd_Edd_var, "10e18 g s-1",
          "\n", "%.4f" %Mdd_Edd_var, "KERRBB UNITS",
          "\n",
          "\n", "10TH PERCENTILE MASS (%.4f SOLAR MASSES)" %percentile_10,
          "\n", "%s" %L_Edd_var_10, "W",
          "\n", "%.4f" %Mdd_Edd_var_10, "10e18 g s-1",
          "\n", "%.4f" %Mdd_Edd_var_10, "KERRBB UNITS",
          "\n-----------------------------------------"
          )

#PRINT VALUES TO A TEXT FILE
with open('Mdd_lims.txt', 'w') as f:
    f.write('Mdd Min & Max Values\n')
    
    for n in range(len(h1)):
        f.writelines("\n-----------------------------------------\n")
        f.writelines('\n%s: %s\n    L_Edd: %s W\n    M_Edd: %s 10e18 g s-1 (90th percentile mass)\n\n    L_Edd: %s W\n    M_Edd: %s 10e18 g s-1 (10th percentile mass)\n' %(n + 1, XRB[n], L_Edd_array[n], Mdd_Edd_array[n],L_Edd_array_10[n], Mdd_Edd_array_10[n]))

        f.writelines("\n    MAX Mdd: %.4f KERRBB UNITS" %Mdd_Edd_array[n])
        f.writelines("\n    MIN Mdd: %.11f KERRBB UNITS" %(1e-6 * Mdd_Edd_array_10[n]))    
        f.writelines("\n    90th PERCENTILE MASS: %.4f Msun\n    10th PERCENTILE MASS: %.4f Msun\n" %(Mass_90_array[n],Mass_10_array[n]))
   
#PRINT TIMER
end = datetime.datetime.now()
print(end - start)