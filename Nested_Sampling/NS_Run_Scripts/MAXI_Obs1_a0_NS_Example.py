import os
import xspec
import bxa.xspec as bxa
from mpi4py import MPI
import scipy.stats

#CUSTOM PRIORS
from Custom_Priors.truncated_gamma import trunc_gamma
from Custom_Priors.truncated_isotropic_sin import trunc_isotropic_sin
from Custom_Priors.truncated_lognorm import trunc_lognorm

#######################
#SAME PRIORS EVERY TIME
#######################

#TRUNCATED GAMMA PRIOR ON COLOR CORRECTION FACTOR
def transform_gamma_hd(quantile):
        alpha = 5.442
        beta = 2.6129
        
        shape = alpha
        scale = 1/beta

        return trunc_gamma(scale = scale, shape = shape, x_low = 1,x_high = 10).ppf(quantile)
        
#TRUNCATED GAMMA PRIOR ON EMISSIVITY INDEX OF INNER DISC
def transform_gamma_I1(quantile):
        alpha = 13.23
        beta = 4.0767
        
        shape = alpha
        scale = 1/beta
        
        return trunc_gamma(scale = scale, shape = shape, x_low = 0,x_high = 10).ppf(quantile)

#TRUNCATED GAMMA PRIOR ON EMISSIVITY INDEX OF OUTER DISC
def transform_gamma_I2(quantile):
        alpha = 13.23
        beta = 4.0767
        
        shape = alpha
        scale = 1/beta
        
        return trunc_gamma(scale = scale, shape = shape, x_low = 0,x_high = 10).ppf(quantile)

###################################### 
#DIFFERENT PRIORS FOR EACH SYSTEM: RGS 
######################################  

def transform_truncated_normal_OIII(quantile):
		myclip_a = 0
		myclip_b = 100
		loc = 0.0231
		scale = 0.0391
		a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale   
		
		return scipy.stats.truncnorm(a,b,loc = loc,scale = scale).ppf(quantile)
		
def transform_truncated_normal_NeII(quantile):
		myclip_a = 0
		myclip_b = 100
		loc = 0.0014
		scale = 0.0047
		a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale   
		
		return scipy.stats.truncnorm(a,b,loc = loc,scale = scale).ppf(quantile)
		
def transform_truncated_normal_NeIII(quantile):
		myclip_a = 0
		myclip_b = 100
		loc = 0.018
		scale = 0.0692
		a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale   
		
		return scipy.stats.truncnorm(a,b,loc = loc,scale = scale).ppf(quantile)
    
################################################
#DIFFERENT PRIORS FOR EACH SYSTEM: SYSTEM PRIORS
################################################     

#SHIFTED LOGNORMAL PRIOR ON MASS
def transform_shifted_log_norm_Mbh(quantile):
        s = 0.192
        loc = 4.4943
        scale = 3.922
        
        return scipy.stats.lognorm(s, loc, scale).ppf(quantile)

#TRUNCATED GAMMA PRIOR ON DISTANCE
def transform_gamma_Dbh(quantile):
        alpha = 31.6603
        beta = 10.6091
        
        shape = alpha
        scale = 1/beta
        
        return trunc_gamma(scale = scale, shape = shape, x_low = 0,x_high = 30).ppf(quantile)

#######################
#THE REST OF THE SCRIPT
#######################

#TELL XSPEC TO BE QUIET
xspec.Xset.chatter = 0

#LOAD RELXILL
xspec.AllModels.lmod("relxill", "/home/s/sdd9/relxill_v2_1")

#PATH TO THE DIRECTORY CONTAINING THIS FILE
path = os.getcwd()

#ASSUMING THE SPECTRAL FILES ARE LOCATED IN A SUB-DIRECTORY CALLED 'SPECTRA'
os.chdir('spectra/rev3531')

#USED TO CHANGE BETWEEN ML FRIENDS SAMPLER AND STEP SAMPLER
#SAFE = ML FRIENDS SAMPLER (e.g. SEE BELOW)
#NUMBER = STEP SAMPLER (e.g. speed = 80)
speed = 'safe'   

#LOAD THE DATA AND IGNORE THE GIVEN ENERGIES
xspec.AllData("rev3531_cor_pileup_PN_spectrum_grp.fits")
xspec.AllData.ignore(f"**-0.8 1.75-2.35 10.0-**")

xspec.AllModels.setEnergies("extend","low,0.01,50 log")
xspec.Fit.statMethod = "cstat"

#CREATE A MODEL 
model = xspec.Model("constant * ismabs (simpl * kerrbb + relxillCp )")   

#SET THE NAME FOR THE MODEL & OTHER USEFUL PARAMETERS (THIS IS ALSO THE FOLDER NS WILL SAVE TO)
telescopes = "EPIC"
obs = "rev3531"
spin = "a0"
step_size = "ML" if speed == "safe" else f"{speed}steps"
number_live_points = 400
filename_extras = "EXAMPLE"
modelname = f"MAXI_J1820p070_{telescopes}_{obs}_{spin}_{step_size}_{number_live_points}lps_{filename_extras}"

#################
#DEFINE THE MODEL 
################# 
xspec.AllModels(1).constant.factor.values = (1, -1)

model.ismabs.H.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.He_II.values = (0, -1)
model.ismabs.C_I.link=f"{model.ismabs.H.index}*2.40e-4*1e6"
model.ismabs.C_II.values = (0, -1)
model.ismabs.C_III.values = (0, -1)
model.ismabs.N_I.link=f"{model.ismabs.H.index}*7.59e-5*1e6"
model.ismabs.N_II.values = (0, -1)
model.ismabs.N_III.values = (0, -1)
model.ismabs.O_I.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.O_II.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.O_III.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.Ne_I.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.Ne_II.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.Ne_III.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.Mg_I.link=f"{model.ismabs.H.index}*2.51e-5*1e6"
model.ismabs.Mg_II.values = (0, -1)
model.ismabs.Mg_III.values = (0, -1)    
model.ismabs.Si_I.link=f"{model.ismabs.H.index}*1.86e-5*1e6"
model.ismabs.Si_II.values = (0, -1)
model.ismabs.Si_III.values = (0, -1)    
model.ismabs.S_I.link=f"{model.ismabs.H.index}*1.23e-5*1e6"
model.ismabs.S_II.values = (0, -1)
model.ismabs.S_III.values = (0, -1)    
model.ismabs.Ar_I.link=f"{model.ismabs.H.index}*2.57e-6*1e6"
model.ismabs.Ar_II.values = (0, -1)
model.ismabs.Ar_III.values = (0, -1)        
model.ismabs.Ca_I.link=f"{model.ismabs.H.index}*1.58e-6*1e6"
model.ismabs.Ca_II.values = (0, -1)
model.ismabs.Ca_III.values = (0, -1)
model.ismabs.Fe.values = (0.10941, 0.001, 0, 0, 100000, 1e+06)
model.ismabs.redshift.values = (0, -1)

model.simpl.Gamma.values = (1.67375, 0.05, 1, 1.1, 4, 5)
model.simpl.FracSctr.values = (0.685, 0.005,0 ,0, 0.4,1)
model.simpl.UpScOnly.values = (1, -1)

model.kerrbb.eta.values = (0, -1)
model.kerrbb.a.values = (0, -1, -0.998, -0.998, 0.998, 0.998)
model.kerrbb.i.values = (62.9188, 0.01, 3, 3, 85, 85)
model.kerrbb.Mbh.values = (8.28397, 0.01, 0, 0, 50, 50)
model.kerrbb.Mdd.values = (0.00488566,0.01,1.85592e-05,1.85592e-05,23.3447,23.3447)
model.kerrbb.Dbh.values = (2.85122,0.01,0,0,30,30)
model.kerrbb.hd.values = (2.21167,0.01,1,1,10,10)
model.kerrbb.rflag.values = (1)
model.kerrbb.lflag.values = (1)
model.kerrbb.norm.values = (1, -1)

model.relxillCp.Incl.link=f"{model.kerrbb.i.index}"
model.relxillCp.a.link=f"{model.kerrbb.a.index}"
model.relxillCp.Rin.values = (-1, -1)
model.relxillCp.Rout.values = (1000, -1)
model.relxillCp.Rbr.values = (18, -1)
model.relxillCp.Index1.values = (1, 0.1, 0, 0, 10, 10)
model.relxillCp.Index2.values = (1, 0.1, 0, 0, 10, 10)
model.relxillCp.z.values = (0, -1)
model.relxillCp.gamma.values = (1.59663,0.01,1.2,1.2,3.4,3.4)
model.relxillCp.logxi.values = (0.698748,0.01,0,0,4.7,4.7)
model.relxillCp.logN.values = (20, -1)
model.relxillCp.Afe.values = (10,0.01,0.5,0.5,10,10)
model.relxillCp.kTe.values = (400, -1)
model.relxillCp.refl_frac.values = (-1, -1)
model.relxillCp.norm.values = (0.000950156,0.01,1.67457e-09,1.67457e-09,0.167457,0.167457)

# --- SETTING PRIORS ---

prior_H = bxa.create_gaussian_prior_for(model, model.ismabs.H, 0.0902, 0.0157)
prior_O_I = bxa.create_gaussian_prior_for(model, model.ismabs.O_I, 80.0559,2.3499)
prior_O_II = bxa.create_gaussian_prior_for(model, model.ismabs.O_II, 2.9786,0.745)
prior_O_III = bxa.create_custom_prior_for(model, model.ismabs.O_III, transform_truncated_normal_OIII)
prior_Ne_I = bxa.create_gaussian_prior_for(model, model.ismabs.Ne_I, 9.6998,1.0126)
prior_Ne_II = bxa.create_custom_prior_for(model, model.ismabs.Ne_II, transform_truncated_normal_NeII)
prior_Ne_III = bxa.create_custom_prior_for(model, model.ismabs.Ne_III, transform_truncated_normal_NeIII)
prior_Fe = bxa.create_gaussian_prior_for(model, model.ismabs.Fe, 3.5528,0.5709)


prior_gamma_S = bxa.create_uniform_prior_for(model, model.simpl.Gamma)
prior_FrSc =bxa.create_uniform_prior_for(model, model.simpl.FracSctr)
                                     
prior_i = bxa.create_gaussian_prior_for(model, model.kerrbb.i, 63, 3)
prior_Mbh = bxa.create_custom_prior_for(model, model.kerrbb.Mbh, transform_shifted_log_norm_Mbh)               
prior_Mdd = bxa.create_loguniform_prior_for(model, model.kerrbb.Mdd)
prior_Dbh = bxa.create_custom_prior_for(model, model.kerrbb.Dbh, transform_gamma_Dbh)
prior_hd = bxa.create_custom_prior_for(model, model.kerrbb.hd, transform_gamma_hd)
   
prior_I1 = bxa.create_custom_prior_for(model, model.relxillCp.Index1, transform_gamma_I1)
prior_I2 =  bxa.create_custom_prior_for(model, model.relxillCp.Index2, transform_gamma_I2)                     

prior_gamma_R = bxa.create_uniform_prior_for(model, model.relxillCp.gamma)
prior_logxi = bxa.create_uniform_prior_for(model, model.relxillCp.logxi)
prior_Afe = bxa.create_uniform_prior_for(model, model.relxillCp.Afe)
prior_norm_R = bxa.create_loguniform_prior_for(model, model.relxillCp.norm)

priors = [
    prior_H, 
    prior_O_I, 
    prior_O_II, 
    prior_O_III, 
    prior_Ne_I, 
    prior_Ne_II, 
    prior_Ne_III, 
    prior_Fe,
    prior_gamma_S, 
    prior_FrSc, 
    prior_i, 
    prior_Mbh, 
    prior_Mdd,
    prior_Dbh, 
    prior_hd,
    prior_I1, 
    prior_I2, 
    prior_gamma_R, 
    prior_logxi, 
    prior_Afe, 
    prior_norm_R
]

#OUTPUT FILENAME
outputfiles_basename = f"{modelname}_fit"

#DEFINE THE SOLVER
solver = bxa.BXASolver(transformations=priors, outputfiles_basename=outputfiles_basename)

#RUN THE NESTED SAMPLER 
results = solver.run(resume=True, speed=speed, n_live_points=number_live_points)

#ON ONLY ONE CORE DO:
if MPI.COMM_WORLD.Get_rank() == 0 or MPI.COMM_WORLD.Get_size() == 1 :
            #SET THE BEST FIT IN XSPEC 
            solver.set_best_fit()
            if os.path.isfile(f"{outputfiles_basename}.xcm") : os.remove(f"{outputfiles_basename}.xcm")
            if os.path.isfile(f"{outputfiles_basename}_model.xcm") : os.remove(f"{outputfiles_basename}_model.xcm")
            #SAVE XCM FILES
            xspec.Xset.save(f"{outputfiles_basename}.xcm", info="a")
            xspec.Xset.save(f"{outputfiles_basename}_model.xcm", info="m")