#In[]:
# importing functions
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.interpolate import interp1d
import math
from scipy.integrate import solve_ivp

def table(f,t):
    tablef=np.zeros(t.shape[0])
    for i in np.arange(t.shape[0]):
        tablef[i]=f(t[i])
    return tablef
# In[ ]:
# Cosmological parameters and other CLASS parameters##############################################
common_settings = {# we need to set the output field to something although
                   # the really releveant outpout here will be set with 'k_output_values'
                   #'output':'mPk',
                   # value of k we want to polot in [1/Mpc]
                   #'k_output_values':k,
                   # LambdaCDM parameters
                   'h':0.678,
                   'omega_b':0.022032,
                   'omega_cdm':0.12038,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'tau_reio':0.0925,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   'YHe':0.246,
                   # other options and settings
                   'compute damping scale':'yes', # needed to output the time of damping scale crossing
                   'gauge':'newtonian'}
##############
#
# call CLASS
#
M = Class()
M.set(common_settings)
M.compute()

#background quantities
background = M.get_background()
H_table=background['H [1/Mpc]']
rhodm_table=background['(.)rho_cdm']
rhob_table=background['(.)rho_b']
rhor_table=background['(.)rho_g']+background['(.)rho_ur']
rhotot_table=background['(.)rho_tot']
a_background=1/(1+background['z'])
background_R=  3.0/4*(background['(.)rho_b'])/(background['(.)rho_g'])

a_at_rho_m_over_r = interp1d((rhob_table+rhodm_table)/(rhor_table),a_background)
a_eq = a_at_rho_m_over_r(1.)

rhodm = interp1d(a_background,rhodm_table)
rhob = interp1d(a_background,rhob_table)
rhor = interp1d(a_background,rhor_table)
rho_tot = interp1d(a_background,rhotot_table)
Ht=interp1d(a_background,H_table)
R=interp1d(a_background,background_R)

quantities = M.get_current_derived_parameters(['z_rec'])
a_rec = 1/(1+quantities['z_rec'])

#thermodynamic quantities
ther=M.get_thermodynamics()
photon_mfp=1/ther["kappa' [Mpc^-1]"] #in comoving units coordinates
therm_a=1/(1+ther['z'])
cb2=ther['c_b^2']

def photon_mfp_full(a):
    if a>therm_a[-1]:
        ans=1/np.interp(1/a,1/therm_a,1/photon_mfp)
    else:
        ans=photon_mfp[-1]*(a/therm_a[-1])**2
    return ans

alpha=lambda a: 1/photon_mfp_full(a)/a/R(a)

def cb2_full(a):
    if a>therm_a[-1]:
        ans=1/np.interp(1/a,1/therm_a,1/cb2)
    else:
        ans=cb2[-1]*(therm_a[-1]/a)**1
    return ans

