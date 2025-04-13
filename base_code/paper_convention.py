# In[ ]:
# import necessary modules
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d, PchipInterpolator
import math
from scipy.integrate import solve_ivp
import pandas as pd
from numpy import log, exp
Gamma=math.gamma

def table(f,t):
    tablef=np.zeros(t.shape[0])
    for i in np.arange(t.shape[0]):
        tablef[i]=f(t[i])
    return tablef

from background import rhodm, rhob, rhor, Ht, alpha, photon_mfp_full, cb2_full, a_rec, a_eq, a_background, rho_tot
omegam=0.12038+0.022032
omegab=0.022032/omegam
omegadm=0.12038/omegam

Hsim=Ht#lambda a: Ht(0.01)*np.sqrt(rhob(a)+rhodm(a)+rhor(a))/np.sqrt(rhob(0.01)+rhodm(0.01)+rhor(0.01))
rho_tot2=lambda a:rhob(a)+rhodm(a)+rhor(a)
rhom=lambda a:rhob(a)+rhodm(a)
#In[]
#scale invariant simplified perturbation eqn post recombination
def ddel_baryon(a,delta): #DM and baryon perturbation equations
        dldm,thdm,dlb,thb=delta[0],delta[1],delta[2],delta[3]
        k2phi=-3./2*(a*Hsim(a))**2*(rhob(a)*dlb+rhodm(a)*dldm)/rho_tot(a)
        del_dldm=-(thdm/a**2/Hsim(a))
        del_thdm=-thdm/a+k2phi/a**2/Hsim(a)
        del_thb=-thb/a+k2phi/a**2/Hsim(a)+Hsim(a)*rhom(a)/rho_tot(a)*1 #the S0 term captures the impact of PMF on delta_b. S0 is taken to be a function of kd, which in turn is a function of a.
        del_dlb=-(thb/a**2/Hsim(a))

        return [del_dldm,del_thdm,del_dlb,del_thb]

def ddel_baryon_grav(a,delta,k): #DM and baryon perturbation equations
        phi,dldm,thdm,dlb,thb=delta[0],delta[1],delta[2],delta[3],delta[4]
        del_phi=-phi/a+(-k**2*phi-3./2*(a*Ht(a))**2*(rhob(a)*dlb+rhodm(a)*dldm)/rho_tot(a))/(3*(a*Ht(a))**2)/a
        del_dldm=-(thdm/a**2/Ht(a) - 3*del_phi)
        del_thdm=-thdm/a+k**2/a**2/Ht(a)*phi
        del_thb=-(1+alpha(a)/Ht(a))*thb/a+k**2/a**2/Ht(a)*phi+Ht(a)*0 #the S0 term captures the impact of PMF on delta_b. S0 is taken to be a function of kd, which in turn is a function of a.
        del_dlb=-(thb/a**2/Ht(a) - 3*del_phi)

        return [del_phi,del_dldm,del_thdm,del_dlb,del_thb]

delta_b0=0
delta_dm0=0
theta_b0=0#theta_b_pmf[i_100]/(norm)/0.817
theta_dm0=0#-(delta_dm_pmf[i_100]-delta_dm_pmf[i_100-1])/(a_pmf[i_100]-a_pmf[i_100-1])*0.01**2*Ht(0.01)/(norm)/1.62

deltasolve0=[delta_dm0,theta_dm0,delta_b0,theta_b0] #initial condition assumes everythong zero as perturbations are sourced by PMFs. I start my computation from a_mfp at which point silk damping should have erasedd any earlier evolution.
deltasolve=solve_ivp(lambda a,y: ddel_baryon(a,y),[a_rec,1],deltasolve0,method='BDF',dense_output=True,atol=1e-6,rtol=1e-5)

xi_b=lambda a: deltasolve.sol(a)[2]
xi_dm=lambda a: deltasolve.sol(a)[0]
xi_m=lambda a: omegab*deltasolve.sol(a)[2]+omegadm*deltasolve.sol(a)[0]
fb_ratio=lambda a: deltasolve.sol(a)[2]/(omegab*deltasolve.sol(a)[2]+omegadm*deltasolve.sol(a)[0])
#In[]:
#matter power spectrum from PMFs
def param(nB):
    if nB==-2.5:
        GnB=0.446299
        kappa_b=1.47/1.8*1.6
        kappa_dm=kappa_b/2
        kappa_m=kappa_b*0.7
        p=1.9
    elif nB==-2:
        GnB=1.178097
        kappa_b=1.57 #This specifies the precise location of Jeans scale for baryons
        kappa_dm=0.77 #This specifies the precise location of Jeans scale for dark matter
        kappa_m=1.17 #This specifies the precise location of Jeans scale for total matter
        p=1.95
    elif nB==-2.9:
        GnB=0.0684258
        kappa_b=0.92/1.8*1.6
        kappa_dm=kappa_b/2
        kappa_m=kappa_b*0.7
        p=1.9
    elif nB==-2.95:
        GnB=0.0284859
        kappa_b=0.92/1.8*1.6
        kappa_dm=kappa_b/2
        kappa_m=kappa_b*0.7
        p=1.9
    elif nB==-1.6:
        GnB=5.185
        kappa_b=1.6
        kappa_dm=kappa_b/2
        kappa_m=kappa_b*0.7
        p=1.9
    elif nB==-1.8:
        GnB=1.89
        kappa_b=1.6
        kappa_dm=kappa_b/2
        kappa_m=kappa_b*0.7
        p=1.9
    elif nB==-2.2:
        GnB=0.816046
        kappa_b=1.6
        kappa_dm=kappa_b/2
        kappa_m=kappa_b*0.7
        p=1.9
    else:
        GnB=5
        kappa_b=1.6
        kappa_dm=kappa_b/2
        kappa_m=kappa_b*0.7
    return [GnB,kappa_b,kappa_dm,kappa_m,p]
k1mpc=1/0.678

#PS0_an=lambda k,B1mpc,nB,eta: (2*np.pi**2/k**3)*10**-4*(k/k1mpc)**(2*nB+10)*B1mpc**4*G(nB)*exp(-2*k**2/kd(B1mpc,nB,eta)**2)
def PS0_an(k,B1mpc,nB):
    [GnB,kappa_b,kappa_dm,kappa_m,p]=param(nB)
    ans=(2*np.pi**2/k**3)*0.918*10**-4*(k/k1mpc)**(2*nB+10)*B1mpc**4*GnB
    return ans

def P_fit_nB2(K,kj,p,gamma):
    ans=gamma*(K/kj)**7*(1+(K/kj)**p)**(-7/p)
    return ans

def Pb_an(a,k,B1mpc,nB):
    if nB<-1.5:
        [GnB,kappa_b,kappa_dm,kappa_m,p]=param(nB)
        kj=(0.1*kappa_b*B1mpc)**(-2/(nB+5))*k1mpc
        ans=xi_b(a)**2*PS0_an(k,B1mpc,nB)*(1+(k/kj)**p)**(-(2*nB+10)/p)
    if nB==2:
        gamma=0.73
        kj=(0.1*1.13*B1mpc)**(-2/(nB+5))*k1mpc
        ans=P_fit_nB2(k,kj,2.5,gamma)*(xi_b(a)/xi_b(0.01))**2
    return ans

def Pdm_an(a,k,B1mpc,nB):
    if nB<-1.5:
        [GnB,kappa_b,kappa_dm,kappa_m,p]=param(nB)
        kj=(0.1*kappa_dm*B1mpc)**(-2/(nB+5))*k1mpc
        ans=xi_dm(a)**2*PS0_an(k,B1mpc,nB)*(1+(k/kj)**p)**(-(2*nB+10)/p)
    if nB==2:
        gamma=0.215
        kj=(0.1*0.334*B1mpc)**(-2/(nB+5))*k1mpc
        ans=P_fit_nB2(k,kj,2.5,gamma)*(xi_dm(a)/xi_dm(0.01))**2
    return ans

def Pm_an(a,k,B1mpc,nB):
    if nB<-1.5:
        [GnB,kappa_b,kappa_dm,kappa_m,p]=param(nB)
        kj=(0.1*kappa_m*B1mpc)**(-2/(nB+5))*k1mpc
        ans=xi_m(a)**2*PS0_an(k,B1mpc,nB)*(1+(k/kj)**p)**(-(2*nB+10)/p)
    if nB==2:
        gamma=0.214
        kj=(0.1*0.596*B1mpc)**(-2/(nB+5))*k1mpc
        ans=P_fit_nB2(k,kj,2.5,gamma)*(xi_m(a)/xi_m(0.01))**2
    return ans

kd=lambda B1mpc,nB: (0.1*0.8*B1mpc)**(-2/(nB+5))*k1mpc
PB_an=lambda k,B1mpc,nB: 10**-18*4*np.pi**2/Gamma((nB+3)/2)*B1mpc**2/k1mpc**3*(k/k1mpc)**(nB)*exp(-k**2/kd(B1mpc,nB)**2)

M_D=lambda B1mpc, nB:5*10**12*(1/kd(B1mpc,nB))**3 #provides the damping mass. Our results breakdown for masses below this.

def findB(nB, Box, Grid):
    knyq=np.pi*Grid/Box #in h/Mpc
    B1mpc_ans=(1.5*knyq)**(-(nB+5)/2)*10/0.8
    return B1mpc_ans

def findBox(nB, Box, Grid):
    knyq=np.pi*Grid/Box #in h/Mpc
    B1mpc_ans=(1.5*knyq)**(-(nB+5)/2)*10/0.8
    return B1mpc_ans

#import colibri.cosmology as cc
#C = cc.cosmo()
Pk_m=np.loadtxt('/scratch/pralegan/ngenic_pmf/PM_lcdm.dat')
K_i=Pk_m[:,0]
PM_i=Pk_m[:,1]

def Pm_tot(a,K,B1mpc,nB):
    #kk,pk_class=C.class_Pk(k=K,z=1/a-1)
    pk_extrapolate=np.exp(np.interp(np.log(K),np.log(K_i),np.log(PM_i)))*(a/0.01)**2
    if (B1mpc==0):
        pk_pmf=0*pk_extrapolate
    else:
        pk_pmf=table(lambda x: Pm_an(a,x,B1mpc,nB), K)
    return pk_pmf+pk_extrapolate#pk_class[0]

#HMF
def HMF(zz,B1mpc,nB):
    size=105
    start=-2
    if (B1mpc==0):
        end=2
    else:
        end=log(5*kd(B1mpc,nB))/log(10)
    step=(end-start)/size
    K=10**np.arange(start,end,step)
    logM = np.arange(5.1,15.,0.1)
    
    HMF_pmf = C.halo_mass_function(logM = logM, k = K, pk = 1*Pm_tot(1/(1+zz),K,B1mpc,nB), window = 'th', mass_fun = 'ShethTormen')

    temp=np.stack([10.**logM,10.**logM*HMF_pmf[0]], axis=0)
    
    return temp

#Extrapolate Pk with gravity from z=100
def Pks_z(zz,Pb_i,PDM_i,K):#Pbi_i and PDM_i are power specctra at z=100, zz is the redshifts at which I want the output
    PDM=np.zeros([zz.shape[0],K.shape[0]])
    Pb=np.zeros([zz.shape[0],K.shape[0]])
    Pm=np.zeros([zz.shape[0],K.shape[0]])
    j=0
    for k in K:
        delta_b0=(k**3*Pb_i[j]/(2*np.pi**2))**0.5
        delta_dm0=(k**3*PDM_i[j]/(2*np.pi**2))**0.5
        phi_0=-3/2*(0.01*Ht(0.01))**2/k**2*(rhob(0.01)*delta_b0+rhodm(0.01)*delta_dm0)/rho_tot(0.01)
        theta_b0=-0.01*Ht(0.01)*delta_b0
        theta_dm0=-0.01*Ht(0.01)*delta_dm0

        deltasolve0=[phi_0,delta_dm0,theta_dm0,delta_b0,theta_b0] #initial condition assumes everythong zero as perturbations are sourced by PMFs. I start my computation from a_mfp at which point silk damping should have erasedd any earlier evolution.
        deltasolve=solve_ivp(lambda a,y: ddel_baryon_grav(a,y,k),[0.01,1],deltasolve0,method='BDF',dense_output=True,atol=1e-6,rtol=1e-5)
        
        l=0
        for z in zz:
            if z==0:
                a_z=0.99
            else:
                a_z=1/(1+z)

            i_z=np.where(deltasolve.t>a_z)[0][0]
            PDM[l,j]=deltasolve.y[1,i_z]**2
            Pb[l,j]=deltasolve.y[3,i_z]**2
            Pm[l,j]=((0.022*deltasolve.y[3,i_z]+0.12*deltasolve.y[1,i_z])/0.142)**2
            l=l+1
        j=j+1
    return [Pb,PDM,Pm]
# %%
