#In[]
import Pk_library as PKL
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

folder_name='temp'
folder_sims='ics'
fullpath_snap = os.path.join(folder_sims, folder_name, 'ics')

folder_pylians='pylians_output'
fullpath_pylians = os.path.join(folder_pylians, folder_name)

#In[]
# parameters                                                                                                                                                                                   
snapshot      = fullpath_snap #snapshot name                                                                                                                                        
grid          = 512    #grid size                                                                                                                                                              
particle_type = [0,1]    #use dark matter [1]                                                                                                                                                
do_RSD        = False   #move particles to redshift-space and calculate Pk in redshift-space                                                                                                   
axis          = 1      #RSD placed along the y-axis                                                                                                                                            
cpus          = 8      #number of openmp threads                                                                                                                                               
folder_out    = fullpath_pylians #folder where to write results                                                                                                                                            

# compute power spectrum of the snapshot                                                                                                                                                       
PKL.Pk_Gadget(snapshot, grid, particle_type, do_RSD, axis, cpus, folder_out)

#In[]
#Get back the outputs Pks
Pk_dm_temp=pd.read_fwf(fullpath_pylians+'/Pk_CDM_z=99.000.dat', header=None).values
Pk_dm=np.column_stack((np.asarray(Pk_dm_temp[:,0], dtype=float),np.asarray(Pk_dm_temp[:,1], dtype=float)))

Pk_b_temp=pd.read_fwf(fullpath_pylians+'/Pk_GAS_z=99.000.dat', header=None).values
Pk_b=np.column_stack((np.asarray(Pk_b_temp[:,0], dtype=float),np.asarray(Pk_b_temp[:,1], dtype=float)))

Pk_m_temp=pd.read_fwf(fullpath_pylians+'/Pk_GAS+CDM_z=99.000.dat', header=None).values
Pk_m=np.column_stack((np.asarray(Pk_m_temp[:,0], dtype=float),np.asarray(Pk_m_temp[:,1], dtype=float)))

Pk_dmb_temp=pd.read_csv(fullpath_pylians+'/Pk_GCDM_z=99.000.dat',delim_whitespace=True, header=None).values
Pk_dmb=np.column_stack((np.asarray(Pk_dmb_temp[:,0], dtype=float),np.asarray(Pk_dmb_temp[:,1], dtype=float)))

pk_b2_temp=pd.read_csv('NGENIC-B_Pks_phase/CAMB/Pb_lcdm.dat',sep='\t').values
pk_b2=np.column_stack((np.asarray(pk_b2_temp[:,0], dtype=float),np.asarray(pk_b2_temp[:,1], dtype=float)))
pk_dm2_temp=pd.read_csv('NGENIC-B_Pks_phase/CAMB/PDM_lcdm.dat',sep='\t').values
pk_dm2=np.column_stack((np.asarray(pk_dm2_temp[:,0], dtype=float),np.asarray(pk_dm2_temp[:,1], dtype=float)))
pk_m2_temp=pd.read_csv('NGENIC-B_Pks_phase/CAMB/PM_lcdm.dat',sep='\t').values
pk_m2=np.column_stack((np.asarray(pk_dm2_temp[:,0], dtype=float),np.asarray(pk_dm2_temp[:,1], dtype=float)))

#In[] analytical estimate

current_dir = os.path.abspath(os.getcwd()) #get current directory
# Add base_code directory to sys.path
sys.path.append(os.path.join(current_dir, 'base_code'))

from paper_convention import Pb_an, Pdm_an, Pm_an, table, kd, xi_b, xi_dm, xi_m

B1mpc=0.013
nB=-2.0
lcdm=1

Pb_an_table1=table(lambda x: Pb_an(0.01,x,B1mpc,nB),pk_b2[:,0])*(3/xi_b(0.01))**2
Pb_an_table2=Pb_an_table1#+0.0*2*np.pi**2/pk_b2[:,0]**3*(exp(-2*kd(1,-2.9,0.1)**2/pk_b2[:,0]**2))
pb_tot=lcdm*pk_b2[:,1]+Pb_an_table2

Pdm_an_table1=table(lambda x: Pdm_an(0.01,x,B1mpc,nB),pk_dm2[:,0])*(0.25/xi_dm(0.01))**2
Pdm_an_table2=Pdm_an_table1#+0.0*2*np.pi**2/pk_dm2[:,0]**3*(exp(-2*kd(1,-2.9,0.3)**2/pk_dm2[:,0]**2))*(0.217/2.43)**2
pdm_tot=lcdm*pk_dm2[:,1]+Pdm_an_table2

Pm_an_table1=table(lambda x: Pm_an(0.01,x,B1mpc,nB),pk_dm2[:,0])*(0.67/xi_m(0.01))**2
Pm_an_table2=Pm_an_table1#+0.0*2*np.pi**2/pk_dm2[:,0]**3*(exp(-2*kd(1,-2.9,0.3)**2/pk_dm2[:,0]**2))*(0.217/2.43)**2
pm_tot=lcdm*pk_m2[:,1]+Pm_an_table2

#In[]:
#Plot
tempk=Pk_b[:-1,0].astype("float")
tempp=Pk_b[:-1,1].astype("float")
tempk2=Pk_dm[:-1,0].astype("float")
tempp2=Pk_dm[:-1,1].astype("float")

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
plt.rcParams.update({'font.size': 16})
miny=10**-5
maxy=10

minx=kd(B1mpc,nB)/20
maxx=kd(B1mpc,nB)*20

# plt.loglog(Pk_b1[:-1,0], Pk_b1[:-1,0]**3*Pk_b1[:-1,1]/(2*np.pi**2),'-', color='k', linewidth=2,label=r"100 x 256")
ax1.loglog(tempk, tempk**3*tempp/(2*np.pi**2),'-', color='darkred', linewidth=2,label=r"Ng baryon")
ax1.loglog(tempk2, tempk2**3*tempp2/(2*np.pi**2),'-', color='k', linewidth=2,label=r"Ng DM")
# plt.loglog(Pk_b3[:,0], Pk_b3[:,0]**3*Pk_b3[:,1]/(2*np.pi**2),'-', color='darkorange', linewidth=2,label=r"50 x 512")
ax1.loglog(pk_b2[:,0], pk_b2[:,0]**3*pb_tot/(2*np.pi**2),'--', color='darkorange', linewidth=2,label=r"an baryon")
ax1.loglog(pk_dm2[:,0], pk_dm2[:,0]**3*pdm_tot/(2*np.pi**2),'--', color='blue', linewidth=2,label=r"an DM")

ax2.loglog(Pk_m[:,0], Pk_m[:,0]**3*Pk_m[:,1]/(2*np.pi**2),'-', color='k', linewidth=2,label=r"Ngenic matter")
ax2.loglog(Pk_dmb[:,0], Pk_dmb[:,0]**3*Pk_dmb[:,1]/(2*np.pi**2),'-', color='darkred', linewidth=2,label=r"Ngenic GCDM")
# plt.loglog(Pk_b3[:,0], Pk_b3[:,0]**3*Pk_b3[:,1]/(2*np.pi**2),'-', color='darkorange', linewidth=2,label=r"50 x 512")
ax2.loglog(pk_m2[:,0], pk_m2[:,0]**3*pm_tot/(2*np.pi**2),'--', color='blue', linewidth=2,label=r"analytical matter")

ax1.loglog([0.01,100],[1,1],'--',color='gray')
ax1.loglog([kd(B1mpc,nB),kd(B1mpc,nB)],[miny,maxy],'--',color='gray')
ax1.text(1.1*kd(B1mpc,nB),miny*1.2,"$k_D$")
ax1.set_xlabel("$k$ (in h/Mpc)")
ax1.set_ylabel(r"$k^3P(k)/[2\pi^2]$",labelpad=0)
ax1.set_ylim(miny, maxy)
ax1.set_xlim(minx, maxx)
ax1.legend(ncol=2)

ax2.loglog([0.01,100],[1,1],'--',color='gray')
ax2.set_xlabel("$k$ (in h/Mpc)")
ax2.set_ylim(miny, maxy)
ax2.set_xlim(minx, maxx)
ax2.legend()

#plt.gcf().subplots_adjust(left=0.12, right=0.96, bottom=0.12, top=0.97)
fig.savefig(fullpath_pylians+'/Pks.pdf',format="pdf", bbox_inches="tight")