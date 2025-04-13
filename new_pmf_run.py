import os
import sys
#In[]
#Choose PMF parameters
B1mpc=0.1 #in nanoGauss
nB=-2.0

same_phase=0 #if set to one then PMF and LCDM perturbations have the same phase. That is no isocurvature

baryon_dm_same=0 #if set to 1 then assumes baryon power =dark matter power=matter power

lcdm=0 #if set to 0 then turn off lcdm initial conditions

#Choose box parameters
Box=1 #in Mpc/h
grid=256

mesh=3 #2 mesh means initial glassfile only has baryon and DM. setting to 3 means that initial glass file has baryon,DM and neutrinos, but we remove neutrinos at the end.

email='pralegan@sissa.it' #specify the email where ulysses will send notification of jobs

if grid%64==0:
    fac=int(grid/64)
else:
    raise ValueError("grid size is not multiple of 64. My whole code requires grid size to be a multiple of 64.")

#In[]
############### Derived parameters for different nB values ##########################
current_dir = os.path.abspath(os.getcwd()) #get current directory
# Add base_code directory to sys.path
sys.path.append(os.path.join(current_dir, 'base_code'))

from paper_convention import param
[GnB,kappa_b,kappa_dm,kappa_m,p]=param(nB)

kappa_b=kappa_b
if baryon_dm_same==1:
    kappa_b=kappa_m
    kappa_dm=kappa_m

coeff=0.918*10**-4*GnB #This is the coefficient in front of dimensionless 3Mpl^2*SB/rhom power spectrum.

#In[]#############################################################################################################################################
#I read a basic parameter file and then edit it for non-vanilla changes
with open('NGENIC-B_Pks_phase/parameterfiles/base.param', 'r') as file:
    # Read the entire file
    content = file.read()

modified_content=content.replace('B1mpc            1.0','B1mpc            '+str(B1mpc))

modified_content=modified_content.replace('nB    		 -2.9','nB    		 '+str(nB))

modified_content=modified_content.replace('coeff    	6.281e-6','coeff    	'+str(coeff))

modified_content=modified_content.replace('kappa_b    	0.92','kappa_b    	'+str(kappa_b))

modified_content=modified_content.replace('kappa_dm    	0.46','kappa_dm    	'+str(kappa_dm))

modified_content=modified_content.replace('p_shape    	2.5','p_shape    	'+str(p))

modified_content=modified_content.replace('Box              55000.0','Box              '+str(Box*1000))

#filename is the name of the parameter file in which the new parameters will be saved
filename='B_'+str(B1mpc).replace('.','_')+'_nB_'+str(-nB).replace('.','_')+'_box_'+str(Box)+'_'+str(grid)+'_phasePK'

if lcdm==0:
    modified_content=modified_content.replace('lambda_lcdm_b    1','lambda_lcdm_b    0')
    modified_content=modified_content.replace('lambda_lcdm_dm   1','lambda_lcdm_dm   0')
    filename=filename+'_no_lcdm'

if same_phase==1 or baryon_dm_same==1:
    modified_content=modified_content.replace('same_phase	0','same_phase	1')
    filename=filename+'_same'

if baryon_dm_same==1:
    modified_content=modified_content.replace('lambda_pmf_b     3','lambda_pmf_b     0.67')
    modified_content=modified_content.replace('lambda_pmf_dm    0.25','lambda_pmf_dm     0.67')
    modified_content=modified_content.replace('FileWithInputSpectrumB   ./CAMB/Pb_lcdm.dat','FileWithInputSpectrumB   ./CAMB/PM_lcdm.dat')
    modified_content=modified_content.replace('FileWithInputSpectrumDM   ./CAMB/PDM_lcdm.dat','FileWithInputSpectrumDM   ./CAMB/PM_lcdm.dat')
    filename=filename+'_matter'

if grid!=512:
    modified_content=modified_content.replace('Nmesh           512','Nmesh           '+str(grid))
    modified_content=modified_content.replace('Nsample         512','Nsample         '+str(grid))
    modified_content=modified_content.replace('GlassTileFac      8','GlassTileFac      '+str(fac))
    if grid<256 or grid==1024:
        modified_content=modified_content.replace('NumFilesWrittenInParallel 16','NumFilesWrittenInParallel '+str(fac))

if mesh==3:
    modified_content=modified_content.replace('GlassFile         ./glassfiles/dummy_glass_2comp_64.dat','GlassFile         ./glassfiles/dummy_glass_CDM_B_NU_64_64_64.dat')
    modified_content=modified_content.replace('OmegaDM_2ndSpecies  0.00','OmegaDM_2ndSpecies  0.001')
    modified_content=modified_content.replace('NU_On                0','NU_On                1')
    filename=filename+'_3mesh'

if mesh==1: #This needs to be debugged; currently does not work!!!!!!!!!!!!!
    modified_content=modified_content.replace('GlassFile         ./glassfiles/dummy_glass_2comp_64.dat','GlassFile         ./glassfiles/dummy_glass_CDM_B_NU_64_64_64.dat')
    modified_content=modified_content.replace('OmegaDM_2ndSpecies  0.00','OmegaDM_2ndSpecies  0.001')
    modified_content=modified_content.replace('NU_On                0','NU_On                1')
    filename=filename+'_1mesh'

#writing the modified parameters into a file
with open('NGENIC-B_Pks_phase/parameterfiles/'+filename+'.param', 'w') as file:
    # Write some initial content to the file
    file.write(modified_content)

#In[]###############################################################################################################################################
#Now I edit the powerspectra file
with open('base_code/powerspectra_base.py', 'r') as file:
    # Read the entire file
    content2 = file.read()

modified_content2=content2.replace("folder_name='temp'","folder_name='"+filename+"'")
if grid!=512:
    modified_content2=modified_content2.replace("grid          = 512","grid          = "+str(grid))
    modified_content2=modified_content2.replace("cpus          = 8","cpus          = "+str(fac))

modified_content2=modified_content2.replace("B1mpc=0.013","B1mpc="+str(B1mpc))
modified_content2=modified_content2.replace("nB=-2.0","nB="+str(nB))
modified_content2=modified_content2.replace("lcdm=1","lcdm="+str(lcdm))

#writing the modified parameters into a file
with open('powerspectra_run.py', 'w') as file:
    # Write some initial content to the file
    file.write(modified_content2)
#In[]#########################################################################################################################################################################
#Now I edit the bash script job file

with open('base_code/job_pylians_phase_base.sh', 'r') as file:
    # Read the entire file
    content3 = file.read()

current_dir = os.path.abspath(os.getcwd()) #get current directory
modified_content3=content3.replace('#SBATCH --job-name=pranjal_pylians_512','#SBATCH --job-name='+filename)
modified_content3=modified_content3.replace('#SBATCH --mail-user=pralegan@sissa.it','#SBATCH --mail-user='+email)

current_dir = os.path.abspath(os.getcwd()) #get current directory
modified_content3=modified_content3.replace('cd /scratch/pralegan/ngenic_pmf/phase_ics/NGENIC-B_Pks_phase','cd '+current_dir+'/NGENIC-B_Pks_phase')

if grid==256:
    modified_content3=modified_content3.replace('#SBATCH --ntasks-per-node=8','#SBATCH --ntasks-per-node=4')
    modified_content3=modified_content3.replace('#SBATCH --mem=30000mb','#SBATCH --mem=5000mb')
    modified_content3=modified_content3.replace('#SBATCH --time=00:05:00','#SBATCH --time=00:03:00')
    modified_content3=modified_content3.replace('mpirun -np 8 ./N-GenIC parameterfiles/temp.param','mpirun -np '+str(fac)+' ./N-GenIC parameterfiles/'+filename+'.param')
elif grid==512:
    modified_content3=modified_content3.replace('mpirun -np 8 ./N-GenIC parameterfiles/temp.param','mpirun -np 8 ./N-GenIC parameterfiles/'+filename+'.param')
elif grid==1024:
    modified_content3=modified_content3.replace('#SBATCH -N 1','#SBATCH -N 2')
    modified_content3=modified_content3.replace('#SBATCH --mem=30000mb','#SBATCH --mem=120000mb')
    modified_content3=modified_content3.replace('#SBATCH --time=00:05:00','#SBATCH --time=00:20:00')
    modified_content3=modified_content3.replace('mpirun -np 8 ./N-GenIC parameterfiles/temp.param','mpirun -np '+str(fac)+' ./N-GenIC parameterfiles/'+filename+'.param')
else:
    raise ValueError("For chosen grid size need to add custom bash script parameters")

modified_content3=modified_content3.replace('mkdir pylians_output/temp','mkdir pylians_output/'+filename)

with open('job_pylians_phase_run.sh', 'w') as file:
    # Write some initial content to the file
    file.write(modified_content3)

#In[]:#############################################################################################################################################################3
#for 1024 because we run in cambridge would like to turn on 32bit in makefile
if grid==1024:
    with open('NGENIC-B_Pks_phase/Makefile', 'r') as file:
    # Read the entire file
        content4 = file.read()
    newcontent4=content4.replace('OPT    +=  -DNO64BITID','#OPT    +=  -DNO64BITID')
    if newcontent4!=content4:
        with open('NGENIC-B_Pks_phase/Makefile', 'w') as file:
        # Read the entire file
            file.write(newcontent4)
else:
    with open('NGENIC-B_Pks_phase/Makefile', 'r') as file:
    # Read the entire file
        content4 = file.read()
    newcontent4=content4.replace('#OPT    +=  -DNO64BITID','OPT    +=  -DNO64BITID')
    if newcontent4!=content4:
        with open('NGENIC-B_Pks_phase/Makefile', 'w') as file:
        # Read the entire file
            file.write(newcontent4)

#In[]#####################################################################################################################################
#submit the script to run
import subprocess
subprocess.run(['sbatch', 'job_pylians_phase_run.sh'])
# %%
