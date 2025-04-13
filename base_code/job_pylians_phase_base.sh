#!/usr/bin/env bash
#
#
# ==== SLURM part (resource manager part) ===== #
#   Modify the following options based on your job's needs.
#   Remember that better job specifications mean better usage of resources,
#   which then means less time waiting for your job to start.
#   So, please specify as many details as possible.
#   A description of each option is available next to it.
#   SLURM cheatsheet:
# 
#     https://slurm.schedmd.com/pdfs/summary.pdf
# 
#
# ---- Metadata configuration ----
#
#SBATCH --job-name=pranjal_pylians_512       # The name of your job, you'll se it in squeue.
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL). Sends you an email when the job begins, ends, or fails; you can combine options.
#SBATCH --mail-user=pralegan@sissa.it    # Where to send the mail
#
# ---- CPU resources configuration  ----  |  Clarifications at https://slurm.schedmd.com/mc_support.html
#
#SBATCH -N 1

#SBATCH --ntasks-per-node=8

#SBATCH --hint=nomultithread

#SBATCH --mem=30000mb

#SBATCH --partition=regular1

#SBATCH --time=00:05:00
#SBATCH --output=%x.o%j              # Standard output log in TORQUE-style -- WARNING: %x requires a new enough SLURM. Use %j for regular jobs and %A-%a for array jobs
#SBATCH --error=%x.e%j               # Standard error  log in TORQUE-style -- WARNING: %x requires a new enough SLURM. Use %j for regular jobs and %A-%a for array jobs
#
# ==== End of SLURM part (resource manager part) ===== #
#
#
# ==== Modules part (load all the modules) ===== #
#   Load all the modules that you need for your job to execute.
#   Additionally, export all the custom variables that you need to export.
#   Example:
# 
#module load openmpi/1.8.3/intel/14.0
#module load fftw/3.3.4/gnu/4.9.2
#module load gsl/1.16/gnu/4.9.2
module load python3
module load texlive/2019
module load gsl/2.2/intel/19.0.4.243

#module load fftw/2.1.5/intel/14.0
#module load scipy/0.15.1/numpy/1.9.1/intel/14.0/mkl/11.1/python/2.7.8
#module load numpy/1.9.1/intel/14.0/mkl/11.1/python/2.7.8
#module load python
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/tronconi/software/hdf5-1.6.5_lib/lib
module load intel/2021.2

export PYTHONPATH=$PYTHONPATH:/home/viel/Pylians3/build/lib.linux-x86_64-3.6
export PYTHONPATH="${PYTHONPATH}:/path/to/COLIBRI/"
export PYTHONPATH="${PYTHONPATH}:/home/viel/.local/lib/python3.6/site-packages/latex"
#
# ==== End of Modules part (load all the modules) ===== #
#
#
# ==== Info part (say things) ===== #
#   DO NOT MODIFY. This part prints useful info on your output file.
#
NOW=`date +%H:%M-%a-%d/%b/%Y`
echo '------------------------------------------------------'
echo 'This job is allocated on '$SLURM_JOB_CPUS_PER_NODE' cpu(s)'
echo 'Job is running on node(s): '
echo  $SLURM_JOB_NODELIST
echo '------------------------------------------------------'
echo 'WORKINFO:'
echo 'SLURM: job starting at           '$NOW
echo 'SLURM: sbatch is running on      '$SLURM_SUBMIT_HOST
echo 'SLURM: executing on cluster      '$SLURM_CLUSTER_NAME
echo 'SLURM: executing on partition    '$SLURM_JOB_PARTITION
echo 'SLURM: working directory is      '$SLURM_SUBMIT_DIR
echo 'SLURM: current home directory is '$(getent passwd $SLURM_JOB_ACCOUNT | cut -d: -f6)
echo ""
echo 'JOBINFO:'
echo 'SLURM: job identifier is         '$SLURM_JOBID
echo 'SLURM: job name is               '$SLURM_JOB_NAME
echo ""
echo 'NODEINFO:'
echo 'SLURM: number of nodes is        '$SLURM_JOB_NUM_NODES
echo 'SLURM: number of cpus/node is    '$SLURM_JOB_CPUS_PER_NODE
echo 'SLURM: number of gpus/node is    '$SLURM_GPUS_PER_NODE
echo '------------------------------------------------------'
#
# ==== End of Info part (say things) ===== #
#

# Should not be necessary anymore with SLURM, as this is the default, but you never know...
cd $SLURM_SUBMIT_DIR


# ==== JOB COMMANDS ===== #
#   The part that actually executes all the operations you want to do.
#   Just fill this part as if it was a regular Bash script that you want to
#   run on your computer.
#   Example:
# 
cd /scratch/pralegan/ngenic_pmf/phase_ics/NGENIC-B_Pks_phase
gmake
mpirun -np 8 ./N-GenIC parameterfiles/temp.param
cd ..
mkdir pylians_output/temp
python3 powerspectra_run.py

# ==== END OF JOB COMMANDS ===== #


# Wait for processes, if any.
echo "Waiting for all the processes to finish..."
wait
