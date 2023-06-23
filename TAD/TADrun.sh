#!/bin/bash -l
#
# -- Request that this job run on whichJet
#SBATCH --partition=tjet,ujet,sjet,vjet,xjet,kjet
# -- Specifying the project account
#SBATCH -A vdatum
#SBATCH -J hello
#SBATCH -c 1
#
# -- Specify where to put stdout and stderr
#SBATCH -o hello-%a.out
#SBATCH -D .
#PBS -M elena.tolkova@noaa.gov
#
# -- Specify a maximum wallclock of 1 hours
#SBATCH --time=00:30:00
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=HPL
#SBATCH -v
#SBATCH --array=0-23
#SBATCH --qos batch 
#SBATCH --mem=30000M
#SBATCH --ntasks=1
date
echo "host name is `hostname` : $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT ../fort.63.nc"
cp   ../TXstabilize/fort.63.nc tt/${SLURM_ARRAY_TASK_ID}.nc
./TADcalc_mp.py $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT tt/${SLURM_ARRAY_TASK_ID}.nc
date
