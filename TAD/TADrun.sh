#!/bin/bash -l
#
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
#SBATCH --time=01:00:00
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=HPL
#SBATCH -v
#SBATCH --array=0-23
#SBATCH --qos batch 
#SBATCH --mem=20000M
#SBATCH --ntasks=1
date
echo "host name is `hostname` : $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT fort.63.nc"
cp   ../EastCoast-55/fort.63.nc tt/${SLURM_ARRAY_TASK_ID}.nc
module load python/3.10.8
python TADcalc_mp.py $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT tt/${SLURM_ARRAY_TASK_ID}.nc -1 -1
date
