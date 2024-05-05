#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-briansmithers        #specify the account to use
#SBATCH --job-name=blister-projections            # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
##SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32              # number of cores to allocate
#SBATCH --mem=512G                     # 2000 MB of Memory allocated; set --mem with care
##SBATCH --mem-per-cpu=100G
#SBATCH --time=1-00:00:00                 # Maximum job run time
##SBATCH --array=0-1847%10                  # Number of jobs in array
#SBATCH --mail-user=shuysman@gmail.com
#SBATCH --mail-type=ALL
##SBATCH --output=fire-%A-%a.out
##SBATCH --error=fire-%A-%a.err
#SBATCH --output=blister-projections-%j.out
#SBATCH --error=blister-projections-%j.err

date
hostname -s
module load Mamba
source $HOME/.bashrc
mamba activate r-spatial
Rscript sdm_projection_ensemble.R
date
