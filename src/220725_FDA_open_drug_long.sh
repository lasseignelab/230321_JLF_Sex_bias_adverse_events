#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=jfisher7@uab.edu
#SBATCH --job-name=FARES_long
#SBATCH -n 1
#SBATCH --mem-per-cpu=10000 
#SBATCH --nodes=1
#SBATCH --time=6-06:00:00  
#SBATCH --share 
#SBATCH --partition=long
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.out.txt
cd /data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/script/
module load intel/2017a
module load R/3.6.2-foss-2018a-X11-20180131-bare

Rscript drugeventscript_drug.R  $SLURM_ARRAY_TASK_ID