#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=jfisher7@uab.edu
#SBATCH --job-name=fisher_test
#SBATCH -n 10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=0-02:00:00  
#SBATCH --nodes=1
#SBATCH --share 
#SBATCH --partition=express
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.out.txt

cd /data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/

module load intel/2017a
module load Anaconda3/5.3.1
source activate
conda activate netzoo

echo $SLURM_ARRAY_TASK_ID 


Rscript ./script/230213_fares_fisher_test.R $SLURM_ARRAY_TASK_ID