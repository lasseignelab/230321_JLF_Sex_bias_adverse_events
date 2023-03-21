#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=jfisher7@uab.edu
#SBATCH --job-name=lioness
#SBATCH -n 1
#SBATCH --mem-per-cpu=50000
#SBATCH --nodes=1
#SBATCH --time=0-12:00:00  
#SBATCH --share 
#SBATCH --partition=short 
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.out.txt

cd /data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/

module load intel/2017a
module load Anaconda3/5.3.1
source activate
conda activate netzoo

echo $SLURM_ARRAY_TASK_ID 
echo $2
echo $3
echo $4

Rscript ./script/230112_lioness_function_script_smaller.R $SLURM_ARRAY_TASK_ID $2 $3 $4