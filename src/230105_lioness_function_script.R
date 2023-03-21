.libPaths("/data/user/jfisher7/.conda/envs/netzoo/lib/R/library")

#make sure the environment is clean
rm(list=ls())

#set the working directory 
setwd("/data/project/lasseigne_lab/JLF_scratch/JLF_Aim1_prelim_data_disease/")

# load in the packages 
library(pandaR)

#protein-protein networks
ppi_net <- readRDS("/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/panda_input_data/230104_symbol_numeric_ppi_network_11_5_0.rds")
#transcription factor and gene networks based on motif info
motif_table <- readRDS("/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/panda_input_data/230104_symbol_motif_table.rds")

#get arguements from the slurm job submission
args <- commandArgs(trailingOnly = TRUE)

print("slurm_inputs")
num1<- args[1] #array job number
print(num1)
tissue_file <- args[2] # file path to the tissue panda network 
print(tissue_file)

expression_file <- args[3] #file path to the tissue proccessed gene expression data
print(expression_file)

tissue <- args[4] #file path to save the lioness output
print(tissue)

#read in files
tissue_panda<- readRDS(tissue_file)

expression <- as.data.frame(readRDS(expression_file))


#Working with only 15 samples at time
N <- ncol(expression)
sz <- 15 #small enough for the short queue
sample_list <- split(1:N, ceiling(seq_along(1:N)/sz))
sample_vector<- sample_list[[num1]] 

print("samples")
print(sample_vector)

#lioness function from the netZooR github from 230105
lion_res <- lapply(sample_vector, function(i) {
        print(paste("Computing network for sample ", i))
         N * tissue_panda@regNet - (N - 1) * panda(motif = motif_table, expr = expression[, -i],  ppi =ppi_net ,progress=TRUE, mode="intersection" )@regNet
        })

lion_file<- paste("/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/output/lioness/lioness_jobs_output/", tissue, "_", num1, "_lioness_res.rds", sep = "")
saveRDS(lion_res, lion_file )

sessionInfo()