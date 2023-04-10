#SR_TAU_CELL
.libPaths("/data/user/jfisher7/.conda/envs/SR_TAU_CELL/lib/R/library")

#make sure the environment is clean
rm(list=ls())

#set the working directory 
dir_path <- "/data/project/lasseigne_lab/JLF_scratch/230321_JLF_Sex_bias_adverse_events/"
setwd(dir_path)

# load in the packages and functions
source(paste0(dir_path, "src/sex_bias_drug_functions.R"))
library(MASS)
library(parallel)
library(BiocParallel)
library(stringr)
library(dplyr)

# load in the data 
patient_safety <- readRDS( paste0(dir_path, "data/patient_safety_filtered_pt.rds"))
drug_ae_df2 <- readRDS(paste0(dir_path, "data/drug_ae_df.rds"))

#formatting data for function

patient_safety_males <- patient_safety[patient_safety$gender == 1, ]
patient_safety_females <- patient_safety[patient_safety$gender == 2, ]

additional_col <- matrix(nrow = nrow(drug_ae_df2), ncol=8)

drug_ae_df <- cbind(drug_ae_df2, additional_col)

#get arguements from the slurm job submission
args <- commandArgs(trailingOnly = TRUE)

print("slurm_inputs")
num1<- args[1] #array job number
print(num1)


#Working with only 15 samples at time
N <- nrow(drug_ae_df)
sz <- 3000 #small enough for the short queue
sample_list <- split(1:N, ceiling(seq_along(1:N)/sz))
sample_vector <- sample_list[[num1]] 

print("samples")
print(sample_vector)

#sex_bias_adverse_event_test
print(date())
parallel_tests <- mclapply(sample_vector, sex_bias_adverse_event_test, drug_ae_df, mc.cores= 10)
print(date())


group_res <- t(as.data.frame(parallel_tests))
file<- paste0(dir_path, "data/fisher_array_results/fishers_res_group_", num1, ".rds")
saveRDS(group_res, file)

sessionInfo()