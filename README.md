# Sex-biased gene expression and gene regulatory networks of sex-bias adverse event drug targets and drug metabolism genes

## Abstract
Background:

Previous pharmacovigilance studies and a retroactive review of cancer clinical trial studies have identified that women were more likely to experience drug adverse events (i.e., any unintended effects of medication), and men were more likely to experience adverse events that resulted in hospitalization or death. These sex-bias adverse events (i.e., SBAEs) are due to many factors not entirely understood, including hormones, tissue differences, pharmacokinetics/pharmacodynamics, and sex chromosome dosage effects.

Methods:
 
 Based on previous studies of sex-bias gene expression and gene-regulatory networks, we evaluated the role of sex- and tissue-specific gene expression and regulation in SBAEs. Therefore, we utilized sex- and tissue-specific gene expression, and constructed sex-specific tissue gene-regulatory networks with the Genotype-Tissue Expression Project  (GTEx) gene expression profiles. Then, we identified drugs associated with SBAEs and their targets and metabolism enzymes to determine if SBAE are associated with sex-bias gene expression or sex-specific gene-regulatory network properties and predicted regulatory relationships. 
 
Results: 

We identified liver-specific gene-regulatory differences for drug metabolism genes between males and females, which could explain observed sex differences in pharmacokinetics and pharmacodynamics. In addition, we found that ~85% of SBAE drug targets had sex-bias gene expression or were core genes of sex- and tissue-specific network communities. This number of SBAE drug targets with sex-bias gene expression was significantly higher compared to randomly selected drug targets. Additionally, the number of SBAE drug targets  that were core genes of sex- and tissue-specific network communities was also significantly higher compared to randomly selected drug targets. Lastly, we provide our results of sex-bias drug-adverse event pairs, drug targets, and drug metabolism enzymes as a resource for the research community.  

Conclusions:

Overall, we provide evidence that many  SBAEs are associated with drug targets and drug metabolism genes that are differentially expressed and regulated between males and females. These SBAE drug metabolism enzymes and drug targets could be used to help explain and potentially identify SBAEs. 


## Authors

- [Jennifer L. Fisher](https://www.github.com/JenFisher7), [Amanda D. Clark](https://github.com/adc0032), [Emma F. Jones](https://github.com/emmafjones) , & [Brittany N. Lasseigne](https://github.com/blasseigne)

Department of Cell, Developmental and Integrative Biology, Heersink School of Medicine, University of Alabama at Birmingham, Birmingham, AL, 35294, USA.

## Dockers and Conda Environments

In addition to the scripts here, the Docker images used for this analysis is publicly available on Docker Hub ([jenfisher7/rstudio_sex_bias_drugs)](https://hub.docker.com/r/jenfisher7/rstudio_sex_bias_drugs). For the Fisher's extact calculations, a conda environment was used (SR_TAU_CELL_environment.yml). At the bottom of this ReadMe are detailed insturctions on how to set up the dockers on local machines and conda enviroments on Cheaha, UAB's High Performance Computing system. This workflow assumes that you cloned the github for this project. 

## Package versions and computing system information
[Comprehensive computer and package version information for docker and conda](https://github.com/lasseignelab/230321_JLF_Sex_bias_adverse_events/blob/main/src/230512_computer_R_info.pdf)

## Scripts

Note for all .Rmd files there is a knitted .html file with outputs of code. Some of the code takes awhile to run (i.e., more than ~1 hour) so the ouput is not in the knit file. 

**GTEx data download and exploration**
- 210806_GTEx_Age_Sex.Rmd


**Network analysis scripts**
- 230103_network_qsmooth.Rmd
  
  - 12 hour-short job on Cheaha
  
- 230104_panda_network_inputs.Rmd

  - 2 hour-express job on Cheaha
  
- 230105_panda_set_up.Rmd

   - 12 hour-short job on Cheaha

- 230109_sex_specific_tissue_networks.Rmd

   - 48 hour-medium job on Cheaha (more time than you need)
   
- 230130_alpaca_sex_tissues.Rmd
   - 12 hour-short job on Cheaha
- 230206_alpaca_results_exploration.Rmd
- 230221_liver_sex_specific_targeting.Rmd
- 230417_liver_drug_edges.Rmd

**FARES scripts**
- 220722_Meddra_mapping.Rmd
- 230124_fares_sbae_drugs.Rmd
  - 2 hour express ; 3 cpu and 80GB per CPU on Cheaha.
<br>
The following 230213_fares_fisher_test scripts use the conda enviroment with a 1 day run on Cheaha. 

- 230213_fares_fisher_test.txt
- 230213_fares_fisher_test.sh
- 230213_fares_fisher_test.R
<br>

- 230216_fisher_res_exploration.Rmd
- 230308_drugbank_info.Rmd
- 230313_sbae_drug_target.Rmd
- 230323_sbae_heatmaps.Rmd
- 230426_more_drug_metabolism_testing.Rmd

## Function Scripts
- sex_bias_drug_functions.R
 - This a script with functions for the analysis scripts above. Most of the functions were written by Jennifer Fisher. However, some of the functions were developed by others. These functions have a comment with them about who wrote the function and if there were any changes to them for our specific study.


## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="150" height="150">

## Funding

- JLF and BNL were supported by R03OD030604; JLF, EFJ, and BNL were supported by R00HG009678; ADC and BNL are supported by U54OD030167; JLF, EJF, and BNL are supported by UAB funds to the Lasseigne Lab. 


## Acknowledgements

 - We would like to thank the all the Lasseigne Lab for providing feedback during this project.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.

## Additional Information for Dockers and Conda Environments 
**How to use the Docker on local mac:**

Download the Docker image from Docker Hub (only done once)
````
docker pull jenfisher7/rstudio_sex_bias_drugs
````

Adjust Docker Desktop settings to CPUs = 12 ; Memory = 64 GB; swap = 1 GB; disk size = 1.6 TB

Mount the github as your working directory to match file paths and to avoid upstream files that could affect results. 
```
docker run -d --rm -p 8787:8787 -e PASSWORD=NBI -v [Replace with your path to the github directory]/230321_JLF_Sex_bias_adverse_events:/home/rstudio/ jenfisher7/rstudio_sex_bias_drugs
```
go to chrome and enter http://localhost:8787/


**How to use the Docker with Singularity:**

Use interactive HPC deskstop on Cheaha. Jobs will require 80GB, 1 CPU, and individual scripts will note the runtime needed. 

Go to the github directory 
````
cd [Replace with your path to the github directory]/230321_JLF_Sex_bias_adverse_events
module load Singularity/3.5.2-GCC-5.4.0-2.26
````
Pull the docker image from Dockerhub and create sif file (only done once)
````
singularity pull docker://jenfisher7/rstudio_sex_bias_drugs
````

create files need to run singularity (only do this step once).
````
mkdir -p run var-lib-rstudio-server
printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf
````
Run the following to start the "docker"
````
#binds the 230321_JLF_Sex_bias_adverse_events directory to your home directory in RStudio
singularity exec --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf --bind ./ rstudio_sex_bias_drugs_latest.sif rserver --server-user=[change to your id] --www-address=127.0.0.1
````
open firefox and go to http://localhost:8787/

Note that scripts that run on cheaha have the following line of code at the top of the script:

````
dir_path <- "/data/project/lasseigne_lab/JLF_scratch/230321_JLF_Sex_bias_adverse_events/"
````
Please change to the correct path to the github directory on Cheaha. 

**How to use the Conda environment:**

This conda enviorment only applies the following scripts: 230213_fares_fisher_test.R, 230213_fares_fisher_test.sh, 230213_fares_fisher_test.txt.
230213_fares_fisher_test.R is the R script used in the 230213_fares_fisher_test.sh slurm batch script. 230213_fares_fisher_test.txt notes the slurm command used in the terminal to run the 230213_fares_fisher_test.sh array job. 

Before running this script you need to create conda enviroment.
````
module load Anaconda3/5.3.1
conda env create -n SR_TAU_CELL -f <path/to/SR_TAU_CELL_environment.yml>
````

Correct the path your conda environment R libraries in the 230213_fares_fisher_test.R script.
````
#SR_TAU_CELL
.libPaths("/data/user/<YOU>/.conda/envs/SR_TAU_CELL/lib/R/library")
````
To run the scripts, you use the commands in 230213_fares_fisher_test.txt in the src directory. This will run the .sh file which will use the R script. With the .sh script, please adjust lines 3, 14 to be your email and directories. In .R file, please adjust lines 2 and 8 to point to the right directores.
These scripts run as an array job so you can run several jobs in parrallel. It is about one day run time. 


