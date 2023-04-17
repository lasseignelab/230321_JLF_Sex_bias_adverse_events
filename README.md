# Sex Bias Adverse Events Project

Future Abstract

## Authors

*** note to authors: This is prelim list of author. nothing is final or decided **
- [Jennifer L. Fisher](https://www.github.com/JenFisher7), [Amanda D. Clark](https://github.com/adc0032), [Emma Jones](https://github.com/emmafjones) , & [Brittany N. Lasseigne](https://github.com/blasseigne)

Department of Cell, Developmental and Integrative Biology, Heersink School of Medicine, University of Alabama at Birmingham, Birmingham, AL, 35294, USA.

## Dockers and Conda Environments

In addition to the scripts here, the Docker images used for this analysis is publicly available on Docker Hub ([jenfisher7/rstudio_sex_bias_drugs](https://hub.docker.com/r/jenfisher7/rstudio_sex_bias_drugs). For the network calculations, a conda environment was used (SR_TAU_CELL_environment.yml). Below are insturctions on how to set up the dockers on local machines and Cheaha and conda enviroments on Cheaha. This workflow assumes that you cloned the github for this project. 


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
conda create -n SR_TAU_CELL -f <path/to/SR_TAU_CELL_environment.yml>
````

Correct the path your conda environment R libraries in the 230213_fares_fisher_test.R script.
````
#SR_TAU_CELL
.libPaths("/data/user/<YOU>/.conda/envs/SR_TAU_CELL/lib/R/library")
````
To run the scripts, you use the commands in 230213_fares_fisher_test.txt in the src directory. This will run the .sh file which will use the R script. With the .sh script, please adjust lines 3, 14 to be your email and directories. In .R file, please adjust lines 2 and 8 to point to the right directores.
These scripts run as an array job so you can run several jobs in parrallel. I suggest running this on the weekend when you don't need cheaha. It is about one day run time. 

## Some saved output might be in the older verison of this project on Cheaha
````
/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/
````


## Scripts

Note for all .Rmd files there is a knitted .html file with outputs of code. Some of the code takes awhile to run (i.e., more than ~1 hour) so the ouput is not in the knit file. 

**GTEx data download and exploration**
- 210806_GTEx_Age_Sex.Rmd

 - This is important script. Update Jen on any issues with script. 

**Network analysis scripts**
- 230103_network_qsmooth.Rmd
  
  - 12 hour-short job on cheaha. important script. Update Jen on any issues with script. 
  
- 230104_panda_network_inputs.Rmd

  - 2 hour-express job on cheaha. important script. Update Jen on any issues with script. 
  
- 230105_panda_set_up.Rmd

   - 12 hour-short job on cheaha. important script. Update Jen on any issues with script. 

- 230109_sex_specific_tissue_networks.Rmd

   - 48 hour-medium job on cheaha (more time than you need). important script. Update Jen on any issues with script. 
 
- 230130_alpaca_sex_tissues.Rmd
   - 12 hour-short job on cheaha. important script. Update Jen on any issues with script. 
- 230206_alpaca_results_exploration.Rmd
- 230221_liver_sex_specific_targeting.Rmd
- 230417_liver_drug_edges.Rmd

**FARES scripts**
- 220722_Meddra_mapping.Rmd
- 230124_fares_sbae_drugs.Rmd
  - 2 hour express ; 3 cpu and 80GB per cpu in cheaha.
<br>
The following 230213_fares_fisher_test scripts use the conda enviroment. 1 day run on cheaha. Suggest running over the weekend. 

- 230213_fares_fisher_test.txt
- 230213_fares_fisher_test.sh
- 230213_fares_fisher_test.R
<br>

- 230216_fisher_res_exploration.Rmd
- 230308_drugbank_info.Rmd
- 230313_sbae_drug_target.Rmd
- 230323_sbae_heatmaps.Rmd

## Function Scripts
- sex_bias_drug_functions.R


## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="150" height="150">

## Funding

List project funding sources.

## Acknowledgements

 - We would like to thank the all the Lasseigne Lab for providing feedback during this project.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.


## Package versions and computing system information

add in the future 

**Old scripts**
remove in the future but noting here for now
- src/vignettes
- src/210603_open_fda_download.sh
- src/220519_FDA_open_drug.sh
- src/220725_FDA_open_drug_long.sh
- src/221031_transfer_sex_diff_test.Rmd
- src/230105_lioness_function_script.R
- src/230105_lioness_function_script.sh
- src/230106_lioness_tissue_runs.txt
- src/230112_lioness_function_script_smaller.R
- src/230112_lioness_function_script_smaller.sh
- src/230126_lioness_outputs_adjustment.Rmd
- src/230126_lioness_outputs_adjustment.html
- src/230130_differential_edges_testing.Rmd
- src/sex_bias_functions.R
- src/test_LV_differences.R

