.libPaths("/data/user/jfisher7/.conda/envs/pca_test/lib/R/library")

setwd("/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/FDA_OPEN/")
#make sure the environment is clean
rm(list=ls())

library(rjson)
library(dplyr)

files <- list.files(path= "/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/FDA_OPEN/" , pattern = ".json")
files_2 <- files[!grepl(".zip", files)]
files_3 <- files_2[!grepl(".sh", files_2)]


###############Function to use###################################################
FDA_JSON_TO_DATA_FRAME<- function(index, paths= files_3, num_drug= 3000){
  #path- path to the json file
  #NUM_DRUG the number of drugs a patent can be on
  #year of FDA cases
  #Quarter of FDA cases
  path<- paths[index]
  result <- fromJSON(file = path)
  print(path)
  # Convert JSON file to a data frame.2\
  #json_data_frame <- as.data.frame(result$meta)
  leng <- length(result$results)
  year <- strtrim(path, 4)
  quarter<- substring(path, 6, 6)
  data<- c("case_numer","year_Q", "receive_date", "source_qualification", "country", "sex", "age", "number_drugs", "serious", "drug", "reaction")
  for (i in 1:leng){ #leng
    case_id<- paste0("case:", i)
    print(case_id)
    case <- result$results[[i]]
    year_Q <- paste(year,"_",quarter, sep="")
    
    rec_date<- case$receivedate
    rec_date <- ifelse(is.null(rec_date), "NA", rec_date)
    
    source_qual <- case$primarysource$qualification
    source_qual <- ifelse(is.null(source_qual), "NA", source_qual)
    
    source_country<- case$primarysource$reportercountry
    source_country <- ifelse(is.null(source_country), "NA", source_country)

    t_F <- "patientsex" %in% names(case$patient)
    #t_F <- as.data.frame(t_F)
    #fact <- ifelse(is.na(t_F[2,2]), FALSE, t_F[2,2] > 1)
    if (  t_F ){
      #patent info
      sex <- ifelse(case$patient$patientsex == "2", "F", "M")
      t_F <- "patientonsetage" %in% names(case$patient)
      if (t_F){
        age <- case$patient$patientonsetage
      }else{
        age<- NA
      }
      reaction <- c()
      serious<- case$serious
      #start for loop for reactions 
      for (k in 1:length(case$patient$reaction)){
        reac <- case$patient$reaction[[k]]$reactionmeddrapt
        reaction<- c(reaction, reac)
        
        
        #start for loop for drugs 
        if (length(case$patient$drug) <= num_drug){
          drug <- c()
          app_number<- c()
          brand_name<- c()
          generic_name <- c()
          manufacturer_name <- c()
          product_type<-c()
          substance_name <- c()
          pharm_epc <- c()
          pharm_moa <- c()
          pharm_unii <- c()
         
          for (j in 1:length(case$patient$drug)){
            drug[j]<- case$patient$drug[[j]]$medicinalproduct
            
            #get the fda_open info 
            #app_number
            t_F <- "application_number" %in% names(case$patient$drug[[j]]$openfda)
             if (t_F){
                app_number[j] <- case$patient$drug[[j]]$openfda$application_number
            }else{
                app_number[j]<- NA
            }
            #brand_name<- c()
            t_F <- "brand_name" %in% names(case$patient$drug[[j]]$openfda)
             if (t_F){
                brand_name[j] <- case$patient$drug[[j]]$openfda$brand_name
            }else{
                brand_name[j]<- NA
            }
            #generic_name <- c()
            t_F <- "generic_name" %in% names(case$patient$drug[[j]]$openfda)
            if (t_F){
                generic_name[j] <- case$patient$drug[[j]]$openfda$generic_name
            }else{
                generic_name[j]<- NA
            }
            #manufacturer_name <- c()
            t_F <- "manufacturer_name" %in% names(case$patient$drug[[j]]$openfda)
            if (t_F){
                manufacturer_name[j] <- case$patient$drug[[j]]$openfda$manufacturer_name
            }else{
                manufacturer_name[j]<- NA
            }
            #product_type<-c()
            t_F <- "product_type" %in% names(case$patient$drug[[j]]$openfda)
            if (t_F){
                product_type[j] <- case$patient$drug[[j]]$openfda$product_type
            }else{
                product_type[j]<- NA
            }
            #substance_name <- c()
            t_F <- "substance_name" %in% names(case$patient$drug[[j]]$openfda)
            if (t_F){
                substance_name[j] <- case$patient$drug[[j]]$openfda$substance_name
            }else{
                substance_name[j]<- NA
            }
            #pharm_epc <- c()
            t_F <- "pharm_class_epc" %in% names(case$patient$drug[[j]]$openfda)
            if (t_F){
                pharm_epc[j] <- case$patient$drug[[j]]$openfda$pharm_class_epc
            }else{
                pharm_epc[j]<- NA
            }
            #pharm_moa <- c()
            t_F <- "pharm_class_moa" %in% names(case$patient$drug[[j]]$openfda)
            if (t_F){
                pharm_moa[j] <- case$patient$drug[[j]]$openfda$pharm_class_moa
            }else{
                pharm_moa[j]<- NA
            }
            #pharm_unii
            t_F <- "unii" %in% names(case$patient$drug[[j]]$openfda)
            if (t_F){
                pharm_unii[j] <- case$patient$drug[[j]]$openfda$unii
            }else{
                pharm_unii[j]<- NA
            }
            
            
            drug_info<- cbind(drug,app_number,brand_name, generic_name,manufacturer_name, product_type, substance_name, pharm_epc, pharm_moa, pharm_unii)
        }
          #add to dataframe 
          
          #case_number <- paste("case_", i, sep="")
          #print(case_number)
          drug_reaction_df <- expand.grid(drug, reaction )
          colnames(drug_info)[1]<- "Var1"
          #print(colnames(drug_info))
          #print(colnames(drug_reaction_df ))
          
          drug_reaction_df <- as.data.frame(drug_reaction_df) %>% inner_join(as.data.frame(drug_info))
          #print(drug_reaction_df )
          #print("trying line 61")
          case_number1<- rep(case$safetyreportid[1], nrow(drug_reaction_df))
          #print("trying line 63")
          year_Q<- rep(year_Q[1], nrow(drug_reaction_df))
          
          #rec_date
          rec_date <- rep(rec_date[1], nrow(drug_reaction_df))
          #source_qual
          source_qual <- rep(source_qual[1], nrow(drug_reaction_df))
          #source_country
          source_country <- rep(source_country[1], nrow(drug_reaction_df))
         # print("trying line 65")
          sex<- rep(sex[1], nrow(drug_reaction_df))
          #print("trying line 67")
          age<- rep(age[1], nrow(drug_reaction_df))
          #print("trying line 69")
          number_drugs<- rep(length(drug), nrow(drug_reaction_df))
          #print("trying line 70")
          serious<- rep(serious[1], nrow(drug_reaction_df) )
          #reaction<- rep(unique(reaction), length(unique(drug)))
          #drug_v2<- rep(unique(drug), number_reaction )
          
          
          #print("trying line 72")
          
           an.error.occured <- FALSE
          tryCatch( { data1 <- cbind(case_number1, year_Q, rec_date,source_qual, source_country,  sex, age, number_drugs) }, error = function(e) {an.error.occured <<- TRUE})
          if(an.error.occured){
            print("error") 
            print(nrow(drug_reaction_df))
            print(case_number1)
            print(year_Q)
            print(rec_date)
            print(source_qual)
            print(source_country)
            print(age)
            print(sex)
            print(number_drugs)
            print(serious)
          }else{
            #print(source_country)
            data1 <- cbind(case_number1, year_Q, rec_date,source_qual, source_country,  sex, age, number_drugs, serious)
            data1<- cbind(data1, drug_reaction_df)
          }
          #print(data1)
          an.error.occured <- FALSE
          tryCatch( { data <- rbind(data, data1) }, error = function(e) {an.error.occured <<- TRUE})
          if(an.error.occured){
            print("error") 
            print(ncol(data))
            print(ncol(data1))
            data <- rbind(data, data1)
          }else{
            data <- rbind(data, data1)
          }
          
        }else{}
        
        #medicinalproduct 
        #drug info 
        #brand_name
      }
    }
    
  }
  #there are duplicate rows
  data<- data[!duplicated(data),]
  data<- data[-1,]
  rownames(data)<- NULL
  print("Done")
  return(data)
}


############### 1 DRUG ###############

#save2 <- mclapply(1:2, FDA_JSON_TO_DATA_FRAME(index, paths=files_3, num_drug=3000), mc.cores = cores)

args <- commandArgs(trailingOnly = TRUE)
print(args)
num1<- as.numeric(args)

file_name<- paste0("/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/FDA_OPEN/220727_FARES_res_", num1, ".rds")

d1_res <- FDA_JSON_TO_DATA_FRAME(index= num1, paths=files_3, num_drug=3000)

saveRDS(d1_res, file_name)
  

print(sessionInfo())
