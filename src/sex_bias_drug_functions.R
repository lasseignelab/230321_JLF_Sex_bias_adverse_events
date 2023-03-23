# These are the functions created for the sex bias adverse events project
#created by Jennifer Fisher

MedDRA_AE_mapping <- function(FARES){
  #inputs:
  #FARES is a data frame with the following columns from the drugevents_drug.R 
  #"case_number1" "year_Q" "sex" "age"  "number_drugs" "serious"  "Var1"   "Var2"  
  
  #output
  #FARES_AE is a data frame with "case_number", "year_Q", "sex", "age" , "number_drugs", "serious" , "drug"  , 
  #"AE_FARES", "llt_id", "llt_name", "pt_id", "pt_name", "soc_id", "soc_name", "soc_short"
  #note if the preferred is only given it doesn't have a lower level term id or name
  #also some terms from FARES are not mapping directly to a term due to phase differences.
  #these are given na for all columns 
  
  #note the meddra mapping created in the 220722_FDA_OPEN_FARES_Exploration is needed for this function
  meddrda_mapping_all<- readRDS("~/data/MedDRA_mapping_llt_pt_soc.rds")
  meddrda_mapping_all$pt_name <- toupper(meddrda_mapping_all$pt_name)
  meddrda_mapping_all$llt_name <- toupper(meddrda_mapping_all$llt_name)
  #want to make sure that the Names match for with preferred term or lower level terms
  
  index<- as.data.frame(matrix(nrow= nrow(FARES), ncol = 7))
  colnames(index)<-  c("llt_id", "llt_name", "pt_id", "pt_name", "soc_id", "soc_name", "soc_short")
  for( i in 1:nrow(FARES) ){ #nrow(FARES)
    
    # determine if it is a lower level term or preferred term 
    test_term <- FARES[i,8] %in% meddrda_mapping_all$pt_name
    
    #if it perffered term match this way
    if(test_term == TRUE){
      #sec1<- c() 
      sec2_num <- grep(FARES[i,8],meddrda_mapping_all$pt_name, fixed = TRUE )[1]
      sec2<- meddrda_mapping_all[sec2_num, 4:8]
      item <- c(NA , NA , as.character(sec2))
      #print(item)
    }
    #some term don't match for either 
    if(FARES[i,8] %in% meddrda_mapping_all$llt_name){
      #if it is lower level term match this way
      item<- as.vector(meddrda_mapping_all[FARES[i,8] == meddrda_mapping_all$llt_name, c(1:2, 4:8)])
      #print(item)
    }else{
      item <- rep(NA, 7)
    }
    #colnames(index)<- NULL 
    #names(item)<- NULL
    #rownames(index) <- NULL
    #print(length(item))
    #print(i)
    if ( ! is.na(as.vector(item)[3]) ){
      suppressWarnings(index[i,]<- as.vector(item))
    }
  }
  FARES_AE<- cbind(FARES, index)
  colnames(FARES_AE) <- c("case_number", "year_Q", "sex", "age" , "number_drugs", "serious" , "drug"  , 
                          "AE_FARES", "llt_id", "llt_name", "pt_id", "pt_name", "soc_id", "soc_name", "soc_short")
  return(FARES_AE)
}


#function to make raw counts to TPM 
Counts_to_tpm <- function(counts, featureLength) {
  #counts the data frame
  #featureLength is a vector with the gene lengths from Recount3
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}


gtex_qsmooth_function<- function(count, metadata, tissue , file_filtered_expression= "/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/gtex_qsmooth/tissue_gtex_filtered_expression_set.rds", file_norm= "/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/gtex_qsmooth/tissue_gtex_qsmooth_expression.rds", file_meta="/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/gtex_qsmooth/tissue_gtex_metadata.rds" ){
  #This function is used to gsmooth normalize based on https://academic.oup.com/biostatistics/article/19/2/185/3949169
  #which was used in the YARN netZoo R package 
  #inputs
  
  #count- the count gtex metrics with all the samples 
  #metadata- the metadata for all the samples 
  #tissue- the specific tissue for subset for the normalization 
  #file_filtered_expression- the the file path and name for the filtered expression 
  #file_norm- the file path and name for the file for the normalized expression 
  #file_meta- the file path and name for the file for the subset of metadata
  
  #output
  #filtered expression (cpm < 1 are removed)
  #normlized table of expression 
  #subsetted metadata
  
  #count- the count gtex metrics with all the samples 
  #metadata- the metadata for all the samples 
  #tissue- the specific tissue for subset for the normalization 
  #file_filtered_expression- the the file path and name for the filtered expression 
  #file_norm- the file path and name for the file for the normalized expression 
  #file_meta- the file path and name for the file for the subset of metadata
  
  #output
  #filtered expression (cpm < 1 are removed)
  #normlized table of expression 
  #subsetted metadata
  
  
  metadata_v2 <- metadata[metadata$gtex.smtsd == tissue,]
  metadata_v2<- metadata_v2[complete.cases(metadata_v2$gtex.smtsd),]
  count_tissue <- count[,colnames(count) %in% rownames(metadata_v2)]
  
  print("check the metadata and count data match")
  print(identical(colnames(count_tissue) , rownames(metadata_v2)))
  
  print("make expression set object")
  
  tissue_set <- ExpressionSet(assayData=count_tissue)
  
  print("filter lowly expressed genes")
  
  tissue_filtered <- filterLowGenes(tissue_set,groups=metadata_v2$gtex.sex)
  
  print("before")
  print(dim(tissue_set))
  print("after")
  print(dim(tissue_filtered))
  
  print("save filtered expression")
  saveRDS(tissue_filtered, file= file_filtered_expression)
  
  
  print("run and save results for qsmooth from YARN R package")
  qsmooth_res <- yarn:::qsmooth(tissue_filtered,groups=metadata_v2$gtex.sex)
  saveRDS(qsmooth_res, file= file_norm)
  
  print("save metadata")
  saveRDS(metadata_v2, file= file_meta)
}

panda_tissue_run<- function(norm_expression,
                            expression_file="/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/output/lioness/proccessed_normalized_expression",  
                            network_file= "/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/output/lioness/tissue_panda"){ 
  #this function assumes that feature info is available in environment   
  
  print("check the gene ids")
  gene_ids<- feature_info$gene_id[feature_info$gene_id %in% rownames(norm_expression)]
  gene_symbols<- feature_info$gene_name[feature_info$gene_id %in% rownames(norm_expression)]
  identical(rownames(norm_expression), gene_ids)
  rownames(norm_expression)<- gene_symbols
  
  print("remove duplicated gene ids")
  express_rm <- keep_max_duplicate(norm_expression)
  saveRDS(express_rm, expression_file )
  
  print( "run panda")
  networks <- panda(expr = as.data.frame(express_rm), motif = motif_table, ppi = ppi_net, progress=TRUE, mode="intersection")
  saveRDS(networks, network_file)
}

keep_max_duplicate<- function(expression_matrix){
  med <- rowMedians(expression_matrix)
  #check
  #print(med[order(-med)][1:5])
  #print(med[order(med)][1:5])
  expres_order<- expression_matrix[order(-med),]
  expres_rm<- expres_order[ !duplicated(rownames(expres_order)),]
  return(expres_rm)
}

sex_networks<-  function(processed_expression, 
                         metadata, 
                         male_network_file, 
                         female_network_file ){
  #processed_expression
  #metadata
  #male_network_file
  #female_network_file
  
  #check to make sure the samples are in the correct order
  print("check to make sure the samples are in the correct order")
  print(identical(colnames(processed_expression), metadata$external_id))
  
  #motif table should only have transcription factors and genes expressed in the tissue 
  motif_table_v2<- motif_table[motif_table[,1] %in% rownames(processed_expression),]
  motif_table_v2<- motif_table_v2[motif_table_v2[,2] %in% rownames(processed_expression),]
  # ppi_net should only have transcription factors and genes expressed in the tissue 
  ppi_net_v2<- ppi_net[ppi_net[,1] %in% rownames(processed_expression),]
  ppi_net_v2<- ppi_net_v2[ppi_net_v2[,2] %in% rownames(processed_expression),]
  
  #divide the male and female samples 
  male<- processed_expression[, metadata$gtex.sex == "M"]
  print(dim(male))
  female <- processed_expression[, metadata$gtex.sex == "F"]
  print(dim(female))
  
  print("creating male network")
  male_network <- panda(expr = as.data.frame(male), motif = motif_table_v2, ppi = ppi_net_v2, progress=TRUE, mode="intersection")
  saveRDS(male_network, male_network_file)
  
  print("creating female network")
  female_network <- panda(expr = as.data.frame(female), motif = motif_table_v2, ppi = ppi_net_v2, progress=TRUE, mode="intersection")
  saveRDS(female_network, female_network_file)
}


create_alpaca_input<- function(file_path_name_female, file_path_name_male, tissue_name){
  
  #make a longer dataframe for the male network
  panda_M <- readRDS(paste0(paste0(dir_path, "/output/alpaca/sex_specific_networks/"),  file_path_name_male, sep = ""))
  regNet_panda_intersection_male <- panda_M@regNet
  topNet_male <- topedges( panda_M , cutoff = 2)
  topNET_male_regnet <- as.data.frame(topNet_male@regNet)
  topNET_male_regnet$TF<- rownames(topNET_male_regnet)
  topNET_male_regnet_long <- topNET_male_regnet  %>%  pivot_longer(!TF, names_to = "gene", values_to = "male") 
  
  #make a longer dataframe for the female network
  panda_F <- readRDS(paste0(paste0(dir_path, "/output/alpaca/sex_specific_networks/"),  file_path_name_female, sep = ""))
  regNet_panda_intersection_female <- panda_F@regNet
  topNet_female <- topedges( panda_F , cutoff = 2)
  topNET_female_regnet <- as.data.frame(topNet_female@regNet)
  topNET_female_regnet$TF<- rownames(topNET_female_regnet)
  topNET_female_regnet_long <- topNET_female_regnet  %>%  pivot_longer(!TF, names_to = "gene", values_to = "female")
  
  #check the tf and genes are in the same order
  print("TF")
  print(identical(topNET_female_regnet_long$TF, topNET_male_regnet_long$TF))
  print("gene")
  print(identical(topNET_female_regnet_long$gene, topNET_male_regnet_long$gene))
  #add together
  topNet_regnet_male_female_long<- cbind(topNET_male_regnet_long, topNET_female_regnet_long$female)
  
  colnames(topNet_regnet_male_female_long)<- c("TF", "gene", "male", "female")
  
  file_name<- paste0(paste0(dir_path, "/output/alpaca/alpaca_adjusted_inputs/"),  tissue_name, "_topNet_regnet_male_female_long.rds")
  
  saveRDS(topNet_regnet_male_female_long, file_name )
}

sex_bias_adverse_event_test<- function( index, drug_event_df ){
  #event and drug
  event<- drug_event_df[index,2]
  drug<- drug_event_df[index,1]
  #a: the number of female patients with target drug-events combinations.
  a_count <- nrow(patient_safety_females[ (grepl(event,patient_safety_females$SE ) & grepl(drug,patient_safety_females$drugs )), ] )
  #b: the number of female patients with target drugs but not target events.
  b_count <- nrow(patient_safety_females[ ( !grepl(event,patient_safety_females$SE ) &  grepl(drug,patient_safety_females$drugs )), ] )
  #c: the number of male patients with target drug-events combinations.
  c_count <- nrow(patient_safety_males[ (grepl(event,patient_safety_males$SE ) & grepl(drug,patient_safety_males$drugs )), ] )
  #d: the number of male patients with target drugs but not target events.
  d_count <- nrow(patient_safety_males[ ( !grepl(event,patient_safety_males$SE ) &  grepl(drug,patient_safety_males$drugs )), ] )
  
  #contingency table
  con_table<- matrix( c(a_count, c_count, b_count, d_count), nrow= 2, dimnames = list(c("Female", "Male"), c("Target Adverse Events", "All other Adverse Events")))
  
  #fisher's exact test
  res<- fisher.test(con_table)
  
  table_info <- c(a_count, b_count, c_count, d_count, res$p.value, res$estimate, res$conf.int[1], res$conf.int[2])
  
  return(table_info)
  
}


#old functions no longer used in main project
lioness_output_adjustment<- function(index){
  #need to determine the number the lioness subsets
  num_files <- length(lioness_output_files[agrep(tissues_wo_ws[index], lioness_output_files)])
  print(num_files)
  lion_file<- paste(paste0(dir_path, "output/lioness/lioness_jobs_output/"), tissues_wo_ws[index], "_", 1:num_files, "_lioness_res.rds", sep = "")
  
  for (j in 1:num_files){
    print(lion_file[j])
    lioness_res <- readRDS(lion_file[j])
    
    
    for (i in 1:length(lioness_res)){
      
      if (j == 1 & i == 1){
        test_df <- as.data.frame(lioness_res[j])
        test_df$Gene <- rownames(test_df)
        test_long <- pivot_longer(test_df, cols=!Gene, names_to = "TF")
        rownames(test_long) <- paste(test_long$TF, test_long$Gene, sep = "_")
        inital_df <- test_long
      }else{
        
        test_df <- as.data.frame(lioness_res[i])
        test_df$Gene <- rownames(test_df)
        test_long <- pivot_longer(test_df, cols=!Gene, names_to = "TF")
        inital_df <- cbind(inital_df, test_long[,3])
      }
      
    }
    print(j)
    print("Done")
  }
  df_file<-paste(paste0(dir_path, "output/lioness/adjusted_lioness_output/"), tissues_wo_ws[index], "_lioness_res_all.rds", sep = "")
  saveRDS(inital_df, df_file)
}