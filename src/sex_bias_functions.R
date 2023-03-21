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
<<<<<<< HEAD
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
=======
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
>>>>>>> bf8151362ff122e4bbf42c52f3a5c20deb8c0bcf
  
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
  
  #divide the male and female samples 
  male<- processed_expression[, metadata$gtex.sex == "M"]
  female <- processed_expression[, metadata$gtex.sex == "F"]
  
  print("creating male network")
  male_network <- panda(expr = as.data.frame(male), motif = motif_table, ppi = ppi_net, progress=TRUE, mode="intersection")
  saveRDS(male_network, male_network_file)
  
  print("creating female network")
  female_network <- panda(expr = as.data.frame(male), motif = motif_table, ppi = ppi_net, progress=TRUE, mode="intersection")
  saveRDS(female_network, female_network_file)
}