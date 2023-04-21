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


gtex_qsmooth_function <- function(count, metadata, tissue , file_filtered_expression= "/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/gtex_qsmooth/tissue_gtex_filtered_expression_set.rds", file_norm= "/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/gtex_qsmooth/tissue_gtex_qsmooth_expression.rds", file_meta="/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/data/gtex_qsmooth/tissue_gtex_metadata.rds" ){
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

panda_tissue_fix <- function(norm_expression,
                            expression_file="/data/project/lasseigne_lab/JLF_scratch/Sex_Bias_Adverse_Events/output/lioness/proccessed_normalized_expression"){ 
  #this function assumes that feature info is available in environment   
  
  print("check the gene ids")
  gene_ids<- feature_info$gene_id[feature_info$gene_id %in% rownames(norm_expression)]
  gene_symbols<- feature_info$gene_name[feature_info$gene_id %in% rownames(norm_expression)]
  identical(rownames(norm_expression), gene_ids)
  rownames(norm_expression)<- gene_symbols
  
  print("remove duplicated gene ids")
  express_rm <- keep_max_duplicate(norm_expression)
  saveRDS(express_rm, expression_file )
  
}

keep_max_duplicate <- function(expression_matrix){
  med <- rowMedians(expression_matrix)
  #check
  #print(med[order(-med)][1:5])
  #print(med[order(med)][1:5])
  expres_order<- expression_matrix[order(-med),]
  expres_rm<- expres_order[ !duplicated(rownames(expres_order)),]
  return(expres_rm)
}

sex_networks <-  function(processed_expression, 
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


create_alpaca_input <- function(file_path_name_female, file_path_name_male, tissue_name){
  
  #make a longer dataframe for the male network
  panda_M <- readRDS(paste0(paste0(dir_path, "/results/alpaca/sex_specific_networks/"), file_path_name_male, sep = ""))
  regNet_panda_intersection_male <- panda_M@regNet
  topNet_male <- topedges( panda_M, cutoff = 2)
  topNET_male_regnet <- as.data.frame(topNet_male@regNet)
  topNET_male_regnet$TF <- rownames(topNET_male_regnet)
  topNET_male_regnet_long <- topNET_male_regnet  %>%  pivot_longer(!TF, names_to = "gene", values_to = "male") 
  
  #make a longer dataframe for the female network
  panda_F <- readRDS(paste0(paste0(dir_path, "results/alpaca/sex_specific_networks/"),  file_path_name_female, sep = ""))
  regNet_panda_intersection_female <- panda_F@regNet
  topNet_female <- topedges( panda_F , cutoff = 2)
  topNET_female_regnet <- as.data.frame(topNet_female@regNet)
  topNET_female_regnet$TF <- rownames(topNET_female_regnet)
  topNET_female_regnet_long <- topNET_female_regnet  %>%  pivot_longer(!TF, names_to = "gene", values_to = "female")
  
  #check the tf and genes are in the same order
  print("TF")
  print(identical(topNET_female_regnet_long$TF, topNET_male_regnet_long$TF))
  print("gene")
  print(identical(topNET_female_regnet_long$gene, topNET_male_regnet_long$gene))
  #add together
  topNet_regnet_male_female_long<- cbind(topNET_male_regnet_long, topNET_female_regnet_long$female)
  
  colnames(topNet_regnet_male_female_long) <- c("TF", "gene", "male", "female")
  
  file_name<- paste0(paste0(dir_path, "/results/alpaca/alpaca_adjusted_inputs/"),  tissue_name, "_topNet_regnet_male_female_long.rds")
  
  saveRDS(topNet_regnet_male_female_long, file_name )
}


drug_ae_events<- function(index){
  #index for the case to expand the drug and adverse event pairs
  
  case_df <- expand.grid(test_drug[[index]], test_ae[[index]])
  case_df<- case_df[ !duplicated(case_df),]
  colnames(case_df)<- c("drug", "adverse_event")
  return(case_df)
}


combine_sub_data_frames<- function(index, list_samples, data_frame_list ){
  group1<-unlist(list_samples[index])
  #print(group1)
  drug_ae_df<- bind_rows(data_frame_list[group1])
  drug_ae_df1<- drug_ae_df[!duplicated(drug_ae_df),]
  return(drug_ae_df1)
}

sex_bias_adverse_event_test <- function( index, drug_event_df ){
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


alpaca_com_plot_function<- function(newfile, com_num, title){
  comm_symbols_genes <- newfile$symbol[newfile$com == com_num]
  set.seed(101)
  pathway_results_male <- gost(query = comm_symbols_genes, 
                               organism = "hsapiens", ordered_query = FALSE, 
                               multi_query = FALSE, significant = TRUE, exclude_iea= FALSE, measure_underrepresentation = FALSE, evcodes = FALSE,
                               user_threshold = 0.05, correction_method = "g_SCS", 
                               domain_scope = "annotated", custom_bg = NULL, 
                               numeric_ns = "", as_short_link = FALSE)
  pathway_results_male_comm_5_res<- pathway_results_male$result
  #print(pathway_results_male_comm_5_res$p_value[1:5])
  
  
  pathway_results_male_comm_5_res$set<- rep(title, nrow(pathway_results_male_comm_5_res))
  
  file_name1<- paste0("~/results/alpaca/alpaca_pathway_analysis/", title, "_pathways.csv")
  write.csv2(as.data.frame(pathway_results_male_comm_5_res[,-14]), file_name1 )
  
  pathway_results_male_comm_5_res<- pathway_results_male_comm_5_res[pathway_results_male_comm_5_res$source %in% c("KEGG","REAC" ),]
  pathway_results_male_comm_5_res$log_p <- -log(pathway_results_male_comm_5_res$p_value)
  
  file_name<- paste("~/results/alpaca/alpaca_pathway_analysis/gpplot_kegg_reactome_", title, ".png", sep = "")
  ggplot(pathway_results_male_comm_5_res, aes(x=term_name, y=log_p))+ geom_point(aes(size = recall)) +coord_flip() + ggtitle(title) + scale_color_viridis(option="plasma")+ theme_bw() +theme(text = element_text(size=10)) + scale_size(range = c(5,10))
  ggsave( file_name)
}

alpaca_output_dw_tf_genes <- function(scores, comms){
  #scores<- male_female_ALPACA_scores
  #comms<- male_female_ALPACA_final_memb
  tosel <- intersect(scores[,1], comms[,1])
  scores <- scores[which(scores[,1] %in% tosel),]
  comms <- comms[which(comms[,1] %in% tosel),]
  scores <- scores[order(scores[,1]),]
  comms <- comms[order(comms[,1]),]
  scores[,1] <- as.character(scores[,1])
  comms[,1] <- as.character(comms[,1])
  all(scores[,1]==comms[,1])
  scores <- cbind(scores, comms[,2])
  row.names(scores) <- scores[,1]
  scores <- scores[,-1]
  colnames(scores) <- c("score", "com")
  alpaca_nodes <- row.names(scores)
  table(substr(alpaca_nodes, nchar(alpaca_nodes), nchar(alpaca_nodes))) # 714 A (tf) 22572 B (gene) in test network
  alpaca_genes <- scores
  alpaca_genes <- alpaca_genes[order(alpaca_genes[,2],decreasing=T),]
  alpaca_genes <- alpaca_genes[,c(2,1)]
  newfile<- matrix(nrow= 1, ncol=2)
  colnames(newfile)<- colnames(alpaca_genes)
  for(n in 1:max(alpaca_genes$com)){
    junk <- alpaca_genes[which(alpaca_genes$com==n),]
    junk <- junk[order(junk$score,decreasing=T),]
    if(nrow(junk)>100){
      junk <- junk[1:100,]
      if(n==1){
        newfile <- junk
      } else {
        newfile <- rbind(newfile, junk)
      }
    }
  }   
  newfile$symbol <- substr(row.names(newfile), 1, nchar(row.names(newfile))-2)
  newfile<- newfile[-1,]
  return(newfile)
}

alpaca_pathway_wrapper<- function(tissue_name, sex){
  print(tissue_name)
  print(sex)
  if (sex == "female"){
    scores_file <- paste0("~/results/alpaca/alpaca_run/", tissue_name, "__ALPACA_scores.txt")
    comms_file <- paste0("~/results/alpaca/alpaca_run/", tissue_name, "__ALPACA_final_memb.txt")
  }else{
    scores_file <- paste0("~/results/alpaca/alpaca_run/", tissue_name, "femalemale__ALPACA_scores.txt")
    comms_file <- paste0("~/results/alpaca/alpaca_run/", tissue_name, "femalemale__ALPACA_final_memb.txt")
  }
  #read in the score file 
  #read in the final memb file 
  file_scores <- read.delim(scores_file, header=FALSE)
  file_comms<- read.delim(comms_file, header=FALSE)
  
  # top 100 genes in the community
  alpaca_core_genes <- alpaca_output_dw_tf_genes(file_scores, file_comms)
  #save the tf and gene files for each community
  core_name<- paste(tissue_name, sex, "core_genes.rds", sep= "_")
  file_core_name<- paste0("~/results/alpaca/alpaca_core_gene_lists/",core_name)
  saveRDS(alpaca_core_genes,file_core_name)
  
  #pathway analysis for each community with 100 genes 
  num_comm<- unique(alpaca_core_genes[,1])
  print(num_comm)
  for( i in 1:length(num_comm)){
    index<- num_comm[i]
    test_name<- paste(tissue_name, sex, index, "community",  sep= "_")
    print(test_name)
    try(alpaca_com_plot_function(alpaca_core_genes, index, test_name))
  }
}

#functions from rrvgo
#rrvgo github (getGoTerm, loadOrgdb, getGoSize, reduceSimMatrix)
getGoTerm <- function(x) {
  sapply(x, function(x) tryCatch(GO.db::GOTERM[[x]]@Term, error=function(e) NA))
}
loadOrgdb <- function(orgdb) {
  if(!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Bioconductor orgdb for ", orgdb, " not found. Consider installing it.",
         call. = FALSE)
  }
  eval(parse(text=paste0(orgdb, "::", orgdb)))
}
getGoSize <- function(terms, orgdb, keytype) {
  if(all(is(orgdb) != "OrgDb")) {
    orgdb <- loadOrgdb(orgdb)
  }
  
  # get all GO terms with genes associated
  go <- suppressMessages(
    AnnotationDbi::select(orgdb,
                          keytype=keytype,
                          columns=c("GO", "ONTOLOGY"),
                          keys=AnnotationDbi::keys(orgdb, keytype=keytype)))
  go <- go[!is.na(go$GO), ]
  go <- go[go$GO %in% terms, ]
  
  # count
  counts   <- table(go$GO)
  go <- go[go$GO %in% terms, ]
  empty    <- terms[!(terms %in% names(counts))]
  nocounts <- setNames(rep(0, length(empty)), empty)
  
  c(counts, nocounts)
}
reduceSimMatrix <- function (simMatrix, scores = NULL, threshold = 0.7, orgdb="org.Hs.eg.db", keytype = "ENTREZID") {
  if (!is.null(scores) && !all(rownames(simMatrix) %in% names(scores))) {
    stop("Scores vector does not contain all terms in the similarity matrix")
  }
  sizes <- getGoSize(rownames(simMatrix), orgdb, keytype)
  if (is.null(scores)) {
    message("No scores provided. Falling back to term's size")
    scores <- sizes
  }
  scores <- scores[match(rownames(simMatrix), names(scores))]
  orows <- match(rownames(simMatrix), names(scores))
  ocols <- match(colnames(simMatrix), names(scores))
  simMatrix <- simMatrix[orows, ocols]
  o <- rev(order(scores, sizes, na.last = FALSE))
  simMatrix <- simMatrix[o, o]
  cluster <- cutree(hclust(as.dist(1 - simMatrix)), h = threshold)
  clusterRep <- tapply(rownames(simMatrix), cluster, function(x) x[which.max(scores[x])])
  data.frame(go = rownames(simMatrix), cluster = cluster, parent = clusterRep[cluster], 
             parentSimScore = unlist(Map(seq_len(nrow(simMatrix)), 
                                         clusterRep[cluster], f = function(i, j) simMatrix[i, 
                                                                                           j])), score = scores[match(rownames(simMatrix), 
                                                                                                                      names(scores))], size = sizes[match(rownames(simMatrix), 
                                                                                                                                                          names(sizes))], term =getGoTerm(rownames(simMatrix)), 
             parentTerm = getGoTerm(clusterRep[cluster]))
}

alpaca_output_dw_all <- function(scores, comms){
  #scores<- male_female_ALPACA_scores
  #comms<- male_female_ALPACA_final_memb
  tosel <- intersect(scores[,1], comms[,1])
  scores <- scores[which(scores[,1] %in% tosel),]
  comms <- comms[which(comms[,1] %in% tosel),]
  scores <- scores[order(scores[,1]),]
  comms <- comms[order(comms[,1]),]
  scores[,1] <- as.character(scores[,1])
  comms[,1] <- as.character(comms[,1])
  all(scores[,1]==comms[,1])
  scores <- cbind(scores, comms[,2])
  row.names(scores) <- scores[,1]
  scores <- scores[,-1]
  colnames(scores) <- c("score", "com")
  alpaca_nodes <- row.names(scores)
  table(substr(alpaca_nodes, nchar(alpaca_nodes), nchar(alpaca_nodes))) # 714 A (tf) 22572 B (gene) in test network
  alpaca_genes <- scores
  alpaca_genes <- alpaca_genes[order(alpaca_genes[,2],decreasing=T),]
  alpaca_genes <- alpaca_genes[,c(2,1)]
  
  alpaca_genes$symbol <- substr(row.names(alpaca_genes), 1, nchar(row.names(alpaca_genes))-2)
  return(alpaca_genes)
}

alpaca_gene_set_analysis <- function(tissue_nam, genes="drugs", gene_set_name= "Drug Metabolism Genes", path="~/results/alpaca/drug_metabolism_genes/", save_com_list= TRUE){
  print(tissue_nam)
  if(genes== "drugs"){
    print("using drug gene list from KEGG")
    genes <- c("ADH1A"	,"ADH1B", "ADH1C", "ADH4",	"ADH5",	"ADH6",	"ADH7",	
               "ALDH1A3",	"ALDH3A1",	"ALDH3B1", "ALDH3B2", "AOX1", "CYP1A2",
               "CYP2A13", "CYP2A6", "CYP2A7",	"CYP2B6",	"CYP2C18",	"CYP2C19",
               "CYP2C8",	"CYP2C9",	"CYP2D6",	"CYP2E1",	"CYP3A4",	"CYP3A43",	"CYP3A5",
               "CYP3A7",	"FMO1",	"FMO2",	"FMO3", "FMO4",	"FMO5",  "GSTA1", "GSTA2", 	"GSTA3",	"GSTA4",
               "GSTA5",	"GSTK1",	"GSTM1",	"GSTM2",	"GSTM3",	"GSTM4",	"GSTM5",	"GSTO1",
               "GSTO2",	"GSTP1",	"GSTT1",	"GSTT2",	"GSTZ1",	"MAOA",	"MAOB",	"MGST1",
               "MGST2",	"MGST3",	"UGT1A1",	"UGT1A10",	"UGT1A3",	"UGT1A4",
               "UGT1A5",	"UGT1A6",	"UGT1A7",	"UGT1A8", "UGT1A9",	"UGT2A1",	"UGT2A3", "UGT2B10",
               "UGT2B11",	"UGT2B15",	"UGT2B17",	"UGT2B28",	"UGT2B4",	"UGT2B7")
  }
  
  #get the data
  scores_file_f <- paste0("~/results/alpaca/alpaca_run/", tissue_nam, "__ALPACA_scores.txt")
  comms_file_f <- paste0("~/results/alpaca/alpaca_run/", tissue_nam, "__ALPACA_final_memb.txt")
  file_scores_f <- read.delim(scores_file_f, header=FALSE)
  file_comms_f <- read.delim(comms_file_f, header=FALSE)
  test<- alpaca_output_dw_all(file_scores_f, file_comms_f)
  
  
  scores_file_m <- paste0("~/results/alpaca/alpaca_run/", tissue_nam, "femalemale__ALPACA_scores.txt")
  comms_file_m <- paste0("~/results/alpaca/alpaca_run/", tissue_nam, "femalemale__ALPACA_final_memb.txt")
  file_scores_m <- read.delim(scores_file_m, header=FALSE)
  file_comms_m <- read.delim(comms_file_m, header=FALSE)
  test2<- alpaca_output_dw_all(file_scores_m, file_comms_m)
  
  if (save_com_list){
    #save this data- alpaca_com_gene_list
    file_name<- paste("~/results/alpaca/alpaca_com_gene_list/", tissue_nam, "_female_alpaca_com_all_gene_df.rds", sep = "")
    saveRDS(test, file_name)
    #save this data- alpaca_com_gene_list
    file_name<- paste("~/results/alpaca/alpaca_com_gene_list/", tissue_nam, "_male_alpaca_com_all_gene_df.rds", sep = "")
    saveRDS(test2, file_name)
  }
  
  
  y_name<- paste("Number of ", gene_set_name)
  #female plot 
  drug_subset <- test[test$symbol %in% genes,]
  drug_subset_list <- factor(drug_subset$com, levels=1:length(unique(test$com)))
  drug_gene_count<- as.data.frame(table(drug_subset_list))
  drug_gene_count$drug_subset_list <-  paste0( "female", drug_gene_count$drug_subset_list)
  
  ggplot(drug_gene_count, aes(x= drug_subset_list, y= Freq, label=Freq)) + geom_bar(stat="identity", fill = "#440154FF", color= "black", alpha=0.7)+
    geom_text(size = 10, position = position_stack(vjust = 0.9)) + xlab("Communities") + ylab(y_name)+theme(text = element_text(size = 30,  face="bold"))
  #save this plot - count_bar_plots
  file_name<- paste(path, "count_bar_plots/", tissue_nam, "_female_alpaca_com_gene_set.png", sep = "")
  ggsave(file_name, width = 20, height=10, units= "in")
  
  #male plot
  drug_subset2 <- test2[test2$symbol %in% drug_genes,]
  drug_subset_list2 <- factor(drug_subset2$com, levels=1:length(unique(test2$com)))
  drug_gene_count2<- as.data.frame(table(drug_subset_list2))
  drug_gene_count2$drug_subset_list <-  paste0( "male", drug_gene_count2$drug_subset_list2)
  
  ggplot(drug_gene_count2, aes(x= drug_subset_list, y= Freq, label=Freq)) + geom_bar(stat="identity", fill = "#21908CFF", color= "black", alpha=0.7)+ 
    geom_text(size = 10, position = position_stack(vjust = 0.9)) + xlab("Communities") + ylab(y_name)+theme(text = element_text(size = 30,  face="bold"))
  #save this plot  - count_bar_plots
  file_name<- paste(path, "count_bar_plots/", tissue_nam, "_male_alpaca_com_gene_set.png", sep = "")
  ggsave(file_name, width = 20, height=10, units= "in")
  
  
  
  drug_subset2$sex<- rep("male", nrow(drug_subset2))
  drug_subset$sex<- rep("female", nrow(drug_subset))
  
  liver_drug_subset<- rbind(drug_subset, drug_subset2)
  liver_drug_subset<- liver_drug_subset[order(liver_drug_subset$symbol),]
  #save the liver_drug_subset- drug_gene_com
  file_name<- paste(path, "gene_com/", tissue_nam, "_alpaca_com_gene_set_df.rds", sep = "")
  saveRDS(liver_drug_subset, file_name)
  
  #sometimes there is not duplicates- removing genes without duplicate
  keep_genes<- names(table(liver_drug_subset$symbol))[table(liver_drug_subset$symbol) >1]
  liver_drug_subset<- liver_drug_subset[liver_drug_subset$symbol %in% keep_genes,]
  
  liver_drug_subset_female<- liver_drug_subset[liver_drug_subset$sex == "female",]
  liver_drug_subset_male<- liver_drug_subset[liver_drug_subset$sex == "male",]
  res <- wilcox.test(liver_drug_subset_female$score, liver_drug_subset_male$score, paired = TRUE)
  p_value <- res$p.value
  
  diff<- log2(mean(liver_drug_subset_female$score)) - log2(mean(liver_drug_subset_male$score))
  
  #adjust the names here
  y_name2<- paste("number of ", gene_set_name)
  title_tissue<- paste("Differential modularity score of ", y_name2, "\nin", tissue_nam, "sex-specific gene regulatory networks", sep = " ")
  label_tissue <- paste("Paired Wilcoxon signed rank test\np-value =", p_value, sep = " ")
  
  ggplot(liver_drug_subset, aes(y=score, x= sex, fill = sex)) + geom_violin() +geom_point() + geom_text(x="female", y=0.002, label= label_tissue, size= 15) +ylab("Differential modularity score") + xlab("Sex") + ggtitle(title_tissue) + scale_fill_manual(values= c( "#440154FF","#21908CFF")) +theme(text = element_text(size = 30,  face="bold"))
  #save the plot - sex_differential_modularity_plots
  file_name<- paste(path, "sex_differential_modularity_plots/", tissue_nam, "_sex_differential_modularity_plot.png", sep = "")
  ggsave(file_name, width = 20, height=10, units= "in")
  
  #return p_value 
  obj<- list(p_value, diff)
  return(obj)
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

gene_com_list<- function(tissue, sex){
  com_df <- readRDS(paste("~/results/alpaca/alpaca_com_gene_list/", tissue, "_", sex,  "_alpaca_com_all_gene_df.rds", sep =""))
  com_numbers<- unique(com_df$com)
  com_sub_list<- list()
  for(i in 1:length(com_numbers)){
    com_sub_list[[i]]<- com_df$symbol[com_df$com ==  com_numbers[i]]
  }
  names(com_sub_list)<- paste0(tissue, "@" ,sex, "@", com_numbers)
  return(com_sub_list)
}

gene_core_list<- function(tissue, sex){
  com_df <- readRDS(paste("~/results/alpaca/alpaca_core_gene_lists/", tissue, "_", sex,  "_core_genes.rds", sep =""))
  com_numbers<- unique(com_df$com)
  com_sub_list<- list()
  for(i in 1:length(com_numbers)){
    com_sub_list[[i]]<- com_df$symbol[com_df$com ==  com_numbers[i]]
  }
  names(com_sub_list)<- paste0(tissue, "@" ,sex, "@", com_numbers)
  return(com_sub_list)
}

go_term_heatmap<- function(female_pathways, male_pathways, threshold = 0.95, file_name){
  #input
  #female_pathways, female_pathways- the gprofiler result data frame for each method with the last column (named set)
  #threshold- is the threshold of how far up or down the go ontology tree to go. default to 0.95
  #file name for the parent go term results
  
  #output 
  #heatmap
  
  #get terms
  female_bp<- female_pathways$term_id[ female_pathways$source == "GO:BP" ]
  male_bp <- male_pathways$term_id[ male_pathways$source == "GO:BP"]
  
  #run the go term semantic similarity 
  go1 <-  unique(c(female_bp, male_bp ))
  go_gbm_up_sim <- mgoSim(go1, go1, semData=hsGO, measure="Wang", combine=NULL)
  
  #get the parent terms from rrvgo for the pathways
  res <- reduceSimMatrix(go_gbm_up_sim, threshold = threshold)
  res_v2<- res[match(colnames(go_gbm_up_sim), res$go),]
  
  
  
  #determine which pathway is enriched in the different methods for the heatmap annotation
  female_drugs <- grepl(paste(female_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  male_drugs<- grepl(paste(male_bp,collapse="|"), colnames(go_gbm_up_sim))
  
  #handing the case where there are no pathways enriched
  if(length(female_bp) ==0){ female_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  if(length(male_bp) ==0){ male_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  
  #pick the number of colors note limit is about 9 
  bp_color<- brewer.pal(n = length(unique(res_v2$parentTerm)), name = "Paired")
  names(bp_color)<- unique(res_v2$parentTerm)
  
  #create heatmap
  row_ha <- HeatmapAnnotation(Female=female_drugs, Male= male_drugs, GO_BP_Group= res_v2$parentTerm , col = list(Female = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"), Male = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"), GO_BP_Group= bp_color ), annotation_name_gp= gpar(fontsize = 16,  fontface = "bold"), annotation_legend_param = list(Female= list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")),  Male = list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")),  GO_BP_Group = list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold"))))
  
  col_fun <- colorRamp2(c(0,  1), c( "black", "yellow"))
  
  #SAVE THE PARENT MATCHES
  res_v3<- cbind(res_v2,female_drugs, male_drugs)
  #return(res_v3)
  write.csv(res_v3, file_name)
  heatmap <- Heatmap(go_gbm_up_sim, nam= "GO Term Similarity (Wang)", col = col_fun, show_column_names = FALSE,  show_row_names = FALSE, top_annotation = row_ha, clustering_distance_rows= "euclidean", clustering_distance_columns=  "euclidean", clustering_method_rows = "ward.D2" ,  clustering_method_columns="ward.D2", heatmap_legend_param = list( title_gp = gpar(fontsize = 15 , fontface = "bold"), labels_gp = gpar(fontsize = 12, fontface = "bold")))
  
  draw(heatmap, annotation_legend_side = "bottom")
  
  # 
  # column_dend = hclust(dist(t(go_gbm_up_sim)), method = "ward.D2")
  # 
  # heatmap_v2 <- Heatmap(matrix(nrow = 0, ncol = ncol(go_gbm_up_sim)),
  #         cluster_columns = column_dend, 
  #         show_column_names = FALSE,  
  #         show_row_names = FALSE, 
  #         top_annotation = row_ha) 
  # draw(heatmap_v2, annotation_legend_side = "bottom")
  # 
}

drug_target_diff_analysis<- function(tissue, path){
  #input 
  #tissue- name of the tissues (match file names)
  #path- directory path to sex specific networks
  
  #outputs 
  # plot of the degree of drug metabolism genes vs other genes
  # p_value- wilcox of degree difference between drug metabolism genes and other genes
  # diff- median of the degree difference of drug metabolism genes (female-male)
  male_file<- paste0("~/results/alpaca/sex_specific_networks/", tissue, "_M_panda.rds")
  M_panda <- readRDS(male_file)
  female_file<- paste0("~/results/alpaca/sex_specific_networks/", tissue, "_F_panda.rds")
  F_panda <- readRDS(female_file)
  
  diff_genes <- calcDegreeDifference( F_panda , M_panda, type="gene", filter = TRUE)
  
  
  
  diff_drug_genes <- diff_genes[names(diff_genes) %in% drug_genes]
  diff_other_genes <- diff_genes[!names(diff_genes) %in% drug_genes]
  
  diff_all_genes<- rbind(data.frame(gene_list =rep( "Drug\nMetabolism", length(diff_drug_genes)), Gene_Degree_Diff= as.numeric(diff_drug_genes)), data.frame(gene_list =rep( "Other", length(diff_other_genes)), Gene_Degree_Diff= as.numeric(diff_other_genes)))
  
  res <- wilcox.test(diff_drug_genes , diff_other_genes  )
  label_tissue <- paste("Wilcoxon signed\nrank test\np-value =\n",  res$p.value, sep = " ")
  title_tissue<- paste("Degree difference of drug metabolismgenes\nbetween", tissue, "sex-specific gene regulatory networks", sep = " ")
  diff_all_genes$gene_list<- factor(diff_all_genes$gene_list, levels= c("Other","Drug\nMetabolism" ))
  
  
  ggplot(diff_all_genes, aes(x=Gene_Degree_Diff, y= gene_list, fill = gene_list)) + geom_violin() +geom_point() +geom_text(y="Drug\nMetabolism", x=-250, label= label_tissue, size= 10) +ggtitle(title_tissue) +xlab("Degree Difference (FEMALE-MALE)") + ylab("Genes")  + scale_fill_viridis( alpha=0.7, option= "H" ,discrete = TRUE)+theme(text = element_text(size = 30,  face="bold"))+ theme(legend.position = "none")
  
  file_name<- paste(path, "sex_degree_difference_plots/", tissue, "_sex_degree_difference_plot.png", sep = "")
  ggsave(file_name, width = 20, height=10, units= "in")
  
  p_value<- res$p.value
  diff<- median(diff_drug_genes)
  #return p_value 
  obj<- list(p_value, diff)
  return(obj)
}

calc_targeting<- function(tissue, path){
  #input 
  #tissue- name of the tissues (match file names)
  #path- directory path to sex specific networks
  male_file<- paste0("~/results/alpaca/sex_specific_networks/", tissue, "_M_panda.rds")
  M_panda <- readRDS(male_file)
  female_file<- paste0("~/results/alpaca/sex_specific_networks/", tissue, "_F_panda.rds")
  F_panda <- readRDS(female_file)
  
  f_degree_genes <- calcDegree( F_panda , type="gene", filter = TRUE, trim=TRUE)
  m_degree_genes <- calcDegree( M_panda, type="gene", filter = TRUE, trim=TRUE)
  
  #save degree in here
  
  #panda_in_degree
  file<- paste(path, "panda_in_degree/", tissue, "_female_sex_degree_gene.rds", sep = "")
  saveRDS(f_degree_genes, file)
  file<- paste(path, "panda_in_degree/", tissue, "_male_sex_degree_gene.rds", sep = "")
  saveRDS(m_degree_genes, file)
  all_degree <-  f_degree_genes  + m_degree_genes
  
  f_prop_degree <- f_degree_genes / all_degree
  file<- paste(path, "panda_in_degree/", tissue, "_female_sex_degree_gene_prop.rds", sep = "")
  saveRDS(f_prop_degree, file)
  
  m_prop_degree <- m_degree_genes / all_degree
  file<- paste(path, "panda_in_degree/", tissue, "_male_sex_degree_gene_prop.rds", sep = "")
  saveRDS(m_prop_degree, file)
  
  #find sex divergent here
  #the proportion of sex-biased edges in the male- and female-biased directions is between 0.4 and 0.6
  sex_div <- names(f_prop_degree)[ f_prop_degree > 0.4 &  f_prop_degree< 0.6]
  
  #find female-bias here 
  f_bias<- names(f_prop_degree)[ f_prop_degree >0.6]
  #the proportion of sex-biased edges in the female direction is greater than 0.6
  
  #find male-bias here
  #male-biased genes (the proportion of sex-biased edges in the male direction is greater than 0.6),
  m_bias<- names(m_prop_degree)[ m_prop_degree >0.6]
  #create list 
  list_bias<- list(sex_div, f_bias, m_bias)
  
  #save 
  file<- paste(path, "panda_in_degree/", tissue, "_sex_bias_gene_list.rds", sep = "")
  saveRDS(list_bias, file)
  
  
  #out- degree
  f_degree_tf <- calcDegree( F_panda , type="tf", filter = TRUE, trim=TRUE)
  m_degree_tf <- calcDegree( M_panda, type="tf", filter = TRUE, trim=TRUE)
  
  #save degree in here
  file<- paste(path, "panda_out_degree/", tissue, "_female_sex_degree_tf.rds", sep = "")
  saveRDS(f_degree_tf, file)
  file<- paste(path, "panda_out_degree/", tissue, "_male_sex_degree_tf.rds", sep = "")
  saveRDS(m_degree_tf, file)
  all_degree <-  f_degree_tf  + m_degree_tf
  
  f_prop_degree <- f_degree_tf / all_degree
  file<- paste(path, "panda_out_degree/", tissue, "_female_sex_degree_tf_prop.rds", sep = "")
  saveRDS(f_prop_degree, file)
  
  m_prop_degree <- m_degree_tf / all_degree
  file<- paste(path, "panda_out_degree/", tissue, "_male_sex_degree_tf_prop.rds", sep = "")
  saveRDS(m_prop_degree, file)
  
  #find sex divergent here
  #the proportion of sex-biased edges in the male- and female-biased directions is between 0.4 and 0.6
  sex_div <- names(f_prop_degree)[ f_prop_degree > 0.4 &  f_prop_degree< 0.6]
  
  #find female-bias here 
  f_bias<- names(f_prop_degree)[ f_prop_degree >0.6]
  #the proportion of sex-biased edges in the female direction is greater than 0.6
  
  #find male-bias here
  #male-biased tf (the proportion of sex-biased edges in the male direction is greater than 0.6),
  m_bias<- names(m_prop_degree)[ m_prop_degree >0.6]
  #create list 
  list_bias<- list(sex_div, f_bias, m_bias)
  
  #save 
  file<- paste(path, "panda_out_degree/", tissue, "_sex_bias_tf_list.rds", sep = "")
  saveRDS(list_bias, file)
  
}  

go_term_heatmap_v2<- function(female_act_pathways, female_repress_pathways, male_act_pathways, male_repress_pathways, threshold = 0.95, file_name){
  #input
  #limma_pathways, deseqq2_pathways, TFL_pathways- the gprofiler result data frame for each method with the last column (named set)
  #indicating if the pathway was from the up or down group
  #list_type- up or down regulated 
  #threshold- is the threshold of how far up or down the go ontology tree to go. default to 0.95
  #file name for the parent go term results
  
  #output 
  #heatmap
  
  #get terms
  female_act_bp<- female_act_pathways$term_id[ female_act_pathways$source == "GO:BP" ]
  female_repress_bp <- female_repress_pathways$term_id[ female_repress_pathways$source == "GO:BP"]
  male_act_bp<- male_act_pathways$term_id[ male_act_pathways$source == "GO:BP" ]
  male_repress_bp <- male_repress_pathways$term_id[ male_repress_pathways$source == "GO:BP"]
  
  #run the go term semantic similarity 
  go1 <-  unique(c(female_act_bp, female_repress_bp, male_act_bp, male_repress_bp ))
  go_sim <- mgoSim(go1, go1, semData=hsGO, measure="Wang", combine=NULL)
  
  #get the parent terms from rrvgo for the pathways
  res <- reduceSimMatrix(go_sim, threshold = threshold)
  res_v2<- res[match(colnames(go_sim), res$go),]
  
  
  
  #determine which pathway is enriched in the different methods for the heatmap annotation
  female_act_list <- grepl(paste(female_act_bp,collapse="|"), colnames(go_sim)) 
  female_repress_list <- grepl(paste(female_repress_bp,collapse="|"), colnames(go_sim)) 
  male_act_list <- grepl(paste(male_act_bp,collapse="|"), colnames(go_sim)) 
  male_repress_list <- grepl(paste(male_repress_bp,collapse="|"), colnames(go_sim)) 
  
  
  #handing the case where there are no pathways enriched
  if(length(female_act_list) ==0){ female_act_list <- rep(FALSE,ncol(go_sim) )} 
  if(length(female_repress_list) ==0){ female_repress_list <- rep(FALSE,ncol(go_sim) )} 
  if(length( male_act_list) ==0){  male_act_list <- rep(FALSE,ncol(go_sim) )} 
  if(length(male_repress_list) ==0){ male_repress_list <- rep(FALSE,ncol(go_sim) )} 
  
  
  #pick the number of colors note limit is about 9 
  bp_color<- brewer.pal(n = length(unique(res_v2$parentTerm)), name = "Paired")
  names(bp_color)<- unique(res_v2$parentTerm)
  
  #create heatmap
  #"Female Repressor Edges"=f_repress_edge, "Male Activator Edges"=m_act_edge, "Male Repressor Edges"=m_repress_edge)
  row_ha = HeatmapAnnotation("Female Activator Edges"=female_act_list, "Female Repressor Edges"= female_repress_list, "Male Activator Edges"=  male_act_list, "Male Repressor Edges" = male_repress_list, "Common Parent GO Term"= res_v2$parentTerm , col = list("Female Activator Edges"= c("TRUE" = "black", "FALSE" = "white"),"Female Repressor Edges" = c("TRUE" = "black", "FALSE" = "white"),"Male Activator Edges" = c("TRUE" = "black", "FALSE" = "white"),"Male Repressor Edges" = c("TRUE" = "black", "FALSE" = "white"),
                                                                                                                                                                                                                                                                  "Common Parent GO Term"= bp_color ), annotation_name_gp= gpar(fontsize = 16,  fontface = "bold"), 
                             annotation_legend_param = list("Female Repressor Edges"= list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")), 
                                                            "Female Activator Edges" = list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")), 
                                                            "Male Activator Edges"= list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")),"Male Repressor Edges"= list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")),
                                                            "Common Parent GO Term" = list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold"))))
  col_fun = colorRamp2(c(0,  1), c( "black", "yellow"))
  
  #SAVE THE PARENT MATCHES
  res_v3<- cbind(res_v2,female_act_list, female_repress_list, male_act_list, male_repress_list )
  write.csv(res_v3, file_name)
  heatmap_v1 <- Heatmap(go_sim, nam= "GO Term Similarity (Wang)", col = col_fun, show_column_names = FALSE,  show_row_names = FALSE, top_annotation = row_ha,  
                        clustering_distance_rows= "euclidean",
                        clustering_distance_columns=  "euclidean",
                        clustering_method_rows = "ward.D2" ,
                        clustering_method_columns="ward.D2", heatmap_legend_param = list( title_gp = gpar(fontsize = 15 , fontface = "bold"), labels_gp = gpar(fontsize = 12, fontface = "bold")))
  #draw(heatmap_v1, annotation_legend_side = "bottom")
  #heatmap
  Heatmap(go_sim, nam= "GO Term Similarity (Wang)", col = col_fun, show_column_names = FALSE,  show_row_names = FALSE, top_annotation = row_ha,  
          clustering_distance_rows= "euclidean",
          clustering_distance_columns=  "euclidean",
          clustering_method_rows = "ward.D2" ,
          clustering_method_columns="ward.D2", heatmap_legend_param = list( title_gp = gpar(fontsize = 15 , fontface = "bold"), labels_gp = gpar(fontsize = 12, fontface = "bold")))
  
  column_dend = hclust(dist(t(go_sim)), method = "ward.D2")
  
  heatmap_v2 <- Heatmap(matrix(nrow = 0, ncol = ncol(go_sim)),
                        cluster_columns = column_dend, 
                        show_column_names = FALSE,  
                        show_row_names = FALSE, 
                        top_annotation = row_ha) 
  draw(heatmap_v2, annotation_legend_side = "bottom")
  
}

#old functions no longer used in main project
lioness_output_adjustment <- function(index){
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

panda_tissue_run <- function(norm_expression,
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

