FROM bioconductor/bioconductor_docker:RELEASE_3_16

# location for mounting
#RUN mkdir /home/rstudio/data


RUN R -e 'BiocManager::install(c("DESeq2", "pasilla", "SummarizedExperiment", "TFEA.ChIP", "ensembldb", "EnsDb.Hsapiens.v75", "AnnotationDbi", "org.Hs.eg.db", "cogena", "BgeeDB", "projectR", "GEOquery", "recount3", "ComplexHeatmap", "factoextra", "limma", "ChemmineR", "customCMPdb", "fmcsR", "PharmacoGx" ))'

RUN R -e 'install.packages(c("readr","tidyverse","ggforce", "plotly", "devtools", "ProliferativeIndex", "gprofiler2","gplots","RColorBrewer", "car", "VennDiagram", "DrInsight", "FactoMineR", "viridis", "ComplexUpset","tidymodels", "caret", "stats", "MASS", "gbm", "rpart", "RWeka", "LogitBoost", "pamr", "naivebayes", "knn", "klaR",  "xgboost", "glm", "kernlab",  "LiblineaR", "discrim", "ranger", "randomForest", "partykit" , "nnet", "earth", "igraph", "mashr", "broom.mixed", "dotwhisker", "skimr", "visNetwork" ))'

RUN R -e 'devtools::install_github("wgmao/PLIER")'
RUN R -e 'devtools::install_github("yduan004/signatureSearch", build_vignettes=FALSE)'
RUN R -e 'BiocManager::install("signatureSearchData")'
RUN R -e 'devtools::install_github("netZoo/netZooR")' 
RUN R -e 'devtools::install_github("yduan004/drugbankR")'
RUN R -e 'utils::install.packages("markdown", repos = "https://packagemanager.rstudio.com/cran/__linux__/jammy/latest")'