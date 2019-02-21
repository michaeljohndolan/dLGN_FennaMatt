#Examine cluster 8 from Hammond (2018) accross multiple different timepoints. Will use this data to parse the Liger 
#alignment of Fenna and Hammond datasets. 

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)

#Initial df processing to extract all timepoint,sex and genotype data
Tim.Mgl.assign<-readRDS("Google Drive (mdolan@broadinstitute.org)/Hammond2018_microglia_DGE/Round_2_40.cluster.assign.RDS")
Tim.Mgl.assign<-as.data.frame(Tim.Mgl.assign)
Tim.Mgl.assign$Barcodes<-row.names(Tim.Mgl.assign)
row.names(Tim.Mgl.assign) <- NULL
strsplit(Tim.Mgl.assign$Barcodes, "_")