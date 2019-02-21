#Analysis of microglial data in Fenna/Matt dataset 
library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)
library(moments) 
library(NNLM)
library(ggplot2)
source("SeuratExtrafunctions.R")
library(data.table)
library(broom)
library(tidyr)
library(devtools)
#install_github('MacoskoLab/liger')
library(liger)

#Set paths and load processed cells 
path<-"/Users/michaeljohndolan/Google Drive (mdolan@broadinstitute.org)/"
setwd(path)

#Load up the data 
#hammond<-readRDS(file = "Hammond2018_microglia_DGE/Round_2_40.filtered.scaled.dge.RDS")
hammond<-CreateSeuratObject(hammond)
mgls<-readRDS("FennaMatt_dLGN_scRNAseq/microglia_Reprocessed.rds")
hammond<-readRDS(file = "FennaMatt_dLGN_scRNAseq/HammondP4P5_processed.rds")

#Extract P4/P5 cells from Hammond dataset 
hammond<-SetAllIdent(hammond, id = "orig.ident") 
hammond.P4.P5<-SubsetData(object = hammond, do.clean = TRUE, ident.use = c("P4", "P5"))

#Process Tim's P4/P5 dataset 
mito.genes<-grep(pattern = "^mt.", x = rownames(x = hammond.P4.P5@raw.data), value = TRUE)
percent.mito<-Matrix::colSums(hammond.P4.P5@raw.data[mito.genes, ])/Matrix::colSums(hammond.P4.P5@raw.data)
hammond.P4.P5<-AddMetaData(object = hammond.P4.P5, metadata = percent.mito, col.name = "percent.mito")
hammond.P4.P5<-FilterCells(object = hammond.P4.P5, subset.names = c("nGene", "nUMI","percent.mito"), 
                     low.thresholds = c(200, 500, -Inf), high.thresholds = c(2500, 15000, 0.1))
hammond.P4.P5<-NormalizeData(object = hammond.P4.P5, normalization.method = "LogNormalize", scale.factor = 10000)
hammond.P4.P5<-FindVariableGenes(object = hammond.P4.P5, mean.function = ExpMean, dispersion.function = LogVMR, 
                          do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
hammond.P4.P5<-ScaleData(object = hammond.P4.P5, genes.use = hammond.P4.P5@var.genes, display.progress = TRUE, 
                    vars.to.regress = c("percent.mito",'nUMI'), do.par = TRUE, num.cores = 4)
hammond.P4.P5<-RunPCA(object = hammond.P4.P5, pc.genes = hammond.P4.P5@var.genes,pcs.compute = 30 , do.print = TRUE, pcs.print = 1:5, 
                genes.print = 5)
PCElbowPlot(object = hammond.P4.P5)
hammond.P4.P5<-FindClusters(object = hammond.P4.P5, reduction.type = "pca", dims.use = 1:15, 
                       resolution = 0.6, print.output = 0, save.SNN = TRUE)
hammond.P4.P5<-RunTSNE(object = hammond.P4.P5, dims.use = 1:15, do.fast = TRUE)
TSNEPlot(object = hammond.P4.P5, pt.size = 0.05, do.label = TRUE)
saveRDS(hammond.P4.P5, "FennaMatt_dLGN_scRNAseq/HammondP4P5_processed.rds")

#Run liger to combine the two datasets 
ligerex<-createLiger(list(dLGN, Hammond2018), combined.seurat = F, use.tsne = F)

#Normalize, scale (but not center) and find variable genes in the shared dataset 
ligerex = normalize(ligerex)
#ligerex = selectGenes(ligerex, var.thresh = 0.1) Is this necessary if we did it and take union?
ligerex = scaleNotCenter(ligerex)

#Perform the factorization
ligerex = optimizeALS(ligerex, k = 20) 
ligerex = quantileAlignSNF(ligerex) #SNF clustering and quantile alignment

#Visualize the alignment
ligerex = runTSNE(ligerex)
plotByDatasetAndCluster(ligerex) #Can also pass in different set of cluster labels to plot

