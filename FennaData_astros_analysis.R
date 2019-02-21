#Analysis of astrocyte data in Fenna/Matt dataset 
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
gene.Sum<-function(obj, norm.to.cell=F) {
  obj<-as.matrix(obj@data)
  obj<-rowSums(obj)
  Genes<-as.character(names(obj))
  Values<-unname(obj)
  obj<-data.frame(Genes=Genes,Values=Values)
} #Takes the row sum for all the genes. Can normalize to number of expressing cells. 

#Set paths and load processed cells 
path<-"/Users/mdolan/Google Drive (mdolan@broadinstitute.org)/FennaMatt_dLGN_scRNAseq/"
setwd(path)
figPath<-"/Users/mdolan/Google Drive (mdolan@broadinstitute.org)/FennaMatt_dLGN_scRNAseq/Figures/"
object<-readRDS("all_Processed.rds")
#astros<-readRDS("astrocytes_Reprocessed.rds")

#Subset and save the microglia for further analysis. 
astros<-SubsetData(object = object, ident.use=c("6", "7"), subset.raw = TRUE, do.clean = TRUE)

#Rerun the analysis pipeline on the mgl data only
#Normalize the data and find the variable genes
astros<-NormalizeData(object = astros, normalization.method = "LogNormalize", scale.factor = 10000)
astros<-FindVariableGenes(object = astros, mean.function = ExpMean, dispersion.function = LogVMR, 
                        do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
astros<-ScaleData(object = astros, genes.use = astros@var.genes, display.progress = TRUE, 
                vars.to.regress = c("percent.mito",'nUMI'), do.par = TRUE, num.cores = 4)
astros<-RunPCA(object = astros, pc.genes = astros@var.genes,pcs.compute = 30 , do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5)
PCElbowPlot(object = astros)
astros<- FindClusters(object = astros, reduction.type = "pca", dims.use = 1:19, 
                    resolution = 0.6, print.output = 0, save.SNN = TRUE)
astros<- RunTSNE(object = astros, dims.use = 1:19, do.fast = TRUE)
TSNEPlot(object = astros, pt.size = 1, do.label = TRUE)
saveRDS(astros, "astrocytes_Reprocessed.rds")

#Run some basic QC on the astros alone 
FeaturePlot(object = astros, features.plot = c("nUMI", "percent.mito", "nGene"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.1)

#Compare different clusters and timepoint 
astros<-SetAllIdent(astros, id = "timep") 
TSNEPlot(object = astros, pt.size = 1, do.label = TRUE)
astros<-SetAllIdent(astros, id = "res.0.6") 
table(astros@meta.data$timep, astros@meta.data$res.0.6)

#Run NMF on the astros. 
astros<-RunNMF(astros,factors.compute = 20,log.norm = T)
pdf(paste0(figPath,"NMF_astros.pdf"))
TSNEPlot(astros, do.label = T, pt.size = 1)
astros<- CurateNMF.seurat(astros, make.tsne = FALSE, feature.plot = TRUE, reduction.embedding = 'tsne',reduction.use = 'nmf',do.reorder = T)
astros<-SetAllIdent(astros, id = "timep") 
TSNEPlot(astros, do.label = T, pt.size = 1)
dev.off()

#Identify marker genes for the different subpopulations
astro.markers<-FindAllMarkers(object = astros, only.pos = T)





