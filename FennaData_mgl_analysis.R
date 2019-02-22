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
gene.Sum<-function(obj, norm.to.cell=F) {
  obj<-as.matrix(obj@data)
  obj<-rowSums(obj)
  Genes<-as.character(names(obj))
  Values<-unname(obj)
  obj<-data.frame(Genes=Genes,Values=Values)
} #Takes the row sum for all the genes. Can normalize to number of expressing cells. 
remove.zero.genes <- function(exp) {
  all.genes.sums <- rowSums(exp@data) 
  genes.use <- names(all.genes.sums [which(all.genes.sums  >0)])
  exp@raw.data <- exp@raw.data[genes.use, ]
  exp@data <- exp@data[genes.use, ]
  exp
}

#Set paths and load processed cells 
path<-"/Users/mdolan/Google Drive (mdolan@broadinstitute.org)/FennaMatt_dLGN_scRNAseq/"
setwd(path)
figPath<-"/Users/mdolan/Google Drive (mdolan@broadinstitute.org)/FennaMatt_dLGN_scRNAseq/Figures/"
object<-readRDS("all_Processed.rds")
#mgls<-readRDS("microglia_Reprocessed.rds")

#Subset and save the microglia for further analysis. See mgl script. 
mgls<-SubsetData(object = object, ident.use="9", subset.raw = TRUE, do.clean = TRUE)

#How many genes are 0 accross the whole population of mgls cells
table(rowSums(mgls@data)==0) #6k! 
mgls<-remove.zero.genes(mgls)
table(rowSums(mgls@data)==0)  #0 

#Rerun the analysis pipeline on the mgl data only
#Normalize the data and find the variable genes
mgls<-NormalizeData(object = mgls, normalization.method = "LogNormalize", scale.factor = 10000)
mgls<-FindVariableGenes(object = mgls, mean.function = ExpMean, dispersion.function = LogVMR, 
                        do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
mgls<-ScaleData(object = mgls, genes.use = mgls@var.genes, display.progress = TRUE, 
                  vars.to.regress = c("percent.mito",'nUMI'), do.par = TRUE, num.cores = 4)
mgls<-RunPCA(object = mgls, pc.genes = mgls@var.genes,pcs.compute = 30 , do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
PCElbowPlot(object = mgls)
mgls<- FindClusters(object = mgls, reduction.type = "pca", dims.use = 1:19, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
mgls<- RunTSNE(object = mgls, dims.use = 1:19, do.fast = TRUE)
TSNEPlot(object = mgls, pt.size = 1, do.label = TRUE)
saveRDS(mgls, "microglia_Reprocessed.rds")

#Run some basic QC on the mgls alone 
FeaturePlot(object = mgls, features.plot = c("nUMI", "percent.mito", "nGene"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.1)

#Examine the genes highly loading onto the PCs.
#Cluster 2 is proliferative. Cluster 3 is a weird bunch of likely doublets/misclustered/ATM-engulfing cells. 
#Cluster 1 and 0 seem clearly microglial
all.marker.list<-FindAllMarkers(mgls)
Cluster1.marker.list<-filter(all.marker.list, cluster=="1")
View(Cluster1.marker.list)

#Run NMF on the dataset to try find alternative groupings. Cluster 1 has expression of several gnes that may match 
#an aging cluster in Tim Hammond's dataset. Plot the tSNE by clusters created by NMF. 
mgls<-RunNMF(mgls,factors.compute = 20,log.norm = T)
pdf(paste0(figPath,"NMF_microglia.pdf"))
TSNEPlot(mgls, do.label = T, pt.size = 1.5)
mgls<- CurateNMF.seurat(mgls, make.tsne = FALSE, feature.plot = TRUE, reduction.embedding = 'tsne',reduction.use = 'nmf',do.reorder = T)
TSNEPlot(mgls, do.label = T, pt.size = 1.5)
dev.off()

#Plot marker genes for each different major cluster (0,1,2). Also plot as feature or dotplots
VlnPlot(mgls, features.plot = c("Ccl4", "Ccl3", "Nfkbia", "Atf3", "Plin2", "Sdc4", "Mt1"), point.size.use = 0.5)
VlnPlot(mgls, features.plot = c("Top2a", "Mki67", "Hmmr", "Cenpf", "Prc1", "Cenpe", "Smc4"), point.size.use = 0.5)
VlnPlot(mgls, features.plot = c("Siglech", "Csf1r", "Maf"), point.size.use = 0.5)

#Look at colocalization of different genes for cluster number 1
FeaturePlot(object = mgls, features.plot = c("Ccl3", "Spp1")
            , cols.use = c("grey","red", "blue", "green"), overlay = TRUE, no.legend = FALSE)

#Examine the correlations with Ccl3
corr.gene<-"Ccl3"
matrix<-mgls@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod[corr.gene,])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations<-na.omit(correlations)
View(correlations)

#How do the number of cells per cluster compare accross different timepoints and plot these. 
proportions<-as.data.frame(table(mgls@meta.data[,c(5,7)]))
proportions<-dcast(data = proportions, res.0.6~timep)
proportions<-melt(proportions)
names(proportions)<-c("Cluster", "Timepoint", "NumberCells")
proportions<-arrange(proportions, desc(Timepoint))

g<-ggplot(proportions, aes(x =  reorder(Timepoint), y = NumberCells, fill = Cluster))
g<-g+geom_bar(stat = "identity", color="Black")
g<-g+xlab("Timepoint")+ylab("Number of Cells")
g

percent.proportions<-mutate(proportions, TotalCells=P10+P5, P5.percent=(P5/TotalCells)*100
                    ,P10.percent=(P10/TotalCells)*100 )

#Find some differences that predict P5 or P10 mgls 
mgls<-SetAllIdent(mgls, id = "timep") 
mgls.ROC<-FindAllMarkers(object = mgls)
mgls<-SetAllIdent(mgls, id = "res.0.6") 

#Generate PDF of different microglial analyses
pdf(paste0(figPath,"analysis_microglia.pdf"))
mgls<-SetAllIdent(mgls, id = "res.0.6") 
TSNEPlot(object = mgls, pt.size = 0.75, do.label = FALSE)
mgls<-SetAllIdent(mgls, id = "timep") 
TSNEPlot(object = mgls, pt.size = 0.65, do.label = FALSE) #P5 v P10
mgls<-SetAllIdent(mgls, id = "res.0.6") 
g #Plot the proportions of different cell-types at different times. 
FeaturePlot(object = mgls, features.plot = c("nUMI", "percent.mito", "nGene"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.3)
VlnPlot(mgls, features.plot = c("Ccl4", "Ccl3", "Nfkbia", "Atf3", "Plin2", "Sdc4", "Mt1"), point.size.use = 0.5)
VlnPlot(mgls, features.plot = c("Top2a", "Mki67", "Hmmr", "Cenpf", "Prc1", "Cenpe", "Smc4"), point.size.use = 0.5)
VlnPlot(mgls, features.plot = c("Siglech", "Csf1r", "Maf", "Tgfbr1", "4632428N05Rik", "Lpar6"), point.size.use = 0.5)
#Look at activation signal across different timepoints, differentially expressed at timepoints suggests it's real signal
VlnPlot(mgls, features.plot = c("Fos", "Egr1", "Ier2", "Ccl3", "Ccl12"), point.size.use = 0.5)
mgls<-SetAllIdent(mgls, id = "timep") 
VlnPlot(mgls, features.plot = c("Fos", "Egr1", "Ier2", "Ccl3", "Ccl12", "Jun"), point.size.use = 0.5)
dev.off()


