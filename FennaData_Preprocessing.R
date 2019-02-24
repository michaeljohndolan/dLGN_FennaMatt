#Code to merge the samples for each timepoint and process, analyze and save all the data. 
#This was run on google cloud 

#Load up required libraries and custom functions
library(Seurat)
library(Matrix)
library(dplyr)
QC.plotter<-function(object) {
  data<-object@meta.data
  median.mito<-median(data$percent.mito); print(paste0("Median mito is: ", median.mito))
  median.nUMI<-median(data$nUMI); print(paste0("Median nUMI is: ", median.nUMI))
  
  g<-ggplot(data, aes(x=nUMI, y=percent.mito))
  g<-g+geom_point(size=0.5, alpha=0.05)
  g<-g+coord_cartesian(expand = F)
  g<-g+geom_hline(yintercept = median.mito, color="red")
  g<-g+geom_vline(xintercept = median.nUMI, color="green")
  g
}
name.extract<-function(name) {
  name<-sapply(strsplit(name, "\\Q.\\E"), "[", 1)
  name
}

#Set count matrix data path (tsv files, not included in repo as very large)
path<-"/Users/mdolan/Google Drive (mdolan@broadinstitute.org)/FennaMatt_dLGN_scRNAseq/"
setwd(path)
samples<-list.files(pattern = "*.txt")

#Loop through each file iteratively and open then merge
#Note the min.cells may change as you filter out subsequently
for(i in 1:length(samples)) {
  if(i==1) {
    #Will initialize a Seurat object for the first file 
    object<-read.table(samples[i],header = TRUE,quote = "",skip = 0,row.names = 1)
    object<-CreateSeuratObject(raw.data = object , min.cells = 3, min.genes = 200, 
                               project =name.extract(samples[i])) 
  }

  if(i>1) {
    temp<-read.table(samples[i],header = TRUE,quote = "",skip = 0,row.names = 1)
    temp<-CreateSeuratObject(raw.data = temp , min.cells = 3, min.genes = 200, 
                             project =name.extract(samples[i]))
    object<-MergeSeurat(object1 = object, object2 = temp, do.normalize = FALSE,add.cell.id1 = i-1 ,add.cell.id2 = i)
    rm(temp)
  } #Will merge subsequent Seurat samples
}

#First part of cell QC, calculate the percent mito genes expressed for each cell and add this to the metadata slot. 
mito.genes <- grep(pattern = "^mt.", x = rownames(x = object@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ])/Matrix::colSums(object@raw.data)
object<-AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")

#Plot cellQC parameters from the metadata slot with violin plot and geneplot. 
VlnPlot(object = object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, point.size.use=0.005)
par(mfrow = c(1, 2))
GenePlot(object = object, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.05)
GenePlot(object = object, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.05)
median(object@meta.data$percent.mito)
median(object@meta.data$nUMI)
median(object@meta.data$nGene)

#Filter the dataset by cells
nrow(object@meta.data)
object<- FilterCells(object = object, subset.names = c("nGene", "nUMI","percent.mito"), 
                         low.thresholds = c(200, 500, -Inf), high.thresholds = c(2500, 15000, 0.1))
nrow(object@meta.data)

#Normalize the data and find the variable genes
object <- NormalizeData(object = object, normalization.method = "LogNormalize", scale.factor = 10000)
object<-FindVariableGenes(object = object, mean.function = ExpMean, dispersion.function = LogVMR, 
                              do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#Scale the data accross the different genes and regress out variations due to nUMI (high expression) and per.cent mito
# (genes correlated with mt expression). Might play around with this b/c presumable microglia during pruning use up 
# a lot of energy. 
object <- ScaleData(object = object, genes.use = object@var.genes, display.progress = TRUE, 
                        vars.to.regress = c("percent.mito",'nUMI'), do.par = TRUE, num.cores = 4)

#Run the PCA. Keep PCs 1-30. 
object<- RunPCA(object = object, pc.genes = object@var.genes,pcs.compute = 30 , do.print = TRUE, pcs.print = 1:5, 
                    genes.print = 5)
PCElbowPlot(object = object)

#Cluster the data with 19 PCs
object <- FindClusters(object = object, reduction.type = "pca", dims.use = 1:19, 
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)

## Run and plot tSNE
object <- RunTSNE(object = object, dims.use = 1:19, do.fast = TRUE)
TSNEPlot(object = object, pt.size = 0.05, do.label = TRUE)

#Annotate the different rows by timepoint and sex 
object@meta.data$timep<-"timep"
object@meta.data$sex<-"sex"
for(i in 1:nrow(object@meta.data)) {
  id<-object@meta.data[i,]$orig.ident
  if(grepl(x = id,pattern = "P5", fixed = TRUE)==TRUE) object@meta.data[i,]$timep<-"P5"
  else object@meta.data[i,]$timep<-"P10"
  if(grepl(x = id,pattern = "f", fixed = TRUE)==TRUE) object@meta.data[i,]$sex<-"female"
  else object@meta.data[i,]$sex<-"male"
}

#Save the final object with tSNE calculated. Move to second script for detailed analysis of individual cell-types. 
saveRDS(object, "all_Processed.rds")
