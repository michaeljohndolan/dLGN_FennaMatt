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
remove.zero.genes <- function(exp) {
  all.genes.sums <- rowSums(exp@data) 
  genes.use <- names(all.genes.sums [which(all.genes.sums  >0)])
  exp@raw.data <- exp@raw.data[genes.use, ]
  exp@data <- exp@data[genes.use, ]
  exp
}
selectGenes_median = function (object, alpha.thresh = 0.99, var.thresh = 0.1, combine = "union", 
                               keep.unique = F, capitalize = F, do.plot = T, cex.use = 0.3) 
{
  if (length(var.thresh) == 1) {
    var.thresh <- rep(var.thresh, length(object@raw.data))
  }
  genes.use <- c()
  for (i in 1:length(object@raw.data)) {
    if (capitalize) {
      rownames(object@raw.data[[i]]) <- toupper(rownames(object@raw.data[[i]]))
      rownames(object@norm.data[[i]]) <- toupper(rownames(object@norm.data[[i]]))
    }
    trx_per_cell <- colSums(object@raw.data[[i]])
    gene_expr_mean <- rowMeans(object@norm.data[[i]])
    gene_expr_var <- sparse.var(object@norm.data[[i]])
    nolan_constant <- median((1/trx_per_cell))
    alphathresh.corrected <- alpha.thresh/nrow(object@raw.data[[i]])
    genemeanupper <- gene_expr_mean + qnorm(1 - alphathresh.corrected/2) * 
      sqrt(gene_expr_mean * nolan_constant/ncol(object@raw.data[[i]]))
    genes.new <- names(gene_expr_var)[which(gene_expr_var/nolan_constant > 
                                              genemeanupper & log10(gene_expr_var) > log10(gene_expr_mean) + 
                                              (log10(nolan_constant) + var.thresh[i]))]
    if (do.plot) {
      plot(log10(gene_expr_mean), log10(gene_expr_var), 
           cex = cex.use, xlab = "Gene Expression Mean (log10)", 
           ylab = "Gene Expression Variance (log10)")
      points(log10(gene_expr_mean[genes.new]), log10(gene_expr_var[genes.new]), 
             cex = cex.use, col = "green")
      abline(log10(nolan_constant), 1, col = "purple")
      legend("bottomright", paste0("Selected genes: ", 
                                   length(genes.new)), pch = 20, col = "green")
      title(main = names(object@raw.data)[i])
    }
    if (combine == "union") {
      genes.use <- union(genes.use, genes.new)
    }
    if (combine == "intersection") {
      genes.use <- intersect(genes.use, genes.new)
    }
  }
  if (!keep.unique) {
    for (i in 1:length(object@raw.data)) {
      genes.use <- genes.use[genes.use %in% rownames(object@raw.data[[i]])]
    }
  }
  object@var.genes <- genes.use
  return(object)
}
#Set paths and load processed cells 
path<-"/Users/michaeljohndolan/Google Drive (mdolan@broadinstitute.org)/"
setwd(path)

#Load up the data 
hammond<-readRDS(file = "Hammond2018_microglia_DGE/Round_2_40.filtered.raw.dge.RDS")
#hammond<-CreateSeuratObject(hammond)
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

#Run liger to combine the two datasets. First remove the excess 0 genes in Tim's data. 
table(rowSums(hammond@data)==0)
hammond<-remove.zero.genes(hammond)
table(rowSums(hammond@data)==0)

#Initialize Liger object 
ligerex<-seuratToLiger(list(mgls, hammond), combined.seurat = F, use.tsne = F)

#Normalize, scale (but not center) and find variable genes in the shared dataset. This did nto work, used a workaround
#taking the union of Seurat-determined variable genes 
ligerex = normalize(ligerex)
ligerex = selectGenes(ligerex)
ligerex = scaleNotCenter(ligerex)

#Determine what K and L to use for the optimizeALS function
#suggestK(ligerex) # plot entropy metric to find an elbow that can be used to select the number of factors
#suggestLambda(ligerex, 40) # plot alignment metric to find an elbow that can be used to select the value of lambda

#Perform the factorization
ligerex = optimizeALS(ligerex, k = 40, lambda = 8) 
ligerex = quantileAlignSNF(ligerex) #SNF clustering and quantile alignment

#Visualize the alignment
ligerex = runTSNE(ligerex)
plotByDatasetAndCluster(ligerex, pt.size = 0.2, text.size = 5) #Can also pass in different set of cluster labels to plot

#Examine how marker gene epression looks in the dataset 
plotGene(ligerex, gene = "Spp1")
plotGene(ligerex, gene = "Ccl3")
plotGene(ligerex, gene = "Atf3")
plotGene(ligerex, gene = "Ccl4")

#Examine the different markers and genes highly loading onto each factor with
#word clouds. Factor number 6 has lots of the markers with Ccl4, Ccl3, Atf3 in common
wordclouds = plotWordClouds(ligerex, return.plots = T)
wordclouds[[6]]

#What cluster does factor 6 relate to? Cluster 13
plotClusterFactors(ligerex)
plotClusterProportions(ligerex)

#Need to link the barcodes to cluster 13. Read in different barcodes and link to the liger object 

#NB need to keep these mgl and hammond files separate for the clusters 
mgl.barcode<-select(mgls@meta.data, orig.ident, res.0.6, timep)
mgl.barcode$barcode<-rownames(mgl.barcode)
rownames(mgl.barcode)<-NULL
hammond.barcode<-select(hammond@meta.data, orig.ident, res.0.6)
hammond.barcode$timep<-"P4P5"
hammond.barcode$barcode<-rownames(hammond.barcode)
rownames(hammond.barcode)<-NULL
indiv.barcodes<-rbind(mgl.barcode, hammond.barcode) #Combine into a single barcode df for both experiments 

#Create  a liger barcode df and merge these #FIX THIS SOMETHING WRONG HAPPENED 
ligerex.barcodes<-as.data.frame(ligerex@clusters)
ligerex.barcodes$barcode<-rownames(ligerex.barcodes)
rownames(ligerex.barcodes)<-NULL

total.barcodes<-merge(x = mgl.barcode, y = ligerex.barcodes, by="barcode")
test<-filter(total.barcodes, timep=="P4P5")














#Examine each factor and the gene loading, and print to a PDF
pdf(file = "Hammond_dLGN_align_factors.pdf")
plotFactors(ligerex)
dev.off()


