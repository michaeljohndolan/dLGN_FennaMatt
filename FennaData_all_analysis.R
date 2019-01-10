#Code to examine the clusters identified in Fenna's dataset 
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
gene.Sum<-function(obj, norm.to.cell=F) {
  obj<-as.matrix(obj@data)
  obj<-rowSums(obj)
  Genes<-as.character(names(obj))
  Values<-unname(obj)
  obj<-data.frame(Genes=Genes,Values=Values)
} #Takes the row sum for all the genes. Can normalize to number of expressing cells. 
cal.std.residuals<-function(df) {
  rownames(df)<-df$Genes
  lm<-lm(data =df[,2:3])
  tidy(lm)
  augment(lm)
} #Uses functions from broom to calculate a linear model 
plot.diff<-function(aug, outlier=outlier) {
  Value1<-names(aug)[2]
  Value2<-names(aug)[3]
  names(aug)[2:3]<-c("val1", "val2")
  g<-ggplot(data =aug, aes(x=val1, y=val2))
  g<-g+geom_point(alpha=0.3, size=0.5)
  g<-g+geom_smooth(method = "lm")
  g<-g+geom_text(size=2
                 ,aes(label=ifelse(.std.resid>outlier |.std.resid<(outlier*-1) ,as.character(.rownames),''))
                 ,hjust=1, vjust=0)
  g<-g+labs(x=Value1, y=Value2 ,title="") #Titles and labels 
  return(g)
} #Plot the comparision datasets with a regression line and label the outliers 
calculate.diff.expression1<-function(obj, outlier=7) {
  #Calcuate the rowsums and normalized rowsums for each gene in a new object. Use normalized data but not scaled and centered
  P5m<-SubsetByPredicate(object = obj, vars.use = c("timep", "sex"),
                         predicate = "timep=='P5' & sex=='male'")
  P5f<-SubsetByPredicate(object = obj, vars.use = c("timep", "sex"),
                         predicate = "timep=='P5' & sex=='female'")
  P10m<-SubsetByPredicate(object = obj, vars.use = c("timep", "sex"),
                          predicate = "timep=='P10' & sex=='male'")
  P10f<-SubsetByPredicate(object = obj, vars.use = c("timep", "sex"),
                          predicate = "timep=='P10' & sex=='female'")
  P5<-SubsetByPredicate(object = obj, vars.use = c("timep"),
                        predicate = "timep=='P5'")
  P10<-SubsetByPredicate(object = obj, vars.use = c("timep"),
                         predicate = "timep=='P10'")
  P5m.genesum<-gene.Sum(P5m)
  P5f.genesum<-gene.Sum(P5f)
  P10m.genesum<-gene.Sum(P10m)
  P10f.genesum<-gene.Sum(P10f)
  P5.genesum<-gene.Sum(P5)
  P10.genesum<-gene.Sum(P10)
  rm(P5m);rm(P5f);rm(P10m);rm(P10f);rm(P5);rm(P10) #Do this when running on laptop for space
  
  #Create comparision datasets 
  P5_P10<-merge(x = P5.genesum, y = P10.genesum, by = "Genes", no.dups = TRUE
                , suffixes=c(".P5", ".P10"))
  P5m_P10m<-merge(x = P5m.genesum, y = P10m.genesum, by = "Genes", no.dups = TRUE
                  , suffixes=c(".P5m", ".P10m"))
  P5f_P10f<-merge(x = P5f.genesum, y = P10f.genesum, by = "Genes", no.dups = TRUE
                  , suffixes=c(".P5f", ".P10f"))
  P5m_P5f<-merge(x = P5m.genesum, y = P5f.genesum, by = "Genes", no.dups = TRUE
                 , suffixes=c(".P5m", ".P5f"))
  P10m_P10f<-merge(x = P10m.genesum, y = P10f.genesum, by = "Genes", no.dups = TRUE
                   , suffixes=c(".P10m", ".P10f"))
  
  #Calculate the std residuals based on a linear model. 
  P5_P10.aug<-cal.std.residuals(P5_P10)
  P5m_P10m.aug<-cal.std.residuals(P5m_P10m)
  P5f_P10f.aug<-cal.std.residuals(P5f_P10f)
  P5m_P5f.aug<-cal.std.residuals(P5m_P5f)
  P10m_P10f.aug<-cal.std.residuals(P10m_P10f)
  
  return(plot.diff(P5_P10.aug, outlier=outlier))
  #plot.diff(P5m_P10m.aug)
  #plot.diff(P5f_P10f.aug)
  #plot.diff(P5m_P5f.aug)
}

path<-"/Users/mdolan/Google Drive (mdolan@broadinstitute.org)/FennaMatt_dLGN_scRNAseq/"
figPath<-"/Users/mdolan/Google Drive (mdolan@broadinstitute.org)/FennaMatt_dLGN_scRNAseq/Figures/"
object<-readRDS("all_Processed.rds")

#Open up a PDF file to print figures 
pdf(paste0(figPath,"analysis_allCells.pdf"))
TSNEPlot(object = object, pt.size = 0.1, do.label = TRUE)

#Quality control on the overall dataset 
FeaturePlot(object = object, features.plot = c("nUMI", "nGene", "percent.mito"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.2)

#Using some well-described markers, identify different cell-types in this dataset. Microglia are clearly 9. Macrophages 13
FeaturePlot(object = object, features.plot = c("C1qa", "Fcrls", "Mrc1", "Tmem119"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.2)
#Clusters 6-7 appear to be astrocytes.
FeaturePlot(object = object, features.plot = c("Gfap"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.5)
#Clusters 1,5,4,0 and 3 are neuronal
FeaturePlot(object = object, features.plot = c("Snap25"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.5)
#Cluster 11 and 8 are oligodendrocyte while 2 and 10 are polydendrocytes 
FeaturePlot(object = object, features.plot = c("Plp1", "Pdgfra"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.2)

#Examine the distributions of cells at different timepoints and the distributions of cell-types 
object<-SetAllIdent(object = object, id = "timep") 
#object<-SetAllIdent(object = object, id = "res.0.6")
TSNEPlot(object = object, pt.size = 0.1, do.label = FALSE)

#Close printing to PDF 
dev.off()

#What are good markers for macrophages versus microglia 
macro_v_mgls<-FindMarkers(object, ident.1 = 9, ident.2 = 13)
View(macro_v_mgls)

#Run NMF on the dataset to identify clear factors. Need to run this on the cloud
object<-RunNMF(object,factors.compute = 18,log.norm = T)
pdf(paste0(figPath,"NMF_allCells.pdf"))
object<- CurateNMF.seurat(object, make.tsne = F, feature.plot = T, reduction.embedding = 'tsne',reduction.use = 'nmf',do.reorder = T)
TSNEPlot(object, do.label = T, pt.size = 0.5)
dev.off()

##Diff expression between timepoints
#Calculate a differential expression PDF with all the details for each cell-type (run Mgl code first)
mgls<-readRDS("microglia_Reprocessed.rds")
pdf(paste0(figPath,"diff_analysis.pdf"))
calculate.diff.expression1(obj = object, outlier = 7)
calculate.diff.expression1(obj = mgls, outlier = 7)
dev.off()
