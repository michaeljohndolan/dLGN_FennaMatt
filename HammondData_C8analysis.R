#Examine cluster 8 from Hammond (2018) accross multiple different timepoints. Will use this data to parse the Liger 
#alignment of Fenna and Hammond datasets. 
#Note not sure if animal column is correct as there are cells in these factors accross multiple timepoints

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)

#Initial df processing 
Tim.Mgl.assign<-readRDS("Google Drive (mdolan@broadinstitute.org)/Hammond2018_microglia_DGE/Round_2_40.cluster.assign.RDS")
Tim.Mgl.assign<-as.data.frame(Tim.Mgl.assign)
Tim.Mgl.assign$Barcodes<-row.names(Tim.Mgl.assign)
row.names(Tim.Mgl.assign) <- NULL

#Extract the timepoint, sex and genotype information and create 2nd df consisting only of cluster8
Tim.Mgl.assign$Timepoint<-sapply(strsplit(Tim.Mgl.assign$Barcodes, "_"), "[", 1)
Tim.Mgl.assign$animal<-sapply(strsplit(Tim.Mgl.assign$Barcodes, "_"), "[", 3)
Tim.Mgl.assign<-select(Tim.Mgl.assign, -Barcodes)
Tim.Mgl.assign<-filter(Tim.Mgl.assign, animal!="C1")
Tim.Mgl.assign<-filter(Tim.Mgl.assign, animal!="LPC")
Tim.Mgl.assign<-filter(Tim.Mgl.assign, animal!="SALINE")

#Create total number of mgls per timepoint (aggregated)
total.mgls<-table(Tim.Mgl.assign$Timepoint) #Get the total number of mgls
total.mgls<-as.data.frame(total.mgls)
names(total.mgls)<-c("Timepoint", 'total.Freq')

#Create total number of C8 mgls per timepoint (aggregated)
c8<-filter(Tim.Mgl.assign, Tim.Mgl.assign=="13")
c8<-select(.data = c8, Tim.Mgl.assign, Timepoint)
c8$Tim.Mgl.assign<-as.character(c8$Tim.Mgl.assign)
c8<-table(c8)
c8<-as.data.frame(c8)

#Calculate the percent of total and merge P4/P5 as P5 has too few cells
Percent.c8<-merge(x = total.mgls, y=c8, by="Timepoint")
Percent.c8<-mutate(Percent.c8, percent=(Freq/total.Freq)*100)
Percent.c8$Timepoint<-as.character(Percent.c8$Timepoint)
Percent.c8<-rbind(Percent.c8, list("P4P5", 14992+2336, 13, 424+107,0))
Percent.c8<-mutate(Percent.c8, percent=(Freq/total.Freq)*100)
Percent.c8<-filter(Percent.c8, Timepoint!="P4")
Percent.c8<-filter(Percent.c8, Timepoint!="P5")

#Reorder the factors
Percent.c8$Timepoint <- factor(Percent.c8$Timepoint, levels = c("E14", "P4P5", "P30", "P100", "Old"))

#Plot this data on simple chart 
g<-ggplot(data=Percent.c8, aes(x=Timepoint, y=percent))
g<-g+geom_bar(stat="identity")
g<-g+ylab("Percentage Cluster 8")
g

#ONCE I WORK OUT THE ANIMAL CODES/IDs
#First calculate the total number of microglia for each animal 
table(Tim.Mgl.assign$Tim.Mgl.assign,Tim.Mgl.assign$animal)


