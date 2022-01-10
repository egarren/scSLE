require(EWCE)
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
library(Seurat)
library(plyr)
library(readr)

load("graphed.RData")

alldata<-SLE.obj.combined 
alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000) #log normalize
matrix <- as.data.frame(alldata@assays$RNA@data) #extract counts from Seurat
matrix2 <- tibble::rownames_to_column(matrix, "Gene Symbol") #genes to column (so saved in TSV)
write_tsv(matrix2,'expression.tsv') #save tsv

compass.meta<- cbind(data.frame(cell_id=rownames(alldata@meta.data),
                          cell_type=alldata@meta.data$my.clusters2),
                     alldata@meta.data[,c("mouse_ID","condition","gender")])

write.csv(compass.meta,file="compass.meta.csv")
