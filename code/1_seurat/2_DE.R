rm(list=ls())
library(Seurat)
library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(biomaRt)
library(plotly)
library(scales)
library(EnhancedVolcano)
library(data.table)
library(ggpubr)
library(limma)
library(VennDiagram)
library(viridis)
library(pheatmap)
library(phylotools)
library(ggforce)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}


load("sce.RData")

## Cluster Analysis
rm(list=setdiff(ls(), c("SLE.obj.combined","metadata","genes.meta")))
Idents(SLE.obj.combined)<-"seurat_clusters"
DefaultAssay(SLE.obj.combined) <- "RNA"
SLE.obj.combined<- NormalizeData(SLE.obj.combined) #if SCTransform data
# Cluster markers
SLE.obj.markers <- FindAllMarkers(object = SLE.obj.combined, test.use = "MAST")
write.csv(SLE.obj.markers,"cluster.markers.csv")
#comparing cluster frequencies
SLE.obj.combined@meta.data$my.clusters <- Idents(SLE.obj.combined)  # Store cluster identities in object@meta.data$my.clusters
cluster.counts<-table(SLE.obj.combined@meta.data$my.clusters, SLE.obj.combined@meta.data$condition)
# Finding conserved markers of clusters
clusters<-names(table(SLE.obj.combined$my.clusters))

#set metadata and subsets
mice<-names(table(SLE.obj.combined$mouse_ID))
Idents(SLE.obj.combined) <- "condition"
SLE.combined<-subset(SLE.obj.combined, idents='SLE.yaa')
WT.combined<-subset(SLE.obj.combined, idents='B6')
WT.mouseID<-metadata[metadata$condition=="B6",]$mouse_ID
SLE.mouseID<-metadata[metadata$condition=="SLE.yaa",]$mouse_ID

##Condition DE
DefaultAssay(SLE.obj.combined) <- "RNA"
Idents(SLE.obj.combined) <- "condition" #setting idents to condition metadata
SLE.obj.combined.autoimmune.response <- FindMarkers(SLE.obj.combined, ident.1 = "SLE.yaa", ident.2 = "B6", 
                                                  min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
avg.SLE.obj.combined <- as.data.frame(log1p(AverageExpression(SLE.obj.combined)$RNA))
avg.SLE.obj.combined$gene <- rownames(avg.SLE.obj.combined)
top30 <- head(SLE.obj.combined.autoimmune.response, n = 30)
SLE.obj.combined.autoimmune.response$log2FC<-log2(exp(SLE.obj.combined.autoimmune.response$avg_log2FC))
write.csv(SLE.obj.combined.autoimmune.response ,"SLEvsWT.csv")
SLE.obj.combined$celltype.condition <- paste(SLE.obj.combined$my.clusters, SLE.obj.combined$condition, sep = "_") #adding metadata identifier


###Within Cluster DE
for (i in clusters){
  Idents(SLE.obj.combined) <- "celltype.condition" #setting idents to new metadata column
  df <- FindMarkers(SLE.obj.combined, ident.1 = paste0(i,"_SLE.yaa"), ident.2 = paste0(i,"_B6"),
                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
  df2 <- head(df, n = 30)
  assign(paste0("cluster",i,"top30"),df2)
  df$log2FC<-log2(exp(df$avg_log2FC))
  write.csv(df,paste0("cluster",i,".autoimmune.markers.csv"))
  assign(paste0("cluster",i,".autoimmune.response"),df)
  Idents(SLE.obj.combined) <- "my.clusters"
  temp<- subset(SLE.obj.combined, idents = i)
  Idents(temp) <- "condition"
  df <- as.data.frame(log1p(AverageExpression(temp)$RNA))
  df$gene <- rownames(df)
  assign(paste0("avg.cluster",i),df)
  assign(paste0("cluster",i),temp)
}
save.image("analyzed.RData")
