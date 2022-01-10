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


load("analyzed.RData")
dir.create("./Bcells")
setwd("./Bcells")

#remove non-SC clusters (change based on which subset)
Idents(SLE.obj.combined) <- "my.clusters"
# cluster.filter.keep<-c("1","3")
# cluster.filter.keep<-c("2","5","7","8")
# cluster.filter.keep<-c("0","5","7","9")
SLE.B.obj<-subset(SLE.obj.combined, idents=cluster.filter.keep)
SLE.B.obj$SCT<-NULL
SLE.B.obj$integrated<-NULL
SLE.B.obj$umap<-NULL
SLE.B.obj$tsne<-NULL

##Integrating by condition
#split by condition 
SLE.B.obj.list <- SplitObject(SLE.B.obj, split.by = "condition")
# SCTransform
for (i in 1:length(SLE.B.obj.list)) {
  SLE.B.obj.list[[i]] <- SCTransform(SLE.B.obj.list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt","HSP.score1",
                                                                          "percent.Rpl","percent.Rps")) #", "CC.Difference",
}
#integrate and anchor
SLE.B.obj.features <- SelectIntegrationFeatures(object.list = SLE.B.obj.list, nfeatures = 3000)
SLE.B.obj.list <- PrepSCTIntegration(object.list = SLE.B.obj.list, anchor.features = SLE.B.obj.features)
SLE.B.obj.anchors <- FindIntegrationAnchors(object.list = SLE.B.obj.list, normalization.method = "SCT",anchor.features = SLE.B.obj.features)
SLE.B.obj.combined <- IntegrateData(anchorset = SLE.B.obj.anchors, normalization.method = "SCT")
SLE.B.obj.combined <- RunPCA(SLE.B.obj.combined)
rm(list=setdiff(ls(), c("SLE.B.obj.combined", "metadata","genes.meta")))

ElbowPlot(object = SLE.B.obj.combined)
ggsave2("elbow.png",device="png")
SLE.B.obj.combined <- RunUMAP(SLE.B.obj.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
SLE.B.obj.combined <- FindNeighbors(SLE.B.obj.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
SLE.B.obj.combined <- FindClusters(SLE.B.obj.combined, resolution = 0.13) #adjust resolution (bigger=more clusters), initially used 0.2
SLE.B.obj.combined <- RunTSNE(object = SLE.B.obj.combined, dims.use = 1:25, do.fast = TRUE) #change number of PCs to use, change perplexity (https://distill.pub/2016/misread-tsne/)
DimPlot(SLE.B.obj.combined, reduction = "umap",label=T)

#examine clusters
Idents(SLE.B.obj.combined)<-"seurat_clusters"
DefaultAssay(SLE.B.obj.combined) <- "RNA"
SLE.B.obj.combined<- NormalizeData(SLE.B.obj.combined) #if SCTransform data
# Cluster markers
SLE.B.obj.markers <- FindAllMarkers(object = SLE.B.obj.combined, test.use = "MAST")
write.csv(SLE.B.obj.markers,"cluster.markers.csv")
top4 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
top1 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
DefaultAssay(SLE.B.obj.combined) <- "integrated"
DoHeatmap(object = SLE.B.obj.combined, features = top4$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=4, height=6,device="png")
DefaultAssay(SLE.B.obj.combined) <- "RNA"
SLE.B.obj.combined<- NormalizeData(SLE.B.obj.combined) #if SCTransform data

# Cluster markers
SLE.B.obj.markers <- FindAllMarkers(object = SLE.B.obj.combined, test.use = "MAST")
write.csv(SLE.B.obj.markers,"cluster.markers.csv")
top4 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
top1 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
DefaultAssay(SLE.B.obj.combined) <- "integrated"
DoHeatmap(object = SLE.B.obj.combined, features = top4$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=4, height=6,device="png")
DefaultAssay(SLE.B.obj.combined) <- "RNA"


for (i in c("condition","Phase","mouse_ID","gender","batch","seurat_clusters")){
  Idents(SLE.B.obj.combined)<-i
  DimPlot(SLE.B.obj.combined, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(SLE.B.obj.combined, features= i,split.by = "condition",pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap.png"),width=10, height=5,device="png")
}
Idents(SLE.B.obj.combined)<-"seurat_clusters"
save.image("merged.RData")

clusters<-as.data.frame(table(SLE.B.obj.combined@meta.data$seurat_clusters))
letters<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
for(i in 1:nrow(clusters)){
  SLE.B.obj.combined@meta.data$cluster_letter[SLE.B.obj.combined@meta.data$seurat_clusters%in%clusters$Var1[i]]<-letters[i]
}
SLE.B.obj.combined.sce <- as.SingleCellExperiment(SLE.B.obj.combined)
save(SLE.B.obj.combined.sce,file="sce.clusters.RData")

## Cluster Analysis
rm(list=setdiff(ls(), c("SLE.B.obj.combined","metadata","genes.meta")))
Idents(SLE.B.obj.combined)<-"seurat_clusters"
DefaultAssay(SLE.B.obj.combined) <- "RNA"
# Cluster markers
SLE.B.obj.markers <- FindAllMarkers(object = SLE.B.obj.combined, test.use = "MAST")
write.csv(SLE.B.obj.markers,"cluster.markers.csv")
top4 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
top1 <- SLE.B.obj.markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
DefaultAssay(SLE.B.obj.combined) <- "integrated"
DoHeatmap(object = SLE.B.obj.combined, features = top4$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=4, height=6,device="png")
DefaultAssay(SLE.B.obj.combined) <- "RNA"
#comparing cluster frequencies
SLE.B.obj.combined@meta.data$my.clusters <- Idents(SLE.B.obj.combined)  # Store cluster identities in object@meta.data$my.clusters
cluster.counts<-table(SLE.B.obj.combined@meta.data$my.clusters, SLE.B.obj.combined@meta.data$condition)
# Finding conserved markers of clusters
clusters<-names(table(SLE.B.obj.combined$my.clusters))

#set metadata and subsets
mice<-names(table(SLE.B.obj.combined$mouse_ID))
Idents(SLE.B.obj.combined) <- "condition"
SLE.combined<-subset(SLE.B.obj.combined, idents='SLE.yaa')
WT.combined<-subset(SLE.B.obj.combined, idents='B6')
WT.mouseID<-metadata[metadata$condition=="B6",]$mouse_ID
SLE.mouseID<-metadata[metadata$condition=="SLE.yaa",]$mouse_ID

##Condition DE
DefaultAssay(SLE.B.obj.combined) <- "RNA"
Idents(SLE.B.obj.combined) <- "condition" #setting idents to condition metadata
SLE.B.obj.combined.autoimmune.response <- FindMarkers(SLE.B.obj.combined, ident.1 = "SLE.yaa", ident.2 = "B6", 
                                                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
avg.SLE.B.obj.combined <- as.data.frame(log1p(AverageExpression(SLE.B.obj.combined)$RNA))
avg.SLE.B.obj.combined$gene <- rownames(avg.SLE.B.obj.combined)
top30 <- head(SLE.B.obj.combined.autoimmune.response, n = 30)
SLE.B.obj.combined.autoimmune.response$log2FC<-log2(exp(SLE.B.obj.combined.autoimmune.response$avg_log2FC))
write.csv(SLE.B.obj.combined.autoimmune.response ,"SLEvsWT.csv")
SLE.B.obj.combined$celltype.condition <- paste(SLE.B.obj.combined$my.clusters, SLE.B.obj.combined$condition, sep = "_") #adding metadata identifier


###Within Cluster DE
for (i in clusters){
  Idents(SLE.B.obj.combined) <- "celltype.condition" #setting idents to new metadata column
  df <- FindMarkers(SLE.B.obj.combined, ident.1 = paste0(i,"_SLE.yaa"), ident.2 = paste0(i,"_B6"), 
                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
  if (!is(df, "try-error")){
    df2 <- head(df, n = 30)
    assign(paste0("cluster",i,"top30"),df2)
    df$log2FC<-log2(exp(df$avg_log2FC))
    write.csv(df,paste0("cluster",i,".autoimmune.markers.csv"))
    assign(paste0("cluster",i,".autoimmune.response"),df)
    Idents(SLE.B.obj.combined) <- "my.clusters"
    temp<- subset(SLE.B.obj.combined, idents = i)
    Idents(temp) <- "condition"
    df <- log1p(AverageExpression(temp)$RNA)
    df$gene <- rownames(df)
    assign(paste0("avg.cluster",i),df)
    assign(paste0("cluster",i),temp)
  } 
}
save.image("analyzed.RData")