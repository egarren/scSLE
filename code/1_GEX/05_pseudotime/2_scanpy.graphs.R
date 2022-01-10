rm(list=ls())
library(plyr)
library(Seurat)
library(dplyr)
library(Matrix)
library(knitr)
library(viridis)
library(ggplot2)
library(cowplot)


load("graphed.RData")


#import scanpy diffmap 
scanpy.diffmap <- read.table("scanpy.diffmap.txt", sep = "\t")
scanpy.diffmap <- scanpy.diffmap[,c(2,3)]
colnames(scanpy.diffmap) <- c("scanpy.DC1", "scanpy.DC2")
rownames(scanpy.diffmap) <- SLE.obj.combined@meta.data$cellID
scanpy.diffmap$scanpy.DC1<-scanpy.diffmap$scanpy.DC1*1000
scanpy.diffmap$scanpy.DC2<-scanpy.diffmap$scanpy.DC2*1000
scanpy.diffmap <- as.matrix(scanpy.diffmap)
SLE.obj.combined[["scanpy.diffmap"]] <- CreateDimReducObject(
  embeddings = scanpy.diffmap, key = "scanpyDC_", assay = DefaultAssay(SLE.obj.combined))
scanpy.dpt <- read.table("scanpy.dpt.txt", sep = "\t")
rownames(scanpy.dpt) <- scanpy.dpt[,1]
scanpy.dpt <- scanpy.dpt[,2, drop = FALSE]
colnames(scanpy.dpt) <- "scanpy.pseudotime"
SLE.obj.combined <- AddMetaData(SLE.obj.combined, metadata = scanpy.dpt)
SLE.obj.combined@meta.data$scanpy.pseudo.rank <- rank(SLE.obj.combined@meta.data$scanpy.pseudotime) 
save(SLE.obj.combined,file="scanpy.pseudotime.RData")

#plot
Idents(SLE.obj.combined) <- "my.clusters2"
DimPlot(SLE.obj.combined, reduction= "scanpy.diffmap")+ labs(x="DC_1",y="DC_2")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank()) #
ggsave2("scanpy.diffmap.png",width=6, height=4,device="png")
FeaturePlot(SLE.obj.combined, features= "scanpy.pseudotime", cols= viridis(100, begin = 0), reduction = "scanpy.diffmap")+ 
  labs(x="DC_1",y="DC_2",color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("scanpy.diffmap.pseudo.png",width=5.5, height=4,device="png")
FeaturePlot(SLE.obj.combined, features= "scanpy.pseudotime", cols= viridis(100, begin = 0))+labs(color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("umap.scanpy.pseudo.png",width=5.5, height=4,device="png")
VlnPlot(SLE.obj.combined, features = "scanpyDC_1",pt.size=0)+ NoLegend()+labs(y="DC_1")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.scanpyDC1.png",width=4, height=4,device="png")
VlnPlot(SLE.obj.combined, features = "scanpy.pseudotime",pt.size=0)+ NoLegend()+labs(y="Pseudotime")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.scanpypseudo.png",width=4, height=4,device="png")
VlnPlot(SLE.obj.combined, features = "scanpy.pseudo.rank",pt.size=0)+ NoLegend()+labs(y="Pseudotime Rank")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.scanpypseudo.rank.png",width=4, height=4,device="png")
FeatureScatter(SLE.obj.combined, feature1 = "scanpyDC_1", feature2 = "PC_1")+labs(color="")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("scanpy.DC.PC.png",width=6, height=4,device="png")


