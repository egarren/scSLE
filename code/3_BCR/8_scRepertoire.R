rm(list=ls())
# Load required packages
library(Seurat)
library(scRepertoire)
library(cowplot)


#load cdr3 data
metadata<-read.csv("SLE.metadata.csv", header=T) #define path
metadata[] <- lapply(metadata, as.character)
mice<-metadata$mouse_ID
contig_list<-list()
for (i in mice){
  df<-read.csv(file=paste0("./cellranger/",i,"/filtered_contig_annotations.csv"),header=T,stringsAsFactors = F)
  contig_list[[i]]<-df
}
head(contig_list[[4]])
combined<-combineBCR(contig_list,samples=names(contig_list),ID=metadata$condition)
# example <- addVariable(combined, name = "genotype",  variables = metadata$condition)

#visualize
quantContig(combined, cloneCall="gene+nt", scale = T,group="ID")#exportTable = T)
ggsave2("unique.png",width=5, height=4,device="png")
abundanceContig(combined, cloneCall = "gene", scale = F,group="ID")
ggsave2("abundance.png",width=5, height=4,device="png")
lengthContig(combined, cloneCall="aa", chains = "combined",group="ID") 
ggsave2("aa.length.png",width=8, height=4,device="png")
lengthContig(combined, cloneCall="nt", chains = "single",group="ID") 
ggsave2("nt.length.png",width=8, height=4,device="png")
clonalHomeostasis(combined, cloneCall = "gene")
ggsave2("gene.homeostasis.png",width=5, height=4,device="png")
clonalHomeostasis(combined, cloneCall = "aa")
ggsave2("aa.homeostasis.png",width=5, height=4,device="png")
clonalProportion(combined, cloneCall = "gene") 
ggsave2("gene.proportion.png",width=5, height=4,device="png")
clonalProportion(combined, cloneCall = "nt") 
ggsave2("nt.proportion.png",width=5, height=4,device="png")
clonalOverlap(combined, cloneCall = "gene+nt", method = "morisita")
ggsave2("overlap.png",width=5, height=4,device="png")
png("size.distribution.png")
clonesizeDistribution(combined, cloneCall = "gene+nt", method="ward.D2")
dev.off()


#diversity
clonalDiversity(combined, cloneCall = "gene", group = "samples",n.boots = 100)
ggsave2("diversity.sample.png",width=8, height=4,device="png")
clonalDiversity(combined, cloneCall = "gene", group = "ID")
ggsave2("diversity.condition.png",width=8, height=4,device="png")

#Seurat
load("graphed.RData")

seurat<-SLE.obj.combined
#rename barcodes
combined2<-list()
for(i in 1:length(combined)){
  df<-combined[[i]]
  df$barcode<-gsub("_B6|_SLE.yaa","",df$barcode)
  combined2[[i]]<-df
}
#combine w/Seurat
seurat <- combineExpression(combined2, seurat, 
                            cloneCall="gene", groupBy = "none", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=2, Medium=3, Large=4, Hyperexpanded=5))



colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(seurat, group.by = "condition") + NoLegend() +
  scale_color_manual(values=colorblind_vector(2))
DimPlot(seurat, group.by = "cloneType",split.by="condition") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")
ggsave2("umap.clonesize.png",width=15, height=6,device="png")
clonalOverlay(seurat, reduction = "umap", freq.cutpoint = 5, bins = 10, facet = "condition") + guides(color = FALSE)
ggsave2("umap.overlay.png",width=12, height=6,device="png")
head(sort(table(seurat@meta.data$CTaa), decreasing=T))
seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CMRYGNYWYFDVW_CLQHGESPFTF", "NA_CQQYNSYPLTF"))
DimPlot(seurat, group.by = "highlight")
ggsave2("umap.selectclones.png",width=7, height=6,device="png")
occupiedscRepertoire(seurat, x.axis = "cluster")
ggsave2("repertoire.cluster.png",width=7, height=6,device="png")
head(sort(table(seurat@meta.data$CTgene), decreasing=T))
alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("condition", "my.clusters2", "cloneType"), 
                   color = "IGHV1-26.IGHJ2..IGHM_NA") + 
  scale_fill_manual(values = c("grey", colorblind_vector(1)))
ggsave2("flow.clone.png",width=7, height=6,device="png")
png("flow.cluster.png",width=7,height=6,units="in",res=300)
alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("condition", "my.clusters2", "cloneType"), 
                   color = "my.clusters2") 
dev.off()
StartracDiversity(seurat, type = "cloneType", sample = "mouse_ID", by = "overall")
ggsave2("diversity.png",width=7, height=6,device="png")

#cluster analysis
combined2 <- expression2List(seurat, group = "my.clusters2")
length(combined2) #now listed by cluster
clonalDiversity(combined2, cloneCall = "nt")
ggsave2("diversity.cluster.png",width=7, height=6,device="png")
clonalHomeostasis(combined2, cloneCall = "nt")
ggsave2("homeosasis.cluster.png",width=7, height=6,device="png")
clonalProportion(combined2, cloneCall = "nt")
ggsave2("proportion.cluster.png",width=7, height=6,device="png")
clonalOverlap(combined2, cloneCall="aa", method="morisita")
ggsave2("overlap.cluster.png",width=6, height=4.5,device="png")
png("size.distribution.cluster.png",width=7,height=6,units="in",res=300)
clonesizeDistribution(combined2, cloneCall = "aa", method="ward.D2")
dev.off()





