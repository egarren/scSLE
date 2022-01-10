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
library(data.table)
library(ggalluvial)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#load data (from https://www.immport.org/shared/study/SDY997)
raw.count<-read.table("SDY997_EXP15176_celseq_matrix_ru10_molecules.tsv",sep="\t",header=T)
rownames(raw.count)<-raw.count$gene
raw.count$gene<-NULL
raw.count[is.na(raw.count)] <- 0
raw.meta<-read.table("SDY997_EXP15176_celseq_meta.tsv",sep="\t",header=T)
rownames(raw.meta)<-raw.meta$cell_name
raw.meta$cell_name<-NULL
obj<-CreateSeuratObject(counts = raw.count,meta.data=raw.meta,min.cells = 3, min.features  = 200,)

##QC and Filter
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
Idents(obj)<-"disease"
VlnPlot(object = obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0, ncol = 3)
ggsave2("vln.QC.png",device="png")
hist(obj@meta.data$nFeature_RNA,breaks=200)
hist(obj@meta.data$nCount_RNA,breaks=200,xlim=c(0,5000))
hist(obj@meta.data$percent.mt,breaks=200)
obj
table(obj@meta.data$disease)
obj<- subset(x = obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & 
                 percent.mt >  -Inf & percent.mt < 30 & type=="Leukocyte") 
obj
table(obj@meta.data$disease)

##Integrating by condition
obj@meta.data$condition<-obj@meta.data$disease
obj.list <- SplitObject(obj, split.by = "condition")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- SCTransform(obj.list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt")) #", "CC.Difference",
}
obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",anchor.features = obj.features)
obj.combined <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")
obj.combined <- RunPCA(obj.combined)

##Clustering
ElbowPlot(object = obj.combined)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
obj.combined <- FindClusters(obj.combined, resolution = 0.15) #adjust resolution (bigger=more clusters), initially used 0.2
obj.combined <- RunTSNE(object = obj.combined, dims.use = 1:25, do.fast = TRUE) #change number of PCs to use, change perplexity (https://distill.pub/2016/misread-tsne/)
DimPlot(obj.combined, reduction = "umap",split.by="condition")
ggsave2("umap.condition.png",width=9, height=5,device="png")
for (i in c("condition","sample","seurat_clusters")){
  Idents(obj.combined)<-i
  DimPlot(obj.combined, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".prefilter.umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt",
            "reads","molecules","genes_detected","percent_mt_molecules")){
  FeaturePlot(obj.combined, features= i,split.by = "condition",pt.size=0.1, order=T)
  ggsave2(paste0(i,".prefilter.umap.png"),width=10, height=5,device="png")
}

#Cluster DE
DefaultAssay(obj.combined) <- "RNA"
obj.combined<- NormalizeData(obj.combined) #if SCTransform data
obj.markers <- FindAllMarkers(object = obj.combined, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25, test.use = "MAST")
top8 <- obj.markers %>% group_by(cluster) %>% top_n(8, avg_logFC)
top4 <- obj.markers %>% group_by(cluster) %>% top_n(4, avg_logFC)
top2 <- obj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top1 <- obj.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
DefaultAssay(obj.combined) <- "integrated"
DoHeatmap(object = obj.combined, features = top8$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=10, height=10,device="png")
DefaultAssay(obj.combined) <- "RNA"
Idents(obj.combined) <- "seurat_clusters"
gene.list<-c("CD8A","CD4","CXCR5","CD40LG","CD3E","NKG7","CD79A","CD27","CD19","CD38","FOXP3","IGHD",
             "CCR7","ITGAM","CD14","CD68","CXCR4","ITGAX","CR1","CXCL13")
p<-FeaturePlot(object = obj.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",ncol=4, combine=F, pt.size=0.001, order=T)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p,ncol=5)
ggsave2("umap.genes.png",width=8,height=8,device="png")
VlnPlot(obj.combined, features = gene.list, group.by = "seurat_clusters",pt.size = 0)
ggsave2("vlnplot.genes.cluster.png",width=8, height=4,device="png")
DotPlot(obj.combined, features = gene.list)+ RotatedAxis()
ggsave2("dotplot.genes.cluster.png",width=8, height=5,device="png")
obj.combined@meta.data$my.clusters <- obj.combined@meta.data$seurat_clusters  # Store cluster identities in object@meta.data$my.clusters
cluster.counts<-table(obj.combined@meta.data$my.clusters, obj.combined@meta.data$condition)

#set metadata and subsets
samples<-names(table(obj.combined$sample))
Idents(obj.combined) <- "condition"
SLE.combined<-subset(obj.combined, idents='SLE')
Control.combined<-subset(obj.combined, idents='Control')
metadata<-unique(obj.combined@meta.data[,c("sample","condition")])
Control.sampleID<-metadata[metadata$condition=="Control",]$sample
SLE.sampleID<-metadata[metadata$condition=="SLE",]$sample

##Condition DE
Idents(obj.combined) <- "condition" #setting idents to condition metadata
obj.combined.autoimmune.response <- FindMarkers(obj.combined, ident.1 = "SLE", ident.2 = "Control", 
                                                  min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
avg.obj.combined <- log1p(AverageExpression(obj.combined)$RNA)
avg.obj.combined$gene <- rownames(avg.obj.combined)
top30 <- head(obj.combined.autoimmune.response, n = 30)
obj.combined.autoimmune.response$log2FC<-log2(exp(obj.combined.autoimmune.response$avg_logFC))
obj.combined$celltype.condition <- paste(obj.combined$my.clusters, obj.combined$condition, sep = "_") #adding metadata identifier

#condition DE (by sample)
Idents(obj.combined) <- "sample"
avg.obj.ms <- log1p(AverageExpression(obj.combined)$RNA)
avg.obj.ms$gene <- rownames(avg.obj.ms)
avg.obj.ms$Control.avg<-rowMeans(avg.obj.ms[,Control.sampleID])
avg.obj.ms$SLE.avg<-rowMeans(avg.obj.ms[,SLE.sampleID])
avg.obj.ms$log2FC<-log2(exp(avg.obj.ms$SLE.avg-avg.obj.ms$Control.avg))
for(i in 1:nrow(avg.obj.ms)) {
  avg.obj.ms$ttest[i] <- my.ttest(avg.obj.ms[i,Control.sampleID],avg.obj.ms[i,SLE.sampleID])
}

###Within Cluster DE
cluster.counts
for (i in 0:3){
  Idents(obj.combined) <- "celltype.condition" #setting idents to new metadata column
  df <- try(FindMarkers(obj.combined, ident.1 = paste0(i,"_SLE"), ident.2 = paste0(i,"_Control"), 
                    min.pct=0,logfc.threshold = -Inf,test.use = "MAST"))
  if(!is(df,"try-error")){
    df2 <- head(df, n = 30)
    assign(paste0("cluster",i,"top30"),df2)
    df$log2FC<-log2(exp(df$avg_logFC))
    assign(paste0("cluster",i,".autoimmune.response"),df)
    Idents(obj.combined) <- "my.clusters"
    temp<- subset(obj.combined, idents = i)
    Idents(temp) <- "condition"
    df <- log1p(AverageExpression(temp)$RNA)
    df$gene <- rownames(df)
    assign(paste0("avg.cluster",i),df)
    assign(paste0("cluster",i),temp)
  }
}

Idents(obj.combined) <- "my.clusters"
DefaultAssay(obj.combined) <- "SCT" 

#rename clusters
new.cluster.ids <- c("CD4","CD8","Mac","NK","B","B-ISG","Mono")
names(new.cluster.ids) <- levels(obj.combined)
obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
obj.combined@meta.data$my.clusters2 <- Idents(obj.combined)  

# ## tSNE/UMAP
# # Visualization
Idents(obj.combined)<-"my.clusters2"
DimPlot(obj.combined, reduction = "umap")
DimPlot(obj.combined, reduction = "umap",pt.size=0.001)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap.png",width=6, height=5,device="png")
ggsave2("umap.highres.tiff",width=7, height=6,dpi=300,device="tiff")
DimPlot(obj.combined, reduction = "umap", split.by = "condition")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap2.png",width=11, height=5,device="png")
#by sample and condition
Idents(obj.combined) <- "condition"
SLE.combined<-subset(obj.combined, idents='SLE')
Control.combined<-subset(obj.combined, idents='Control')
Idents(SLE.combined) <- "my.clusters2"
DimPlot(obj.combined, reduction = "tsne", split.by = "condition")
ggsave2("tsne2.png",width=11.5, height=5.5,device="png")
gene.list<-c("CD8A","CD4","CXCR5","CD40LG","CD3E","NKG7","CD79A","CD27","CD19","CD38","FOXP3","IGHD",
             "CCR7","ITGAM","CD14","CD68","CXCR4","ITGAX","CR1","CXCL13")
p<-FeaturePlot(object = obj.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F, pt.size=0.001, order=T)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist = p,ncol=5)
ggsave2("umap.GOI.genes.png",width=12,height=11,device="png")

## Cluster Analysis
Idents(obj.combined) <- "my.clusters2"
DefaultAssay(obj.combined) <- "integrated" 
obj.combined<-BuildClusterTree(obj.combined)
png("cluster.tree.png",width=4,height=6,units="in",res=300)
PlotClusterTree(obj.combined,font=1)
dev.off()
# # Cluster heatmap
top6 <- obj.markers %>% group_by(cluster) %>% top_n(6, avg_logFC)
DefaultAssay(obj.combined) <- "integrated"
DoHeatmap(object = obj.combined, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=6, height=6,device="png")
DefaultAssay(obj.combined) <- "RNA"
#UMAP by cluster marker
gene.list<-c("CD3E","CD4","CD8A","LYZ","KLRF1","MS4A1","BCL11A","IRF8")
p<-FeaturePlot(object = obj.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=4)
ggsave2("umap.clustermarkers.png",width=12,height=6,device="png")
#vln cluster marker
p<-VlnPlot(obj.combined, features = gene.list,group.by = "my.clusters2",pt.size = 0, combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=4)
ggsave2("vlnplot.clustermarkers.png",width=12, height=6,device="png")
#vln cluster marker
p<-RidgePlot(obj.combined, features = gene.list,group.by = "my.clusters2", combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=4)
ggsave2("ridgeplot.clustermarkers.png",width=15, height=6,device="png")
#dot plot
top2 <- obj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top2$gene<-make.unique(top2$gene,sep="--")
DotPlot(obj.combined, features = rev(top2$gene), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "condition")+ RotatedAxis()
top6$gene<-make.unique(top6$gene,sep="--")
DotPlot(obj.combined, features = top6$gene)+ 
  coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())+ RotatedAxis()
ggsave2("dotplot.png",width=6, height=12,device="png")

# Cluster frequency comparison
cluster.Counts <- table(obj.combined@meta.data$my.clusters2,obj.combined@meta.data$condition)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  theme(legend.title = element_blank(),axis.title=element_blank())+ 
  scale_x_discrete(labels= c("Control","SLE"))#+ theme(legend.position = "none")
ggsave2("cluster.prop.flow.png",width=4, height=4,device="png")
#cluster freq by sample
samples<-names(table(obj.combined@meta.data$sample))
freqlist=list()
for (i in samples) {
  freq<-as.data.frame(table(obj.combined@meta.data[obj.combined@meta.data$sample==i,]$my.clusters2, 
                            obj.combined@meta.data[obj.combined@meta.data$sample==i,]$condition))
  freq$pct<-100*(freq$Freq/sum(freq$Freq))
  freq$sample<-i
  freqlist[[i]]<-freq
}
cluster.freq.samples<-do.call(rbind,freqlist)
setnames(cluster.freq.samples, old=c("Var1","Var2"), new=c("cluster","condition"))
cluster.freq.samples$condition<-as.character(cluster.freq.samples$condition)
cluster.freq.samples$condition<-as.factor(cluster.freq.samples$condition)
ggbarplot(cluster.freq.samples, x = "cluster", y = "pct",add = c("mean_se", "jitter"),palette=c("black","red"),
          color = "condition",xlab=F,position = position_dodge(0.8), legend="right")+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2,label.y=110)+#stat_compare_means(method="anova",label.y=0)+
  labs(y = "Frequency")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
ggsave2("cluster.freq.samples.png",width=4.5, height=3,device="png")

#UMAP comparison
gene.list<-c("PRG4","CXCR4","MX1","TIMP3")
p<-FeaturePlot(obj.combined, features = gene.list, split.by = "condition",
               pt.size=0.001,combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+NoAxes()+
    theme(panel.border = element_rect(colour = "black", size=1),
          plot.title=element_blank(),axis.title.y.right=element_blank(),
          axis.line=element_blank())
}
cowplot::plot_grid(p[[1]],p[[5]],p[[2]],p[[6]],p[[3]],p[[7]],p[[4]],p[[8]],ncol=2)
ggsave2("umap.GOI.bycondition.png",width=6, height=10,device="png")

#Population VlnPlot
DefaultAssay(obj.combined)<-"RNA"
p<-VlnPlot(object = obj.combined, features =gene.list, pt.size=0.001,group.by ="condition",cols=c("grey","red"),combine=F)
for(i in 1:(length(p))) {
  p[[i]] <- p[[i]]+theme(axis.title=element_blank(),axis.text.x=element_blank(),
                         plot.title = element_text(size=30))+NoLegend()
}
cowplot::plot_grid(plotlist = p,ncol=2)
ggsave2("vlnpopulation.GOI.bycondition.png",width=7, height=10,device="png")
p<-VlnPlot(obj.combined, features = gene.list, split.by = "condition", 
           group.by = "my.clusters2",cols=c("grey","red"),pt.size = 0,combine=F)
for(i in 1:(length(p)-1)) {
  p[[i]] <- p[[i]]+theme(axis.title=element_blank(),axis.text.x=element_blank())# NoLegend()+
}
p[[length(p)]] <- p[[length(p)]]+theme(axis.title=element_blank())
cowplot::plot_grid(plotlist = p,ncol=1,rel_heights=c(1,1,1,1.5))
ggsave2("vlnplot.GOI.cluster.png",width=3, height=8,device="png")
DefaultAssay(obj.combined) <- "integrated"
Idents(obj.combined)<-"condition"
DoHeatmap(object = obj.combined, features = rownames(top30), label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("SLEvsControlheatmap.png",width=4, height=6,device="png")
DefaultAssay(obj.combined) <- "RNA"

## Within Cluster Comparison
for(i in ls(pattern="avg")){
  df<-get(i)
  if(!is(df,"try-error")){
    colnames(df)[apply(sapply(c("Control","log2FC.y"), function (y) sapply(colnames(df), 
                                                                       function (x) grepl(y, x))), 1, any)]<-"Control"
    colnames(df)[apply(sapply(c("SLE","log2FC.x"), function (y) sapply(colnames(df), 
                                                                        function (x) grepl(y, x))), 1, any)]<-"SLE"
    df$sig<-"unsig"
    df$sig[abs(df$SLE-df$Control)>0.2]<-"DE"
    gene.list<-rownames(df[abs(df$SLE-df$Control)>0.2&!is.na(df$gene),])
    if(length(gene.list)>30){gene.list<-rownames(df[abs(df$SLE-df$Control)>0.5&!is.na(df$gene),])}
    if(length(gene.list)>30){gene.list<-rownames(df[abs(df$SLE-df$Control)>1&!is.na(df$gene),])}
    df$sig <- factor(df$sig, levels = c("unsig","DE"))
    p2 <- ggplot(df, aes(SLE,Control)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) + 
      ggtitle(i)+theme_classic()+
      labs(x="SLE",y="Control")+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
            axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
      scale_color_manual(values = c("black","red"))
    if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
    p2
    ggsave2(paste0(i,".scatter.png"),width=4, height=4,device="png")
  }
}
#KOvsWT cluster UMAP
cluster.counts
clusters<-0:3
for(i in clusters){
  df<-get(paste0("cluster",i,".autoimmune.response"))
  gene.list<-rownames(head(df, n = 3))
  FeaturePlot(obj.combined, features = gene.list, split.by = "condition", max.cutoff = 3,cols = c("grey", "red"))
  ggsave2(paste0("umap.cluster",i,".DEgenes.png"),width=10, height=14,device="png")
}
#KOvsWT cluster VlnPlot
for(i in clusters){
  df<-get(paste0("cluster",i,".autoimmune.response"))
  gene.list<-rownames(head(df, n = 6))
  cluster.seurat<-get(paste0("cluster",i))
  DefaultAssay(cluster.seurat)<-"RNA"
  VlnPlot(object = cluster.seurat, pt.size=0.02, features =gene.list, ncol=2)
  ggsave2(paste0("cluster",i,".DEgenes.vln.png"),width=4, height=7,device="png")
}
# KOvsWT cluster heatmap
for(i in clusters){
  df<-get(paste0("cluster",i,".autoimmune.response"))
  gene.list<-rownames(head(df, n = 30))
  cluster.seurat<-get(paste0("cluster",i))
  DefaultAssay(cluster.seurat) <- "integrated"
  DoHeatmap(object = cluster.seurat, features = gene.list, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
  ggsave2(paste0("564vsControlheatmap.cluster",i,".png"),width=6, height=6,device="png")
  DefaultAssay(cluster.seurat) <- "RNA"
}
#Volcano Plots
for(i in c(ls(pattern=".autoimmune.response"))){
  df<-get(i)
  if(is.numeric(df$avg_logFC)){
    df$log2FC<-log2(exp(df$avg_logFC))
    gene.list<-rownames(df[abs(df$log2FC)>0.5&df$p_val_adj<0.001,])
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>1&df$p_val_adj<10e-10,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>1.25&df$p_val_adj<10e-15,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>1.5&df$p_val_adj<10e-20,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>1.75&df$p_val_adj<10e-25,])}
    if(length(gene.list)>20){gene.list<-rownames(df[abs(df$log2FC)>2&df$p_val_adj<10e-30,])}
    if(length(gene.list)>30){gene.list<-NA}
    EnhancedVolcano(df,lab = rownames(df),
                    x = 'log2FC',y = 'p_val_adj',title = i,col=c("black","black","black","red3"),
                    selectLab=gene.list,xlab=bquote(~Log[2]~ (frac("SLE","Control"))),
                    pCutoff = 0.001,FCcutoff = 0.5,pointSize = 0.5,labSize = 2,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                    legendPosition="none",drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                    subtitle="", caption="",border="full",cutoffLineWidth=0,
                    gridlines.major=F,gridlines.minor=F,titleLabSize=10
    )
    ggsave2(paste0(i,".volcano.png"),width=3.5, height=4,device="png")
  }
}


###MS compare

brenner.all<-obj.combined.autoimmune.response
brenner.CD4<-cluster0.autoimmune.response
brenner.CD8<-cluster1.autoimmune.response
brenner.mac<-cluster2.autoimmune.response

#load mouse SLE data
load("../graphed.RData")
mouse.all<-SLE.obj.combined.autoimmune.response
mouse.CD4<-cluster1.autoimmune.response
mouse.CD8<-cluster3.autoimmune.response
mouse.mac<-cluster2.autoimmune.response


for(k in c("CD4","CD8","mac","all")){
  brenner.raw<-get(paste0("brenner.",k))
  brenner.raw$hgnc_symbol<-rownames(brenner.raw)
  #import mouse gene symbols
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mgi_list<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = brenner.raw$hgnc_symbol , 
                   mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  m <- match(brenner.raw$hgnc_symbol, mgi_list$HGNC.symbol)
  brenner.raw<-cbind(brenner.raw,mgi_list[m,])
  brenner.raw[brenner.raw$hgnc_symbol %like% ".12706", ]
  #mouse vs human
  mouse.raw<-get(paste0("mouse.",k))
  mouse.raw$MGI.symbol<-rownames(mouse.raw)
  mouse<-data.frame(gene=rownames(mouse.raw),ms.log2FC=mouse.raw$log2FC)
  brenner<-data.frame(gene=brenner.raw$MGI.symbol,h.log2FC=brenner.raw$log2FC,h.gene=brenner.raw$hgnc_symbol)
  
  #venn diagram
  DE.list<-list()
  for(i in c("brenner.raw","mouse.raw")){
    df<-get(i)
    DE<-df[abs(df$log2FC)>0.2&df$p_val_adj<0.05,]
    DE.list[[i]]<-DE$MGI.symbol[!is.na(DE$MGI.symbol)]
  }
  names(DE.list)<-c("Arazi et al.","SLE.yaa")
  venn.diagram(
    x = DE.list,category.names = names(DE.list),filename = paste0(k,".human.ms.venn.png"),
    imagetype="png",height=700,width=700,margin=0.1,
    lwd=1,col=c("darkmagenta", "chartreuse4"),#,"cyan4"
    fill=c(alpha("darkmagenta",0.3), alpha('chartreuse4',0.3)),#, alpha('cyan4',0.3)
    cex=0.3,fontface="bold",fontfamily="sans",ext.line.lwd=0.2,ext.length=0.9,ext.dist=0.05,
    cat.cex=0.3,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
    cat.col = c("darkmagenta", "chartreuse4"),cat.dist = c(0.1, 0.1)#,cat.pos=c(0,180,0)
  )
  intersect(DE.list[[1]],DE.list[[2]]) 
  #correlation
  df<-brenner
  comp<-merge(mouse,df,by="gene")
  comp<-comp[!is.na(comp$h.log2FC),]
  rownames(comp)<-make.unique(as.character(comp$gene))
  comp$sig<-"unsig"
  comp$sig[abs(comp$ms.log2FC)>0.2&abs(comp$h.log2FC)>0.2]<-"co.DE"
  comp$sig[abs(comp$ms.log2FC)>0.2&abs(comp$h.log2FC)<=0.5]<-"ms.DE"
  comp$sig[abs(comp$ms.log2FC)<=0.2&abs(comp$h.log2FC)>0.5]<-"h.DE"
  ms.up<-head(rownames(comp[order(abs(comp$ms.log2FC),decreasing=T),]),n=20)
  h.up<-head(rownames(comp[order(abs(comp$h.log2FC),decreasing=T),]),n=20)
  diff<-rownames(comp[abs(comp$ms.log2FC-comp$h.log2FC)>0.5&!is.na(comp$gene),])
  diff<-head(rownames(comp[order(abs(comp$ms.log2FC-comp$h.log2FC),decreasing=T),]),n=20)
  gene.list<-unique(c(ms.up,h.up,diff))
  comp$sig <- factor(comp$sig, levels = c("unsig","co.DE","ms.DE","h.DE"))
  p2 <- ggplot(comp, aes(ms.log2FC,h.log2FC)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) + 
    labs(x=bquote(~Log[2]~ (frac("SLE.yaa","B6"))),y=bquote(~Log[2]~ (frac("SLE","HC"))))+theme_classic()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
          axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
    geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
    scale_color_manual(values = c("grey","black","red","orange"))
  if(length(gene.list)>0 & length(gene.list)<=60){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
  p2
  ggsave2(paste0(k,".human.ms.scatter.png"),width=5, height=5,device="png")
}

