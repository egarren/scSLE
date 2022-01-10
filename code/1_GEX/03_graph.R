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
library(ggalluvial)

load("analyzed.RData")

## change depending on which subset analyzing
# SLE.obj.combined<-SLE.myeloid.obj.combined
# SLE.obj.markers<-SLE.myeloid.obj.markers
# SLE.obj.combined<-SLE.T.obj.combined
# SLE.obj.markers<-SLE.T.obj.markers
# SLE.obj.combined<-SLE.B.obj.combined
# SLE.obj.markers<-SLE.B.obj.markers
# SLE.obj.combined<-SLE.CD4.obj.combined
# SLE.obj.markers<-SLE.CD4.obj.markers
# SLE.obj.combined<-SLE.CD8.obj.combined
# SLE.obj.markers<-SLE.CD8.obj.markers

Idents(SLE.obj.combined) <- "my.clusters"
DefaultAssay(SLE.obj.combined) <- "SCT" 

#rename clusters
# gene.list<-c("Cd3e","Cd4","Cd8a","Itgam","Cd14","Itgae","Sirpa","Siglech","Ly6c1","Ccr2","Cx3cr1",
#              "Fcgr1","Ptprc","Cd19","Mki67","Tbx21","Klra1","Bcl6","Cxcr5","Ly6g","Cd68","Itgax",
#              "Adgre1","Mertk","Xcr1","Clec9a","Csf1r","Lyz2","Timd4","Siglec1","Sell","Sdc1","Prdm1","Cr2",
#              "Fcer2a","Cd1d1","Ighd","Ighg1","Ighm","Itgax","Fcrl5","Siglecf","Clec9a","Gzmb","Aicda")
# gene.list<-c("Cd14","Sirpa","Siglech","Ccr2","Cx3cr1","S100a8",
#              "Fcgr1","Ly6g","Cd68","Itgax","Cxcr2","Arg2",
#              "Adgre1","Mertk","Csf3r","Clec9a","Csf1r","Lyz2","Siglec1","Cr2",
#              "Fcer2a","Cd1d1","Fcrl5","Siglecf","Clec9a","Timd4","Vcan","Itgal","Ace",
#              "Il1b","Ccr5","C1qb","Fabp5","Saa3","Clec10a","Flt3",
#              "F13a1","Itgal","S100a9","Cd74")
# gene.list<-c("Ptprc","Cd19","Mki67","Tbx21","Ighd","Ighg1","Ighm","Fcrl5",
#              "Ms4a1","Cd38","Cd27","Sdc1","Cr2","Xbp1","Slamf7","Cd86","Bach2","Cxcr5",
#              "Bcl6","Aicda","Cd83","Cxcr4","Sema7a","Prdm1",
#              "Cd1d1","Cd24a","Cd93","Cd5","Spn","Fcer2a")
# gene.list<-c("Cd3e","Cd4","Cd8a","Cd40lg","Cxcr5","Bcl6","Ccl5","Foxp3",
#              "Ctla4","Pdcd1","Sell","Icos","Cxcr4","Ccr6","Ccr7")
# gene.list<-c("Cd4","Cd40lg","Cxcr5","Foxp3",
#              "Ctla4","Pdcd1","Sell","Icos","Ccr7")
# gene.list<-c("Cd8a","Ccr7","Sell","Id3","Gzmk","Tcf7","Nr4a1","Il2rb","Klra3")
# gene.list<-c("Itgax","Itgam","Marco","Cd209a")
p<-FeaturePlot(object = SLE.obj.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F, pt.size=0.001, order=T)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend()+NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist = p,ncol=3)
ggsave2("umap.GOI.genes.png",width=6,height=6.5,device="png")
##change depending on subset
# new.cluster.ids <- c("B","CD4","Mac","CD8","PMN","B-early","DC","PC","NK","B-ISG","pDC","Baso")
# new.cluster.ids<-c("Mac-1","Mono","PMN-1","Mac-2","PMN-2","Mac-3","FDC","pDC","Baso")
# new.cluster.ids<-c("Fo","Activated","MZ","T2-B","PC","T1-B","B-ISG","B1")
# new.cluster.ids<-c("CD4","CD8","CD8-eff","Treg","Cd74","CD4-eff")
# new.cluster.ids<-c("CD4","Treg","CD4-eff","CD4-mem","CD4-exhausted","CD4-naive","CD4-activated","Tfh")
# new.cluster.ids<-c("CD8-naive","CD8-eff","CD8-stem")
names(new.cluster.ids) <- levels(SLE.obj.combined)
SLE.obj.combined <- RenameIdents(SLE.obj.combined, new.cluster.ids)
SLE.obj.combined@meta.data$my.clusters2 <- Idents(SLE.obj.combined)  

# ## tSNE/UMAP
# # Visualization
DimPlot(SLE.obj.combined, reduction = "umap", label=T)
DimPlot(SLE.obj.combined, reduction = "umap",pt.size=0.001, label=T)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                   axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap.png",width=8, height=7,device="png")
ggsave2("umap.highres.tiff",width=8, height=5,dpi=300,device="tiff")
DimPlot(SLE.obj.combined, reduction = "umap", split.by = "condition")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2("umap2.png",width=11, height=5,device="png")
#by mouse and condition
Idents(SLE.obj.combined) <- "condition"
SLE.yaa.combined<-subset(SLE.obj.combined, idents='SLE.yaa')
B6.combined<-subset(SLE.obj.combined, idents='B6')
Idents(SLE.yaa.combined) <- "my.clusters"
DimPlot(SLE.yaa.combined, reduction = "umap", split.by = "mouse_ID")+ NoLegend()+NoAxes()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("umap.SLE.yaa.mouse.png",width=5.5, height=3,device="png")
Idents(B6.combined) <- "my.clusters"
DimPlot(B6.combined, reduction = "umap", split.by = "mouse_ID")+ NoLegend()+NoAxes()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("umap.B6.mouse.png",width=5.5, height=3,device="png")


## Cluster Analysis
Idents(SLE.obj.combined) <- "my.clusters2"
DefaultAssay(SLE.obj.combined) <- "integrated" 
SLE.obj.combined<-BuildClusterTree(SLE.obj.combined)
png("cluster.tree.png",width=4,height=6,units="in",res=300)
PlotClusterTree(SLE.obj.combined,font=1)
dev.off()
# # Cluster heatmap
top6 <- SLE.obj.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
DefaultAssay(SLE.obj.combined) <- "integrated"
DoHeatmap(object = SLE.obj.combined, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap2.png",width=15, height=10,device="png")
DefaultAssay(SLE.obj.combined) <- "RNA"
#UMAP by cluster marker
DefaultAssay(SLE.obj.combined) <- "RNA" 
i="Ciita"
FeaturePlot(object = SLE.obj.combined, features = i, cols = c("grey", "blue"), split.by="condition",reduction = "umap")
ggsave2(paste0("umap.",i,".png"),width=10,height=5,device="png")
VlnPlot(SLE.obj.combined, features = i,group.by = "my.clusters2",split.by="condition",pt.size = 0.1)
ggsave2(paste0("vln.",i,".png"),width=10,height=5,device="png")
# gene.list<-c("Cd19","Cd4","Cd8a","Cd68","Aicda","S100a8","Cr2","Itgax","Mki67","C1qb","Klra8","Ifit3","Itgam","Siglech","Cd63")
# gene.list<-c("Cd19","Cd4","Cd8a","Cd68","Prdm1","S100a8","Vpreb3","Vcam1","Jchain","C1qb","Klra8","Ifit3","Irf7","Siglech","Cd63")
# gene.list<-c("Cd68","S100a8","C1qb","Siglech","Cd63","Ace","Ccr2","Csf3r","Fcrl5","Ltf","Cd5l","Cd21","Ly6g")
# gene.list<-c("Cd19","Mki67","Ighd","Fcrl5","Sdc1","Cr2","Xbp1","Slamf7","Aicda","Cd83","Sema7a","Prdm1",
#              "Cd1d1","Cd24a","Cd93","Cd5","Spn","Fcer2a","Lars2","Mif","Jchain","Vpreb3","Ifit3","Apoe")
# gene.list<-c("Mki67","Cr2","Jchain","Ms4a1","Itgax","Irf7","Cd93","Apoe","Cd93","Xbp1","Cd24a","Fcrl5")
# gene.list<-c("Cd3e","Cd4","Cd8a","Cd40lg","Ccl5","Nkg7","Foxp3","Cd74","Cd79a","Cxcr6")
# gene.list<-c("Cd4","Cd40lg","Ccl5","Foxp3","Igfbp4","Lef1","Stat1","Lars2","Id2","Ccr7","Nr4a1","Pdcd1")
# gene.list<-c("Cd8a","Ccl5","Igfbp4","Lef1","Ccr7","S100a6","Myb","Slamf6","Id3","Il2rb")
p<-FeaturePlot(object = SLE.obj.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=4)
ggsave2("umap.clustermarkers.png",width=18,height=11,device="png")
ggsave2("umap.clustermarkers2.png",width=9,height=7,device="png")
gene.list<-c("Itgax","Fcrl5")
p<-FeaturePlot(object = SLE.obj.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=2)
ggsave2("umap.mz.png",width=8,height=4,device="png")
#vln cluster marker
gene.list<-c("Itgax","Fcrl5")
p<-VlnPlot(SLE.obj.combined, features = gene.list,group.by = "my.clusters2",pt.size = 0.1, combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=3)
ggsave2("vlnplot.clustermarkers.png",width=8, height=4,device="png")
ggsave2("vlnplot.goi.png",width=8, height=4,device="png")
#dot plot
top6$gene<-make.unique(top6$gene,sep="--")
DotPlot(SLE.obj.combined, features = top6$gene)+ 
  coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())+ RotatedAxis()
ggsave2("dotplot.png",width=6, height=12,device="png")
DotPlot(SLE.obj.combined, features = top6$gene)+ 
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title=element_blank())
ggsave2("dotplot2.png",width=12, height=5,device="png")

# Cluster frequency comparison
ggplot(SLE.obj.combined@meta.data, aes(x = my.clusters2,fill = condition)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)),position='dodge')+ 
  scale_y_continuous(labels = percent)+ 
  labs(y="Frequency", x = "Cluster")
ggsave2("clusterfreq.png",width=4, height=2.5,device="png")
cluster.Counts <- table(SLE.obj.combined@meta.data$my.clusters2,SLE.obj.combined@meta.data$condition)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  theme(legend.title = element_blank(),axis.title=element_blank())+ 
  scale_x_discrete(labels= c("B6","SLE.yaa"))#+ theme(legend.position = "none")
ggsave2("cluster.prop.flow.png",width=4, height=4,device="png")
#cluster freq by mouse
mice<-dimnames(table(SLE.obj.combined@meta.data$mouse_ID))[[1]]
freqlist=list()
for (i in mice) {
  freq<-as.data.frame(table(SLE.obj.combined@meta.data[SLE.obj.combined@meta.data$mouse_ID==i,]$my.clusters2, 
                            SLE.obj.combined@meta.data[SLE.obj.combined@meta.data$mouse_ID==i,]$condition))
  freq$pct<-100*(freq$Freq/sum(freq$Freq))
  freq$mouse<-i
  freqlist[[i]]<-freq
}
cluster.freq.mice<-do.call(rbind,freqlist)
setnames(cluster.freq.mice, old=c("Var1","Var2"), new=c("cluster","condition"))
cluster.freq.mice$condition<-as.character(cluster.freq.mice$condition)
cluster.freq.mice$condition[cluster.freq.mice$condition%in% c("B6")]<-"B6"
cluster.freq.mice$condition[cluster.freq.mice$condition%in% c("SLE.yaa")]<-"SLE.yaa"
cluster.freq.mice$condition<-as.factor(cluster.freq.mice$condition)
cluster.freq.mice$condition <- factor(cluster.freq.mice$condition, levels = c("B6","SLE.yaa"))
ggbarplot(cluster.freq.mice, x = "cluster", y = "pct",add = c("mean_se", "jitter"),palette=c("black","red"),
          color = "condition",xlab=F,position = position_dodge(0.8), legend="right")+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2,label.y=60)+#stat_compare_means(method="anova",label.y=0)+
  labs(y = "Frequency")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
ggsave2("cluster.freq.mice.png",width=4.5, height=3,device="png")

#Volcano Plots
df<-SLE.obj.combined.autoimmune.response
df<-SLE.myeloid.obj.combined.autoimmune.response
df<-SLE.B.obj.combined.autoimmune.response
df<-SLE.T.obj.combined.autoimmune.response
df<-SLE.CD4.obj.combined.autoimmune.response
df<-SLE.CD8.obj.combined.autoimmune.response

df<-df[!grepl("Igh|Igkc",rownames(df)),]
df2<-df[df$p_val_adj<0.01,]
gene.list<-head(rownames(df2[order(df2$p_val_adj),]),n=50)
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(df,lab = rownames(df),
                x = 'log2FC',y = 'p_val_adj',title = "",col=c("black","black","black","red3"),
                selectLab=gene.list,xlab=bquote(~Log[2]~ (frac("SLE.yaa","B6"))),
                pCutoff = 0.01,FCcutoff = 0.2,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                legendPosition="none",drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                subtitle="", caption="",border="full",cutoffLineWidth=0,
                gridlines.major=F,gridlines.minor=F,titleLabSize=10
)
ggsave2("SLE.obj.combined.autoimmune.response.volcano2.png",width=7, height=6,device="png")
for(i in c(ls(pattern=".autoimmune.response"))){
  df<-get(i)
  if(is.numeric(df$avg_log2FC)){
    df<-df[!grepl("Igh|Igkc",rownames(df)),]
    df2<-df[df$p_val_adj<0.01,]
    gene.list<-head(rownames(df2[order(df2$p_val_adj),]),n=20)
    EnhancedVolcano(df,lab = rownames(df),
                    x = 'log2FC',y = 'p_val_adj',title = i,col=c("black","black","black","red3"),
                    selectLab=gene.list,xlab=bquote(~Log[2]~ (frac("SLE.yaa","B6"))),
                    pCutoff = 0.01,FCcutoff = 0.2,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                    legendPosition="none",#drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                    subtitle="", caption="",border="full",cutoffLineWidth=0,
                    gridlines.major=F,gridlines.minor=F,titleLabSize=10
    )
    ggsave2(paste0(i,".volcano.png"),width=4, height=4,device="png")
  }
}

#UMAP comparison
goi="Cd74"
FeaturePlot(SLE.obj.combined, features = goi, split.by = "condition", max.cutoff = 3,cols = c("grey", "red"), pt.size=0.01, order=T)
ggsave2(paste0("umap.",goi,".png"),width=10, height=5,device="png")
gene.list<-c("Itgax","Itgam","Marco","Cd209a")
FeaturePlot(SLE.obj.combined, features = gene.list, split.by = "condition", max.cutoff = 3,cols = c("grey", "red"), pt.size=0.01, order=T)
ggsave2("umap2.goi.png",width=8, height=12,device="png")
FeaturePlot(SLE.obj.combined, features = gene.list, split.by = "mouse_ID", max.cutoff = 3,cols = c("grey", "red"))
ggsave2("umap.GOI.bymouse.png",width=19, height=7,device="png")
gene.list<-c("Gm42031","Ly6a","Lag3","Id3")
p<-FeaturePlot(SLE.obj.combined, features = gene.list, split.by = "condition",
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
DefaultAssay(SLE.obj.combined)<-"RNA"
gene.list<-c("B2m","Cebpb","Cr2","Ctsb","Cybb","Fcer1g","Gas5","Gngt2","H2.K1","Itgb1","Lgals3","Ifitm2")
gene.list<-c("Fabp4","Gngt2","Cd74","Pecam1","Cd300e","Serpinb6a","Fcgr3","Gpx1","ccl6","Pltp","Gfl2") #myeloid
gene.list<-c("AW112010","Lag3","Ly6a","Rps24","S100a11","S100a6","Tigit","Maf","Itgb1","Tnfrsf4","Cxcr3","Satb1","Lef1","Eea1","Igfbp4") #T
gene.list<-c("B2m","Cxcr3","H2.K1","Jchain","Ly6a","Stat1","Tmsb4x","Cr2","Tap1","Cd24a","Slpi","Slamf6","Siglecg","Socs3","Ly6c2") #B
gene.list<-c("S100a6","Lag3","Ly6a","Tigit","Srgn","Eea1","Maf","Itgb1","Tnfrsf4","Satb1","Lef1","Igfbp4","Actn1","Rps24","Cd82") #CD4
gene.list<-c("Lef1","Ly6a","AW112010","Gzmk","Ccl5","Lag3","Ptms","Nkg7","H2.Q7") #CD8
p<-VlnPlot(object = SLE.obj.combined, features =gene.list, pt.size=0.001,group.by ="condition",cols=c("grey","red"),combine=F)
for(i in 1:(length(p))) {
  p[[i]] <- p[[i]]+theme(axis.title=element_blank(),axis.text.x=element_blank(),
                         plot.title = element_text(size=30))+NoLegend()
}
cowplot::plot_grid(plotlist = p,ncol=3)
ggsave2("vlnpopulation.GOI.bycondition.png",width=15, height=8,device="png")

##Ig Isotype analysis
DefaultAssay(SLE.obj.combined)<-"RNA"
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol","external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "mgi_symbol",values=list(rownames(SLE.obj.combined@assays[["RNA"]]@meta.features)),
                  mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="ensembl.org"),useCache=F) #useast.
gene.list<-unique(genes.meta$mgi_symbol[genes.meta$gene_biotype=="IG_C_gene"])
p<-VlnPlot(object = SLE.obj.combined, features =gene.list, pt.size=0.001,group.by ="condition",cols=c("grey","red"),combine=F)
for(i in 1:(length(p))) {
  p[[i]] <- p[[i]]+theme(axis.title=element_blank(),axis.text.x=element_blank())+NoLegend()
}
cowplot::plot_grid(plotlist = p,ncol=7)
ggsave2("vlnpopulation.Ig.isotype.bycondition.png",width=12, height=5,device="png")


#Cluster VlnPlot
goi="Fcgr4"
VlnPlot(SLE.obj.combined, features = goi, split.by = "condition", cols=c("grey","red"),group.by = "my.clusters2",pt.size = 0.5)#,legend='right')
ggsave2(paste0("vlnplot.",goi,".cluster.png"),width=8, height=5,device="png")
gene.list<-c("Cd74","Pecam1","Fcgr3","H2.Aa","C1qb","Ciita")
p<-VlnPlot(SLE.obj.combined, features = gene.list, split.by = "condition", 
        group.by = "my.clusters2",cols=c("grey","red"),pt.size = 0.1,combine=F)
for(i in 1:(length(p))) { #-1
  p[[i]] <- p[[i]]+theme(axis.title=element_blank())+NoLegend() #,axis.text.x=element_blank()
}
cowplot::plot_grid(plotlist = p,ncol=3)#,rel_heights=c(1,1,1,1.5))
ggsave2("vlnplot.GOI.cluster.png",width=13, height=4.5,device="png")
gene.list<-c("Itgax","Fcrl5")
p<-VlnPlot(SLE.obj.combined, features = gene.list, split.by = "condition", 
           group.by = "my.clusters2",cols=c("grey","red"),pt.size = 0.1,combine=F)
for(i in 1:(length(p))) { #-1
  p[[i]] <- p[[i]]+theme(axis.title=element_blank())+NoLegend() #,axis.text.x=element_blank()
}
cowplot::plot_grid(plotlist = p,ncol=2)#,rel_heights=c(1,1,1,1.5))
ggsave2("vlnplot.mz.cluster.png",width=8, height=4,device="png")
VlnPlot(SLE.obj.combined, features = "Tbx21", split.by = "condition", 
        group.by = "my.clusters2",cols=c("grey","red"),pt.size = 0.1,combine=F)
ggsave2("vlnplot.tbx21.cluster.png",width=8, height=4,device="png")
# SLE.yaavsB6 heatmap
Idents(SLE.obj.combined)<-"condition"
DefaultAssay(SLE.obj.combined) <- "integrated"
DoHeatmap(object = SLE.obj.combined, features = rownames(top30), label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("SLE.yaavsB6heatmap.png",width=4, height=6,device="png")
DefaultAssay(SLE.obj.combined) <- "RNA"


save.image("graphed.RData")
save(list=c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"),"clusters"),file="DE.dfs.seurat.RData")
saveRDS(SLE.obj.combined, file="seurat_obj.rds")
