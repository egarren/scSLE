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
library(ggseqlogo)
library(immunarch)
library(cowplot)
library(ggrepel)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}
data_concater2 <- function(x){
  x<- levels(factor(x))
  paste(x[1])
}
data_concater3 <- function(x){
  x<- levels(factor(x))
  if(length(x)>1){paste(x[1:2], collapse = "+")}else{paste(x[1])}
}

load("clone.ab.data.RData")
load("graphed.RData")

T_all.combined<-SLE.obj.combined
#paired data
ab.barcode<-Tall.vdj[Tall.vdj$productive=="true"&!is.na(Tall.vdj$productive)&Tall.vdj$chain %in% c("TRA","TRB"),] 
ab.barcode<-as.data.table(ab.barcode)[, lapply(.SD, data_concater3), by=barcode]
ab.barcode<-ab.barcode[ab.barcode$chain=="TRA+TRB",]
ab.barcode<-transform(ab.barcode, cdr3.freq = ave(seq(nrow(ab.barcode)), cdr3, FUN=length))
ab.barcode<-transform(ab.barcode, ms_cdr3.freq = ave(seq(nrow(ab.barcode)), ms_cdr3, FUN=length))
ab.barcode<-transform(ab.barcode, vab.freq = ave(seq(nrow(ab.barcode)), v_gene, FUN=length))
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
rownames(ab.barcode)<-ab.barcode$barcode
T_all.combined<-AddMetaData(T_all.combined, metadata=ab.barcode)

#key objects
Idents(T_all.combined) <- "condition"
SLE.yaa.combined<-subset(T_all.combined, idents='SLE.yaa')
B6.combined<-subset(T_all.combined, idents='B6')
# identifying top clones
SLE.yaa.top.clones<-names(head(sort(table(SLE.yaa.combined@meta.data$cdr3), decreasing = T), n=9))
B6.top.clones<-names(head(sort(table(B6.combined@meta.data$cdr3), decreasing = T), n=9))
save.image("vdj.RData")

#mapping TCR recovery
DimPlot(T_all.combined, reduction = "umap", group.by="chain",split.by = "condition",
        cols=c("#e31a1c","#08519c","#33a02c"))+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
ggsave2("TRmap.png",width=9, height=4,device="png")

##Clonotype recovery
total.recoveryfreq=1-sum(is.na(T_all.combined@meta.data$chain))/nrow(T_all.combined@meta.data)
SLE.yaa.recoveryfreq=1-sum(is.na(SLE.yaa.combined@meta.data$chain))/nrow(SLE.yaa.combined@meta.data)
B6.recoveryfreq=1-sum(is.na(B6.combined@meta.data$chain))/nrow(B6.combined@meta.data)
#recovery by mouse
recovery.df<-as.data.frame(table(T_all.combined@meta.data$mouse_ID))
recovery.df$recovery<-0
recovery.df$condition<-0
for (i in 1:nrow(recovery.df)) {
  mouse<-recovery.df$Var1[i]
  df<-T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==mouse,]
  recovery.df$recovery[i]<-100*(1-sum(is.na(df$chain))/nrow(df))
  recovery.df$condition[i]<-T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==mouse,]$condition[1]
}
ggbarplot(recovery.df, x = "condition", y = "recovery", 
          add = c("mean_se", "jitter"),
          position = position_dodge(0.8),palette="npg")+
  stat_compare_means(aes(group = condition),label = "p.signif", method="t.test")+
  labs(x="condition", y = "% recovery")
ggsave2("TRrecovery.png",width=4, height=3,device="png")
#recovery by cluster
recovery.df<-as.data.frame(table(T_all.combined@meta.data$mouse_ID,T_all.combined@meta.data$my.clusters2))
recovery.df$recovery<-0
recovery.df$condition<-"B6"
for (i in 1:nrow(recovery.df)) {
  mouse<-recovery.df$Var1[i]
  cluster<-recovery.df$Var2[i]
  df<-T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==mouse&T_all.combined@meta.data$my.clusters2==cluster,]
  recovery.df$recovery[i]<-100*(1-sum(is.na(df$chain))/nrow(df))
  cond<-metadata$condition[metadata$mouse_ID==mouse]
  if(cond =="SLE.yaa"){recovery.df$condition[i]<-"SLE.yaa"}
  # recovery.df$condition[i]<-T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==mouse,]$condition[1]
}
recovery.df$condition <- factor(recovery.df$condition, levels = c("B6","SLE.yaa"))
ggbarplot(recovery.df, x = "Var2", y = "recovery",add = c("mean_se", "jitter"),
          color = "condition",position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2.5,label.y=110)+
  labs(y = "% recovery")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
ggsave2("TRrecovery.bycluster.png",width=5, height=3,device="png")

#mapping expanded clones by cdr3aa
Tall.vdj.clones<-subset(x = T_all.combined, subset = cdr3!="None")
FeaturePlot(Tall.vdj.clones, features= "cdr3.freq",pt.size=0.1, order=T,cols=c("grey","darkgreen"))+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
        axis.title=element_blank(),plot.title=element_blank(),legend.position="none")
ggsave2("umap.combined.cdr3.freq.png",width=3, height=3,device="png")
p<-FeaturePlot(Tall.vdj.clones, features= "cdr3.freq",split.by = "condition",pt.size=0.1, 
               order=T,cols=c("grey","darkgreen"),combine=F) #NoLegend()+
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] +NoAxes()+#NoLegend()+
    theme(panel.border = element_rect(colour = "black", size=1),plot.title=element_blank(),
          axis.title.y.right=element_blank(),
          axis.line=element_blank())
}
cowplot::plot_grid(p[[1]],p[[2]],ncol=2)
ggsave2("umap.cdr3.freq.png",width=9, height=3.5,device="png")

##Mapping individual clones within condition
plot.list <- list()
Idents(B6.combined) <- "cdr3"
# plot.list <- list()
for (i in unique(x = names(head(sort(table(B6.combined@meta.data$cdr3), decreasing = T), n=4)))) {
  plot.list[[i]] <- DimPlot(
    object = B6.combined, cols.highlight="forestgreen",
    cells.highlight = Cells(subset(B6.combined, idents=i)),sizes.highlight=1.5) + 
    NoLegend() + NoAxes()+ggtitle(i)+
    theme(plot.title = element_text(size = 5,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
Idents(SLE.yaa.combined) <- "cdr3"
for (i in unique(x = names(head(sort(table(SLE.yaa.combined@meta.data$cdr3), decreasing = T), n=4)))) {
  plot.list[[i]] <- DimPlot(
    object = SLE.yaa.combined, cols.highlight="forestgreen",
    cells.highlight = Cells(subset(SLE.yaa.combined, idents=i)),
    sizes.highlight=1.5) + NoLegend() + NoAxes()+ggtitle(i)+
    theme(plot.title = element_text(size = 5,hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
CombinePlots(plots = plot.list, ncol = 4)
ggsave2("topclones.umap.png",width=10, height=5,device="png")

#top clone distribution to different clusters, by mouse
mice<-unique(T_all.combined@meta.data$mouse_ID)
freqlist=list()
for (i in mice) {
  df<-T_all.combined@meta.data[T_all.combined@meta.data$mouse_ID==i,]
  clones<-names(head(sort(table(df$cdr3), decreasing = T), n=5))
  df2<-df[df$cdr3 %in% clones, ]
  freq<-as.data.frame(table(df$my.clusters2))
  freq$pct<-100*freq$Freq/sum(freq$Freq)
  freq$mouse<-i
  freq$condition<-unique(df$condition)
  freqlist[[i]]<-freq
}
topclone.freq.mice<-do.call(rbind,freqlist)
ggbarplot(topclone.freq.mice, x = "Var1", y = "pct",add = c("mean_se", "jitter"),
          color = "condition",position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2.5,label.y=60)+
  # stat_compare_means(method="anova",label.y=0)+
  labs(y = "Cluster distribution of cells belonging \n to 5 largest clonotypes per mouse")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))#,axis.title.y=element_text(size=10))
ggsave2("topclone.clust.distr.mice.png",width=5, height=4,device="png")

#######CloneDE
#public repertoire df
immdata = repLoad("./rep")
df<-clone.data.ab
imm.list<-list()
for(i in unique(df$mouse_ID)){ #creating immunarch list from collapsed data
  df2<-df[df$mouse_ID==i,]
  df2<-transform(df2,Clones=ave(seq(nrow(df2)),cdr3,FUN=length))
  df2<-as.data.table(df2)[, lapply(.SD, data_concater), by=cdr3]
  df2$Clones<-as.numeric(as.character(df2$Clones))
  df2$Proportion<-df2$Clones/sum(df2$Clones)
  imm.list[[paste0(i,"_")]]<-tibble(
    Clones=df2$Clones, Proportion=df2$Proportion, CDR3.nt=df2$cdr3_nt,CDR3.aa=df2$cdr3,
    V.name=df2$v_gene,D.name=df2$d_gene,J.name=df2$j_gene
  )
}
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
pr.B6 = pubRepFilter(pr, immdata$meta, c(condition = "B6"))
pr.SLE.yaa = pubRepFilter(pr, immdata$meta, c(condition = "SLE.yaa"))
pr.B6[is.na(pr.B6)]<-0
pr.SLE.yaa[is.na(pr.SLE.yaa)]<-0
pr.B6[["avgfreq.B6"]] = rowMeans(public_matrix(pr.B6), na.rm = T)
pr.SLE.yaa[["avgfreq.SLE.yaa"]] = rowMeans(public_matrix(pr.SLE.yaa), na.rm = T)
pr.B6[["sum.B6"]] = rowSums(public_matrix(pr.B6)[,1:2], na.rm = T)
pr.SLE.yaa[["sum.SLE.yaa"]] = rowSums(public_matrix(pr.SLE.yaa)[,1:2], na.rm = T)
pr.res.cdr3 = dplyr::full_join(pr.B6, pr.SLE.yaa, by = "CDR3.aa")
pr.res.cdr3.graph<-pr.res.cdr3
pr.res.cdr3.graph<-as.data.frame(pr.res.cdr3.graph)
pr.res.cdr3.graph[is.na(pr.res.cdr3.graph)]<-0
pr.res.cdr3.graph[["Samples.sum"]] = pr.res.cdr3.graph[["Samples.x"]] + pr.res.cdr3.graph[["Samples.y"]]
pr.res.cdr3.graph[["freq.ratio"]] = apply(pr.res.cdr3.graph[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.cdr3.graph[["log2FC"]]= apply(pr.res.cdr3.graph[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log2((x[2])/(x[1])))
rownames(pr.res.cdr3.graph)<-pr.res.cdr3.graph$CDR3.aa
###expanded clone 564 vs B6 scatter
labels<-unique(pr.res.cdr3.graph[abs(pr.res.cdr3.graph$log2FC)>1&pr.res.cdr3.graph$Samples.sum>1&
                           (pr.res.cdr3.graph$sum.B6>1|pr.res.cdr3.graph$sum.SLE.yaa>1),])
ggplot(pr.res.cdr3.graph,aes(x = sum.SLE.yaa, y =  sum.B6,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(shape=21) + #aes(colour=factor(annotate),fill = factor(annotate)), 
  guides(fill = guide_legend(override.aes = list(size = 7)),color=F)+ #,size=F
  labs(x = "SLE.yaa", y = "B6", size="Samples")+#,fill="Expansion"
  theme(legend.direction = "vertical", legend.box = "horizontal")#+#theme(legend.title = element_blank())+
  ggsave2("overlap.scatter.TCRab.expansion.png",width=5.5, height=3,device="png")

# identifying expanded clones
SLE.yaa.expanded <- names(table(SLE.yaa.combined@meta.data$cdr3))[table(SLE.yaa.combined@meta.data$cdr3) > 1]
SLE.yaa.unexpanded <- names(table(SLE.yaa.combined@meta.data$cdr3))[table(SLE.yaa.combined@meta.data$cdr3) ==1]
B6.expanded <- names(table(B6.combined@meta.data$cdr3))[table(B6.combined@meta.data$cdr3) > 1]
B6.unexpanded <- names(table(B6.combined@meta.data$cdr3))[table(B6.combined@meta.data$cdr3) ==1]

#expanded DE (by condition)
Idents(SLE.yaa.combined) <- "cdr3"
expandedDE.SLE.yaa<-FindMarkers(SLE.yaa.combined, ident.1 = SLE.yaa.expanded, ident.2 = SLE.yaa.unexpanded, 
                             min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.SLE.yaa$gene<-rownames(expandedDE.SLE.yaa)
Idents(B6.combined) <- "cdr3"
expandedDE.B6<-FindMarkers(B6.combined, ident.1 = B6.expanded, ident.2 = B6.unexpanded, 
                            min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
expandedDE.B6$gene<-rownames(expandedDE.B6)
exp.vs.expandedDE <-merge(expandedDE.SLE.yaa[,c("gene","avg_log2FC")], expandedDE.B6[,c("gene","avg_log2FC")], by="gene")
rownames(exp.vs.expandedDE)<-exp.vs.expandedDE$gene

Idents(T_all.combined) <- "cdr3"
SLEexp.vs.B6exp.DE<-FindMarkers(T_all.combined, ident.1 = SLE.yaa.expanded, ident.2 = B6.expanded, 
                                min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
SLEexp.vs.B6exp.DE$gene<-rownames(SLEexp.vs.B6exp.DE)

save(T_all.combined,file="T_all.RData")
save(list=c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"),"clusters"),file="DE.dfs.RData")
save.image("vdj.analyzed.RData")

#Volcano Plots
for(i in c(ls(pattern=".cdr3.autoimmune.response"),ls(pattern="exp.vs.expandedDE"),ls(pattern="expandedDE.ms"),ls(pattern="exp.DE"))){
  df<-get(i)
  if(!is(df,"try-error")){
    if(!("p_val_adj" %in% colnames(df))){
      if("ttest" %in% colnames(df) ){df$p_val_adj<-p.adjust(df$ttest,method="BH")}
      if("p_val" %in% colnames(df) ){df$p_val_adj<-p.adjust(df$p_val,method="BH")}
    }
    if(!("avg_log2FC" %in% colnames(df))){
      if("delta" %in% colnames(df)){df$avg_log2FC<-df$delta}
      if("avg_logFC" %in% colnames(df)){df$avg_log2FC<-log2(exp(df$avg_logFC))}
    }
    if(is.numeric(df$avg_log2FC)&is.numeric(df$p_val_adj)){
      df2<-df[df$p_val_adj<0.01,]
      gene.list<-head(rownames(df2[order(df2$p_val_adj),]),n=15)
      EnhancedVolcano(df,lab = rownames(df),
                      x = 'avg_log2FC',y = 'p_val_adj',title = i,col=c("black","black","black","red3"),
                      selectLab=gene.list,#xlab=bquote(~Log[2]~ (frac("SLE.yaa","B6"))),
                      pCutoff = 0.01,FCcutoff = 0.2,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                      legendPosition="none",#drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                      subtitle="", caption="",border="full",cutoffLineWidth=0,
                      gridlines.major=F,gridlines.minor=F,titleLabSize=10
      )
      ggsave2(paste0(i,".volcano.png"),width=4, height=4,device="png")
    }
  }
}


