rm(list=ls())
library(biomaRt)
library(EnhancedVolcano)
library(cowplot)
library(dplyr)
library(ggpubr)
library(Seurat)
library(data.table)
library(VennDiagram)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

load("scTfh/analyzed.RData") #load DEs from scTfh data (see https://github.com/egarren/scTfh)
scTfh.raw<-T_all.combined.autoimmune.response
scTfh.raw$avg_log2FC<-scTfh.raw$log2FC
scTfh.raw$gene<-rownames(scTfh.raw)
rm(list=setdiff(ls(),"scTfh.raw"))
  
#mouse vs human
load("T_cell/analyzed.RData") #load DEs from CD4 subset
cluster0.autoimmune.response$gene<-rownames(cluster0.autoimmune.response)
SLE.yaa<-data.frame(gene=rownames(cluster0.autoimmune.response),
                  SLE.yaa.log2FC=cluster0.autoimmune.response$avg_log2FC)
chimera<-data.frame(gene=rownames(scTfh.raw),
                    chimera.log2FC=scTfh.raw$avg_log2FC)


#venn diagram
DE.list<-list()
for(i in c("scTfh.raw","cluster0.autoimmune.response")){
  df<-get(i)
  DE<-df[abs(df$avg_log2FC)>0.2&df$p_val_adj<0.05,]
  DE.list[[i]]<-DE$gene
}
names(DE.list)<-c("564Igi chimera","SLE.yaa")
venn.diagram(
  x = DE.list,category.names = names(DE.list),filename = "SLE.chim.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("#440154ff", '#21908dff'),fill=c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("#440154ff", '#21908dff'),cat.dist = c(0.1, 0.1)
)
intersect(DE.list[[1]],DE.list[[2]]) 



#correlation
options(ggrepel.max.overlaps = Inf)
comp<-merge(SLE.yaa,chimera,by="gene")
comp<-comp[!is.na(comp$chimera.log2FC),]
rownames(comp)<-make.unique(as.character(comp$gene))
comp$sig<-"unsig"
comp$sig[abs(comp$SLE.yaa.log2FC)>0.5&abs(comp$chimera.log2FC)>0.2]<-"co.DE"
comp$sig[abs(comp$SLE.yaa.log2FC)>0.5&abs(comp$chimera.log2FC)<=0.2]<-"SLE.yaa.DE"
comp$sig[abs(comp$SLE.yaa.log2FC)<=0.5&abs(comp$chimera.log2FC)>0.2]<-"chimera.DE"
SLE.yaa.up<-rownames(comp[abs(comp$SLE.yaa.log2FC)>0.5&!is.na(comp$gene),])
if(length(SLE.yaa.up)>30){SLE.yaa.up<-rownames(comp[abs(comp$SLE.yaa.log2FC)>0.7&!is.na(comp$gene),])}
if(length(SLE.yaa.up)>30){SLE.yaa.up<-rownames(comp[abs(comp$SLE.yaa.log2FC)>1&!is.na(comp$gene),])}
if(length(SLE.yaa.up)>30){SLE.yaa.up<-rownames(comp[abs(comp$SLE.yaa.log2FC)>1.5&!is.na(comp$gene),])}
chimera.up<-rownames(comp[abs(comp$chimera.log2FC)>0.2&!is.na(comp$gene),])
if(length(chimera.up)>15){chimera.up<-rownames(comp[abs(comp$chimera.log2FC)>0.3&!is.na(comp$gene),])}
if(length(chimera.up)>15){chimera.up<-rownames(comp[abs(comp$chimera.log2FC)>0.4&!is.na(comp$gene),])}
if(length(chimera.up)>15){chimera.up<-rownames(comp[abs(comp$chimera.log2FC)>0.5&!is.na(comp$gene),])}
if(length(chimera.up)>15){chimera.up<-rownames(comp[abs(comp$chimera.log2FC)>0.75&!is.na(comp$gene),])}
diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>0.5&!is.na(comp$gene),])
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>1&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>2&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>3&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>4&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>5&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>6&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>10&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>15&!is.na(comp$gene),])}
if(length(diff)>15){diff<-rownames(comp[abs(comp$SLE.yaa.log2FC-comp$chimera.log2FC)>20&!is.na(comp$gene),])}
gene.list<-unique(c(SLE.yaa.up,chimera.up,diff))
comp$sig <- factor(comp$sig, levels = c("unsig","co.DE","SLE.yaa.DE","chimera.DE"))
p2 <- ggplot(comp, aes(SLE.yaa.log2FC,chimera.log2FC)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) + 
  labs(x=bquote(~Log[2]~ (frac("SLE.yaa","B6"))),y=bquote(~Log[2]~ (frac("564Igi","AID"))))+theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
        axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
  geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
  scale_color_manual(values = c("grey","black","red","orange"))
if(length(gene.list)>0 & length(gene.list)<=60){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
p2
ggsave2("chimera.SLE.yaa.scatter.png",width=5, height=5,device="png")

