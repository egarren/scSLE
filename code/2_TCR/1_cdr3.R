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
library(tools)
library(dplyr)
# library(alakazam)
library(ggalluvial)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(circlize)
library(useful)
library(ggseqlogo)
library(yingtools2)
library(qgraph)
library(vegan)
gm_mean = function(x, na.rm=TRUE){  #geometric mean functions for graphing
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 
gm_ymax<- function(x, na.rm = TRUE)
{
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))+exp(sd(log(x), na.rm = na.rm))
}
barcoder<-function(df, prefix){
  df$ms.barcoder<-prefix
  df$orig.barcode<-df$barcode
  # df$barcode <- gsub("\\-1", "", df$barcode)
  df$barcode <- paste0(prefix, df$barcode)
  df$ms_clone <- paste0(prefix, df$raw_clonotype_id) 
  df$ms_v <- paste0(prefix, df$v_gene) 
  df$ms_j<- paste0(prefix, df$j_gene) 
  df$ms_cdr3 <- paste0(prefix, df$cdr3) 
  df
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
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n),na.rm=T)
  breaks[!duplicated(breaks)]
}



#load cdr3 data
metadata<-read.csv("SLE.metadata.csv", header=T) #define path
metadata[] <- lapply(metadata, as.character)
mice<-metadata$mouse_ID
for (i in mice){
  csv.path<-paste0(i,"/outs/vdj_t/filtered_contig_annotations.csv")
  fasta.path<-paste0(i,"/outs/vdj_t/filtered_contig.fasta")
  df<-merge(read.csv(csv.path, header=T),read.fasta(fasta.path),by.x="contig_id",by.y="seq.name")
  colnames(df)[colnames(df)=="length"]<-"vdj_length"
  assign(paste0(i,".t.vdj"),barcoder(df,prefix=paste0(i,"_")))
  csv.path<-paste0(i,"/outs/vdj_b/filtered_contig_annotations.csv")
  fasta.path<-paste0(i,"/outs/vdj_b/filtered_contig.fasta")
  df<-merge(read.csv(csv.path, header=T),read.fasta(fasta.path),by.x="contig_id",by.y="seq.name")
  colnames(df)[colnames(df)=="length"]<-"vdj_length"
  assign(paste0(i,".b.vdj"),barcoder(df,prefix=paste0(i,"_")))
}

#rename barcodes
Tall.vdj<-rbind(m232.t.vdj,m233.t.vdj,m234.t.vdj,m235.t.vdj)
Ball.vdj<-rbind(m232.b.vdj,m233.b.vdj,m234.b.vdj,m235.b.vdj)
combined.vdj<-rbind(Tall.vdj,Ball.vdj)

#paired data
productive.barcode<-combined.vdj[combined.vdj$productive=="true"&!is.na(combined.vdj$productive),] 
productive.barcode<-as.data.table(productive.barcode)[, lapply(.SD, data_concater3), by=barcode]
productive.barcode<-productive.barcode[productive.barcode$chain %in% c("IGH+IGK","IGH+IGL","TRA+TRB"),]
productive.barcode<-transform(productive.barcode, cdr3.freq = ave(seq(nrow(productive.barcode)), cdr3, FUN=length))
productive.barcode<-transform(productive.barcode, ms_cdr3.freq = ave(seq(nrow(productive.barcode)), ms_cdr3, FUN=length))
productive.barcode<-transform(productive.barcode, vab.freq = ave(seq(nrow(productive.barcode)), v_gene, FUN=length))
rownames(productive.barcode)<-productive.barcode$barcode

load("graphed.RData")

#extract and add Seurat metadata
clone.data<-Tall.vdj[Tall.vdj$productive=="true"&!is.na(Tall.vdj$productive)&Tall.vdj$chain %in% c("TRA","TRB"),] 
SLE.obj.combined@meta.data$orig.barcode<-rownames(SLE.obj.combined@meta.data)
clone.data<-left_join(x = clone.data, y = subset(SLE.obj.combined@meta.data, select=-c(ms.barcoder)),by = c("barcode"="orig.barcode"),keep=F)
clone.data<-left_join_replace(clone.data,metadata,by="ms.barcoder")
clone.data.seurat<-clone.data[!is.na(clone.data$my.clusters),] 
clone.meta<- as.data.table(clone.data)[, lapply(.SD, data_concater), by=cdr3]
clone.meta$aa_length<-nchar(clone.meta$cdr3)
clusters<-unique(clone.data.seurat$my.clusters)
mice<-unique(clone.data$mouse_ID)
B6.mice<-unique(clone.data[clone.data$condition=="B6",]$mouse_ID)
SLE.yaa.mice<-unique(clone.data[clone.data$condition=="SLE.yaa",]$mouse_ID)
save(list=c(ls(pattern="clone."),"metadata","Tall.vdj","clusters",ls(pattern="mice")),file="clone.data.RData")

#TCR venn diagram
df<-clone.data
dz.list<-list()
for(i in unique(df$chain)){
  dz.list[[i]]<-unique(df$barcode[df$chain==i])
}
venn.diagram(
  x = dz.list,category.names = names(dz.list),filename = "TR.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("#440154ff", '#21908dff'),fill=c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("#440154ff", '#21908dff'),cat.dist = c(0.1, 0.1)
)

##CDR3 property analysis by cluster
#WebLogo
#all clones w/n all mice, comparing conditions
plot.list<-list()
for(j in c("B6","SLE.yaa")){
  for(i in c("TRA","TRB")){
    df<-clone.data[clone.data$chain==i & clone.data$condition ==j,]
    freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
    freq.tab$aa_length<-nchar(as.character(freq.tab$Var1))
    for(k in 14){ #10:18
      if(j=="B6"){name="WT"}else{name="SLE.yaa"}
      df2<-freq.tab[freq.tab$aa_length==k,]
      plot.list[[paste0(i,j,k)]]<-ggplot()+geom_logo(as.character(df2$Var1))+
        theme_logo()+ggtitle(paste0(i,"_",name))+
        theme(axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5)) #,"_",k
    }
  }
}
CombinePlots(plot.list,ncol=2,legend="bottom")
ggsave2("cdr3.weblogo.bycondition.png",width=6,height=4,device="png")

#aa length
#all clones w/n all mice, comparing conditions
plot.list<-list()
plot.list2<-list()
for(i in c("TRA","TRB")){
  df.list<-list()
  for(j in c("B6","SLE.yaa")){
    df<-clone.data[clone.data$chain==i & clone.data$condition ==j,]
    freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
    freq.tab$aa_length<-nchar(as.character(freq.tab$Var1))
    if(j=="B6"){freq.tab$condition<-"WT"}else{freq.tab$condition<-"SLE.yaa"}
    df.list[[j]]<-freq.tab #add to list for each condition
  }
  df<-rbindlist(df.list)
  df$condition <- factor(df$condition, levels = c("WT","SLE.yaa"))
  mu <- ddply(df, "condition", summarise, grp.mean=mean(aa_length))
  plot.list[[i]]<-ggplot(df, aes(x=aa_length,color=condition,fill=condition)) + 
    geom_histogram(aes(y=..density..),binwidth=1,position="identity",alpha=0.2) +#geom_density(fill=NA)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+ggtitle(i)+theme_classic()+
    labs(x = "CDR3 Length", y = "") +theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())+
    scale_color_manual(values = c("black", "red"))+
    scale_fill_manual(values = c("black", "red"))
  plot.list2[[i]]<-ggbarplot(df, x = "condition", y = "aa_length", add = c("mean","jitter"),
                             position = position_dodge(0.8), legend="right")+
    stat_compare_means(aes(group = condition))+
    labs(x="condition", y = "CDR3 Length")+ggtitle(i)
}
CombinePlots(plot.list,ncol=2,legend="right")
ggsave2("cdr3.aa_length.histo.png",width=7,height=3,device="png")
CombinePlots(plot.list2,ncol=2,legend="right")
ggsave2("cdr3.aa_length.bar.dot.png",width=8,height=4,device="png")

#Heatmap by mouse v usage
plot.list<-list()
for(i in c("TRA","TRB")){
  df<-clone.data[clone.data$chain==i,]
  tab<-table(df$v_gene,df$mouse_ID)
  Proportions <- scale(tab,scale=colSums(tab),center=FALSE)*100 
  annot.col<-data.frame(BMChimera=metadata$condition)
  annot.col$BMChimera[annot.col$BMChimera %in% "B6"]<-"WT"
  annot.col$BMChimera[annot.col$BMChimera %in% "SLE.yaa"]<-"SLE.yaa"
  rownames(annot.col)<-metadata$mouse_ID
  plot.list[[i]]<-pheatmap(as.data.frame.matrix(Proportions),annotation_col=annot.col,main=i,annotation_names_col=F,
                           annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red")),fontsize_row=5)
}
save_pheatmap_png <- function(x, filename, width=700, height=1200, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(plot.list[["TRA"]], "TRA.v_gene.heatmap.png")
save_pheatmap_png <- function(x, filename, width=700, height=600, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(plot.list[["TRB"]], "TRB.v_gene.heatmap.png")

###VDJtools export
dir.create("./vdjtools")
setwd("./vdjtools")
mice<-metadata$ms.barcoder
#by condition 
file.list=list()
for (i in mice) { 
  db<-clone.data[clone.data$ms.barcoder==i,]
  db<-transform(db, ab.cdr3.freq = ave(seq(nrow(db)), cdr3, FUN=length)) #calculate cdr3 frequencies
  db<- as.data.table(db)[, lapply(.SD, data_concater2), by=cdr3]   #collapse by cdr3
  db$ab.cdr3.freq <- as.numeric(as.character(db$ab.cdr3.freq))
  vdj <- db %>% transmute(count=ab.cdr3.freq,   #make vdjtools table
                          freq=ab.cdr3.freq/sum(ab.cdr3.freq, na.rm=TRUE),
                          cdr3nt=cdr3_nt,
                          cdr3aa=cdr3,
                          v=v_gene,
                          d=d_gene,
                          j=j_gene)
  write.table(vdj, file=paste0(i,"_vdjtools.tab"), quote=FALSE, sep="\t", row.names=FALSE)
  file.list[i]<-paste0(i,"_vdjtools.tab")
}
vdjtools.metadata<-data.frame(file.name=unlist(file.list,use.names=F),sample.id=names(file.list))
vdjtools.metadata<-left_join(x = vdjtools.metadata, y = metadata, by = c("sample.id"="ms.barcoder"),keep=F)
write.table(vdjtools.metadata, file = "vdjtools.metadata.txt", sep = "\t",quote=F,row.names = FALSE)
setwd("../")
#by cluster
dir.create("./vdjtools.cluster")
setwd("./vdjtools.cluster")
vdjtools.metadata2<-data.frame()
for (i in mice) { 
  db<-clone.data[clone.data$ms.barcoder==i,]
  for(j in clusters){
    db<-db[!is.na(db$my.clusters),]
    db2<-db[db$my.clusters==j,]
    db2<-transform(db2, ab.cdr3.freq = ave(seq(nrow(db2)), cdr3, FUN=length)) #calculate cdr3 frequencies
    db2<- as.data.table(db2)[, lapply(.SD, data_concater2), by=cdr3]   #collapse by cdr3
    db2$ab.cdr3.freq <- as.numeric(as.character(db2$ab.cdr3.freq))
    vdj <- db2 %>% transmute(count=ab.cdr3.freq,   #make vdjtools table
                             freq=ab.cdr3.freq/sum(ab.cdr3.freq, na.rm=TRUE),
                             cdr3nt=cdr3_nt,
                             cdr3aa=cdr3,
                             v=v_gene,
                             d=d_gene,
                             j=j_gene)
    write.table(vdj, file=paste0(i,".clust",j,"_vdjtools.tab"), quote=FALSE, sep="\t", row.names=FALSE)
    vdjtools.metadata2[paste0(i,".clust",j),1]<-paste0(i,".clust",j,"_vdjtools.tab")
    vdjtools.metadata2[paste0(i,".clust",j),2]<-paste0(i,".clust",j)
    vdjtools.metadata2[paste0(i,".clust",j),3]<-i
    vdjtools.metadata2[paste0(i,".clust",j),4]<-j
  }
}
colnames(vdjtools.metadata2)<-c("file.name","sample.id","mouse","Cluster")
vdjtools.metadata2<-left_join(x = vdjtools.metadata2, y = metadata, by = c("mouse"="ms.barcoder"),keep=F)
write.table(vdjtools.metadata2, file = "vdjtools.clust.metadata.txt", sep = "\t",quote=F,row.names = FALSE)
setwd("../")

##DeepTCR export
#scTfh data
dir.create("./deepTCR")
setwd("./deepTCR")
dir.create("./SLE.data")
setwd("./SLE.data")
for(i in names(table(clone.data$condition))){
  df<-clone.data[clone.data$condition==i,]
  dir.create(paste0("./",i))
  setwd(paste0("./",i))
  for(j in names(table(df$mouse_ID))){
    df2<-df[df$mouse_ID==j,]
    df.TRA<-df2[df2$chain=="TRA",]
    df.TRB<-df2[df2$chain=="TRB",]
    cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene","d_gene")],
                      df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
    cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters")]),by="barcode")
    cells$clone_ID<-paste0(cells$cdr3.x,cells$v_gene.x,cells$j_gene.x,
                           cells$cdr3.y,cells$v_gene.y,cells$j_gene.y,cells$mouse_ID)
    cells.cloneID<-transform(cells,freq=ave(seq(nrow(cells)),clone_ID,FUN=length))
    cells.cloneID<-as.data.table(cells.cloneID)[, lapply(.SD, data_concater), by=clone_ID]
    cells.cloneID$freq<-as.numeric(as.character(cells.cloneID$freq))
    write.table(cells.cloneID,file=paste0(j,".deepTCR.tsv"),sep="\t",row.names=F,col.names=F,quote=F)
  }
  setwd("../")
}
setwd("../")


#paired data
clone.data.ab<-Tall.vdj[Tall.vdj$productive=="true"&!is.na(Tall.vdj$productive)&Tall.vdj$chain %in% c("TRA","TRB"),] 
clone.data.ab<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater3), by=barcode]
clone.data.ab<-clone.data.ab[clone.data.ab$chain=="TRA+TRB",]
clone.data.ab<-transform(clone.data.ab, cdr3.freq = ave(seq(nrow(clone.data.ab)), cdr3, FUN=length))
clone.data.ab<-transform(clone.data.ab, ms_cdr3.freq = ave(seq(nrow(clone.data.ab)), ms_cdr3, FUN=length))
clone.data.ab<-transform(clone.data.ab, vab.freq = ave(seq(nrow(clone.data.ab)), v_gene, FUN=length))
SLE.obj.combined@meta.data$orig.barcode<-rownames(SLE.obj.combined@meta.data)
clone.data.ab<-left_join(x = clone.data.ab, y = subset(SLE.obj.combined@meta.data, select=-c(ms.barcoder)), by = c("barcode"="orig.barcode"),keep=F)
clone.data.ab<-left_join_replace(clone.data.ab,metadata,by="ms.barcoder")
clone.data.ab.seurat<-clone.data.ab[!is.na(clone.data.ab$my.clusters),] 
clone.meta.ab<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater), by=cdr3]
clone.meta.ab$cdr3.freq<-as.numeric(as.character(clone.meta.ab$cdr3.freq))
# clone.meta.ab$aa_length<-nchar(clone.meta.ab$cdr3)
clone.data.ab.ms_cdr3<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater), by=ms_cdr3]
clone.data.ab.ms_cdr3$cdr3.freq <- as.numeric(as.character(clone.data.ab.ms_cdr3$cdr3.freq))
clone.data.ab.ms_cdr3$ms_cdr3.freq <- as.numeric(as.character(clone.data.ab.ms_cdr3$ms_cdr3.freq))
clone.data.ab$clust_cdr3<-paste0(clone.data.ab$my.clusters2,"_",clone.data.ab$cdr3)
clone.data.ab.clust_cdr3<-as.data.table(clone.data.ab)[, lapply(.SD, data_concater), by=clust_cdr3]
clone.data.ab.clust_cdr3$cdr3.freq <- as.numeric(as.character(clone.data.ab.clust_cdr3$cdr3.freq))
clusters<-unique(clone.data.seurat$my.clusters)
mice<-unique(clone.data$mouse_ID)
B6.mice<-unique(clone.data[clone.data$condition=="B6",]$mouse_ID)
SLE.yaa.mice<-unique(clone.data[clone.data$condition=="SLE.yaa",]$mouse_ID)
save(list=c(ls(pattern="clone."),"metadata","Tall.vdj","clusters",ls(pattern="mice")),file="clone.ab.data.RData")

##Venn diagram
#between mice
df<-clone.data.ab
dz.list<-list()
for(i in unique(df$mouse_ID)){
  dz.list[[i]]<-unique(df$cdr3[df$mouse_ID==i])
}
venn.diagram(
  x = dz.list,filename = "mouse.cdr3.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  # lwd=1,col=c("black", 'red'),fill=c(alpha("black",0.3), alpha('red',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  # cat.col = c("black", 'red'),cat.dist = c(0.1, 0.1),
  ext.pos=180,ext.line.lwd=0.25,ext.dist=-0.2
)
intersect(dz.list[["m234"]],dz.list[["m235"]])
#between conditions
df<-clone.data.ab
dz.list<-list()
for(i in unique(df$condition)){
  dz.list[[i]]<-unique(df$cdr3[df$condition==i])
}
venn.diagram(
  x = dz.list,category.names = c("WT","SLE.yaa"),filename = "condition.cdr3.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("black", 'red'),fill=c(alpha("black",0.3), alpha('red',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("black", 'red'),cat.dist = c(0.1, 0.1),
  ext.pos=180,ext.line.lwd=0.25,ext.dist=-0.2
)
#between clusters
df<-clone.data.ab[clone.data.ab$my.clusters %in% c("0","1","2","3"),]
dz.list<-list()
for(i in c("CD4","CD8","CD8-eff","Treg")){
  dz.list[[i]]<-unique(df$cdr3[df$my.clusters2==i])
}
venn.colors<-hue_pal()(7)[1:4]
venn.diagram(
  x = dz.list,category.names = names(dz.list),filename = "clusters.cdr3.venn.png",
  imagetype="png",height=1200,width=1200,margin=0.1,
  lwd=1,col=venn.colors,fill=alpha(venn.colors,0.3),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.5,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = venn.colors#,cat.dist = c(0.1, 0.1,0.1,0.1)
)

##clone size
# all clones w/n all mice, comparing conditions
df.list<-list()
for(j in unique(clone.data.ab$condition)){
  df<-clone.data.ab[clone.data.ab$condition ==j,]
  freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
  freq.tab$condition<-j
  # if(j=="B6"){freq.tab$condition<-"WT"}else{freq.tab$condition<-"SLE.yaa"}
  df.list[[j]]<-freq.tab #add to list for each cluster
}
df<-rbindlist(df.list)
ggbarplot(df, x = "condition", y = "Freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="none",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),method="wilcox.test")+
  labs(x="", y = "Clone Size") #+scale_y_log10()
ggsave2("bar.dot.condition.cdr3.freq.png",width=3, height=4,device="png")
mu <- ddply(df, "condition", summarise, grp.mean=gm_mean(Freq))
df$condition <- factor(df$condition, levels = c("WT","SLE.yaa"))
ggplot(df, aes(x=Freq,color=condition,fill=condition)) + 
  geom_histogram(aes(y=..density..),bins=10,position="identity",alpha=0.2) +#geom_density(fill=NA)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+theme_classic()+
  labs(x = "Clone Size", y = "Density")+scale_x_log10()+theme(legend.title = element_blank())+
  scale_color_manual(values = c("black", "red"))+
  scale_fill_manual(values = c("black", "red"))
ggsave2("histogram.condition.cdr3.freq.png",width=4, height=3,device="png")
#all clones w/n all mice w/n cluster
df.list <- list()
plot.list<-list()
mycluster.names<-names(sort(table(clone.data.ab.seurat$my.clusters2),decreasing=T))
for (i in mycluster.names){
  for(j in c("B6","SLE.yaa")){
    df<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters2==i&clone.data.ab.seurat$condition==j,] #gating on given cluster
    freq.tab<-as.data.frame(table(df$cdr3))  #ms_cdr3 frequency table for given cluster
    freq.tab$condition<-j
    freq.tab$my.clusters2<-i
    df.list[[paste0(i,j)]]<-freq.tab #add to list for each cluster
  }
  df2<-rbindlist(df.list)
  
}
df<-rbindlist(df.list)
ggbarplot(df, x = "my.clusters2", y = "Freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  # stat_compare_means(aes(group = condition),label = "p.format", method="wilcox.test",label.y=200,size=4)+
  labs(y = "Clone Size")+theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(limits = mycluster.names) #+scale_y_log10()
ggsave2("bar.dot.cluster.cdr3.freq.png",width=4, height=3,device="png")

##top 10 cdr3 percent (diversity), within cluster, by mouse
df.list <- list()
for (i in unique(clone.data.ab.seurat$my.clusters2)){
  df.cluster<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters2==i,] #gating on given cluster
  ms.clone.freq<-as.data.frame(table(df.cluster$ms_cdr3)) #ms_cdr3 frequency table for given cluster
  ms.clone.freq$ms_cdr3<-as.character(ms.clone.freq$Var1) 
  for(k in c("SLE.yaa","B6")){
    meta2<-metadata[metadata$condition==k,]
    meta<-data.frame(mice=names(table(meta2$ms.barcoder))) #blank dataframe of mice (to fill in upcoming loop)
    meta$condition<-k
    # if(k=="B6"){meta$condition<-"WT"}else{meta$condition<-"SLE.yaa"}
    meta$top10.pct<-0
    meta$my.clusters2<-i
    for (j in 1:nrow(meta)) {
      df<-ms.clone.freq[grepl(meta$mice[j],ms.clone.freq$ms_cdr3), ] #frequency info for cdrs3 from given mouse
      df <- df[order(df$Freq,decreasing=T),] #calculating what % of all clones the top 10 clones represent
      pct <- 100*(df$Freq/sum(df$Freq))
      meta$top10.pct[j]<-sum(pct[1:10]) #add to list for each mouse
    }
    df.list[[paste0(i,k)]]<-meta
  }
}
df<-rbindlist(df.list)
ggbarplot(df, x = "my.clusters2", y = "top10.pct", color="condition", add = c("mean_se","jitter"),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2.5,label.y=110)+
  labs(y = "% of repertoire in \n 10 larges clones")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(limits = mycluster.names) #+scale_y_log10()
ggsave2("bar.dot.cluster.cdr3.top10freq.png",width=5, height=3,device="png")

#diversity score (all cells)
df.list <- list()
for (i in unique(clone.data.ab$condition)){
  df.cluster<-clone.data.ab[clone.data.ab$condition==i,] #gating on given cluster
  meta2<-metadata[metadata$condition==i,]
  meta<-data.frame(mice=names(table(meta2$ms.barcoder))) #blank dataframe of mice (to fill in upcoming loop)
  meta$condition<-i
  for(l in c("shannon","simpson","invsimpson")){
    meta[[l]]<-0
    for (j in 1:nrow(meta)) {
      ms.clust<-df.cluster[df.cluster$ms.barcoder==meta$mice[j],]
      meta[[l]][j]<-diversity(table(ms.clust$cdr3),index=l)
    }
  }
  df.list[[i]]<-meta
}
df<-rbindlist(df.list)
ggbarplot(df, x = "condition", y = "shannon", color="condition", add = c("mean_se","jitter"),add.params = list(size = 0.5),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),method="t.test")+
  labs(x="", y = "Diversity") #+scale_y_log10()
ggsave2(paste0("bar.dot.all.cdr3.div.png"),width=4, height=3,device="png")
##Diversity scores (within cluster)
df.list <- list()
for (i in unique(clone.data.ab.seurat$my.clusters2)){
  df.cluster<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters2==i,] #gating on given cluster
  for(k in c("SLE.yaa","B6")){
    meta2<-metadata[metadata$condition==k,]
    meta<-data.frame(mice=names(table(meta2$ms.barcoder))) #blank dataframe of mice (to fill in upcoming loop)
    meta$condition<-k
    meta$my.clusters2<-i
    for(l in c("shannon","simpson","invsimpson")){
      meta[[l]]<-0
      for (j in 1:nrow(meta)) {
        ms.clust<-df.cluster[df.cluster$ms.barcoder==meta$mice[j],]
        meta[[l]][j]<-diversity(table(ms.clust$cdr3),index=l)
      }
    }
    df.list[[paste0(i,k)]]<-meta
  }
}
df<-rbindlist(df.list)
ggbarplot(df, x = "my.clusters2", y = "shannon", color="condition", add = c("mean_se","jitter"),add.params = list(size = 0.5),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2.5,label.y=7)+ #
  labs(y = "Diversity")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(limits = mycluster.names) #+scale_y_log10()
ggsave2(paste0("bar.dot.cluster.cdr3.div.png"),width=4, height=3,device="png")


#Rarefaction plot (https://peat-clark.github.io/BIO381/veganTutorial.html)
for(i in c("0","1","2","3","4","5","all")){
  df.cluster<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters==i,] #gating on given cluster
  if(i=="all"){df.cluster<-clone.data.ab}
  raref.df<-as.data.frame.matrix(table(df.cluster$mouse_ID,df.cluster$cdr3))
  colors<-rep("black",nrow(raref.df))
  colors[rownames(raref.df) %in% SLE.yaa.mice]<-"red"
  # rarecurve(raref.df)
  png(paste0("clust.",i,".raref.png"),width=3.5,height=3.5,units="in",res=200)
  rarecurve(raref.df, step = 20, sample = min(rowSums(raref.df)), col = alpha(colors,0.6), cex = 0.2,
            ylab="Clones",label=F,lwd=4)
  dev.off()
}

#non-metric multidimensional scaling (by condition)
df<-clone.data.ab.seurat
community_matrix<-as.data.frame.matrix(table(df$my.clusters2,df$cdr3))
example_NMDS=metaMDS(community_matrix, k=5,trymax=1000,threshold=0.99)#"The stress, or the disagreement between 2-D configuration and predicted values from the regression"
#A good rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation
png("tcr.NMDS.cluster.stres.png",width=4,height=4,units="in",res=200)
stressplot(example_NMDS,p.col="forestgreen",l.col="black")
dev.off()
# plot(example_NMDS)
colors=hue_pal()(6)
elevation=runif(14,0.5,1.5)
png("tcr.NMDS.cluster.png",width=6,height=6,units="in",res=200)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="black",pch=21,bg=adjustcolor("grey",0.5),pcex=1,air=50000,cex=0.1)
orditorp(example_NMDS,display="sites",col=colors,air=0.0001,cex=1.25)
dev.off()


##Clonotype Distribution Heatmaps
#cluster propotions by CDR3
clone.data.ab.top<-clone.data.ab.seurat[clone.data.ab.seurat$cdr3.freq>2,]
clone.data.ab.top<-clone.data.ab.top[clone.data.ab.top$cdr3!="None",]
Counts <- table(as.character(clone.data.ab.top$my.clusters2),as.character(clone.data.ab.top$cdr3))
Proportions <- scale(Counts,scale=colSums(Counts),center=FALSE)*100 
pheatmap(as.data.frame.matrix(Proportions))#,cutree_cols = 6
#add heatmap metadata
clone.meta.top<- as.data.table(clone.data.ab.top)[, lapply(.SD, data_concater), by=cdr3]
clone.meta.top$cdr3.freq <- as.numeric(as.character(clone.meta.top$cdr3.freq))
annot.col<-data.frame(cdr3=colnames(Proportions))
annot.col<-merge(annot.col,clone.meta.top[,c("cdr3","cdr3.freq","condition")],by="cdr3")
names(annot.col)[names(annot.col) == "condition"] <- "BMChimera"
names(annot.col)[names(annot.col) == "cdr3.freq"] <- "Clone Size"
annot.col$BMChimera[annot.col$BMChimera %in% "B6"]<-"WT"
annot.col$BMChimera[annot.col$BMChimera %in% "SLE.yaa"]<-"SLE.yaa"
annot.col$BMChimera[annot.col$BMChimera %in% "B6+SLE.yaa"]<-"Both"
row.names(annot.col) <- annot.col$cdr3
annot.col$cdr3 <- NULL
p1<-pheatmap(as.data.frame.matrix(Proportions),annotation_col=annot.col,fontsize_col=3,annotation_names_col=F,
             annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red","Both"=alpha("red",0.5))))#cutree_cols = 5,
save_pheatmap_png <- function(x, filename, width=1400, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.cluster.heatmap.png")

##Clonotype Overlap Heatmaps
##between mice
#by number of shared  clones (between mice)
ms.overlap.df<-data.frame()
ms.overlap.list<-list() #list of lists (in case want to extract shared clone info)
for(i in unique(clone.data.ab$mouse_ID)){
  mouse1.clones<-clone.data.ab[clone.data.ab$mouse_ID==i,]
  for (j in unique(clone.data.ab$mouse_ID)){
    if(i!=j){
      mouse2.clones<-clone.data.ab[clone.data.ab$mouse_ID==j,]
      ms.overlap.list[[paste0(i,"+",j)]]<-intersect(mouse1.clones$cdr3,mouse2.clones$cdr3)
      ms.overlap.df[j,i]<-length(ms.overlap.list[[paste0(i,"+",j)]])
    }else{ms.overlap.df[j,i]<-NA}
  }
}
pheatmap(ms.overlap.df)
annot.meta<- as.data.table(clone.data.ab)[, lapply(.SD, data_concater), by=mouse_ID]
annot.col<-data.frame(mouse_ID=colnames(ms.overlap.df))
annot.col<-merge(annot.col,annot.meta[,c("mouse_ID","condition")],by="mouse_ID")
names(annot.col)[names(annot.col) == "condition"] <- "BMChimera"
annot.col$BMChimera[annot.col$BMChimera %in% "B6"]<-"WT"
annot.col$BMChimera[annot.col$BMChimera %in% "SLE.yaa"]<-"SLE.yaa"
row.names(annot.col) <- annot.col$mouse_ID
annot.col$mouse_ID <- NULL
#better coloring
mat_breaks <- quantile_breaks(ms.overlap.df, n = 101)
colfunc <- colorRampPalette(c("#EDF8E9", "#006D2C"))
p1<-pheatmap(ms.overlap.df,annotation_col=annot.col,annotation_names_col=F,show_colnames = F,
             annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red")),
             color=colfunc(length(mat_breaks) - 1),breaks=mat_breaks)
save_pheatmap_png <- function(x, filename, width=750, height=450, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.overlap.ms.heatmap.png")
#circos plot
png("clone.overlap.ms.circos.png",width=10,height=10,units="in",res=200)
par(las=1)
circos.par(start.degree=270,cell.padding=c(0,0,0,0),
           gap.after=c(rep(0,length(B6.mice)-1),15,rep(0,length(SLE.yaa.mice)-1),15))
chordDiagram(data.matrix(ms.overlap.df),annotationTrack = "grid", reduce=0,
             # grid.col=rep(sample(brewer.pal(12,"Paired")),length.out=length(c(B6.mice,SLE.yaa.mice))),#col=col_fun,
             preAllocateTracks = list(track.height = 0.6),
             order=c(B6.mice,SLE.yaa.mice))
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",niceFacing = TRUE, adj = c(0.5, 0))
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()
mtext("WT",side=2,line=-8)
mtext("SLE.yaa",side=4,line=-10)
dev.off()
#by size of shared clones
ms.overlap.size.df<-data.frame()
for(i in unique(clone.data.ab$mouse_ID)){
  mouse1.clones<-clone.data.ab[clone.data.ab$mouse_ID==i,]
  for (j in unique(clone.data.ab$mouse_ID)){
    if(i!=j){
      mouse2.clones<-clone.data.ab[clone.data.ab$mouse_ID==j,]
      ms.overlap.size.df[j,i]<-sum(mouse1.clones$cdr3 %in% ms.overlap.list[[paste0(i,"+",j)]])+
        sum(mouse2.clones$cdr3 %in% ms.overlap.list[[paste0(i,"+",j)]])
    }else{ms.overlap.size.df[j,i]<-NA}
  }
}
pheatmap(ms.overlap.size.df)
#better coloring
mat_breaks <- quantile_breaks(ms.overlap.size.df, n = 101)
colfunc <- colorRampPalette(c("#EDF8E9", "#006D2C"))
p1<-pheatmap(ms.overlap.size.df,annotation_col=annot.col,annotation_names_col=F,show_colnames = F,
             annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red")),
             color=colfunc(length(mat_breaks) - 1),breaks=mat_breaks)
save_pheatmap_png <- function(x, filename, width=750, height=450, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.overlap.size.ms.heatmap.png")
#by percent of clonotypes
ms.overlap.pct<-data.frame()
for(i in mice){
  mouse1.clones<-unique(clone.data.ab[clone.data.ab$mouse_ID==i,]$cdr3)
  for (j in mice){
    if(i!=j){
      mouse2.clones<-unique(clone.data.ab[clone.data.ab$mouse_ID==j,]$cdr3)
      pct.shared<-100*(length(ms.overlap.list[[paste0(i,"+",j)]])/
                         sum(length(mouse1.clones),length(mouse2.clones)))
      ms.overlap.pct[j,i]<-pct.shared
    }else{ms.overlap.pct[j,i]<-NA}
  }
}
p1<-pheatmap(ms.overlap.pct,annotation_col=annot.col,annotation_names_col=F,show_colnames = F,
             annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red")))
mat_breaks <- quantile_breaks(ms.overlap.pct, n = 101)
colfunc <- colorRampPalette(c("#EDF8E9", "#006D2C"))
p1<-pheatmap(ms.overlap.pct,annotation_col=annot.col,annotation_names_col=F,show_colnames = F,
             annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red")),
             color=colfunc(length(mat_breaks) - 1),breaks=mat_breaks)
save_pheatmap_png(p1, "clone.overlap.ms.pct.heatmap.png")
#by percent of cells
ms.overlap.size.pct.df<-data.frame()
for(i in unique(clone.data.ab$mouse_ID)){
  mouse1.clones<-clone.data.ab[clone.data.ab$mouse_ID==i,]
  for (j in unique(clone.data.ab$mouse_ID)){
    if(i!=j){
      mouse2.clones<-clone.data.ab[clone.data.ab$mouse_ID==j,]
      ms.overlap.size.pct.df[j,i]<-(sum(mouse1.clones$cdr3 %in% ms.overlap.list[[paste0(i,"+",j)]])+
        sum(mouse2.clones$cdr3 %in% ms.overlap.list[[paste0(i,"+",j)]]))/(length(mouse1.clones$cdr3)+length(mouse2.clones$cdr3))
    }else{ms.overlap.size.pct.df[j,i]<-NA}
  }
}
pheatmap(ms.overlap.size.pct.df)
#better coloring
mat_breaks <- quantile_breaks(ms.overlap.size.pct.df, n = 101)
colfunc <- colorRampPalette(c("#EDF8E9", "#006D2C"))
p1<-pheatmap(ms.overlap.size.pct.df,annotation_col=annot.col,annotation_names_col=F,show_colnames = F,
             annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red")),
             color=colfunc(length(mat_breaks) - 1),breaks=mat_breaks)
save_pheatmap_png <- function(x, filename, width=750, height=450, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.overlap.size.pct.ms.heatmap.png")

#public clones
public.clones<-unique(unlist(ms.overlap.list,use.names=F))
clone.data.ab.public<-clone.data.ab[clone.data.ab$cdr3 %in% public.clones,]
Counts <- table(as.character(clone.data.ab.public$mouse_ID),as.character(clone.data.ab.public$cdr3))
Proportions <- scale(Counts,scale=colSums(Counts),center=FALSE)*100 
pheatmap(as.data.frame.matrix(Proportions))
#add heatmap metadata
clone.meta.top<- as.data.table(clone.data.ab.public)[, lapply(.SD, data_concater), by=cdr3]
clone.meta.top$cdr3.freq <- as.numeric(as.character(clone.meta.top$cdr3.freq))
annot.col<-data.frame(cdr3=colnames(Proportions))
annot.col<-merge(annot.col,clone.meta.top[,c("cdr3","cdr3.freq","condition")],by="cdr3")
names(annot.col)[names(annot.col) == "condition"] <- "BMChimera"
names(annot.col)[names(annot.col) == "cdr3.freq"] <- "Clone Size"
annot.col$BMChimera[annot.col$BMChimera %in% "B6"]<-"WT"
annot.col$BMChimera[annot.col$BMChimera %in% "SLE.yaa"]<-"SLE.yaa"
annot.col$BMChimera[annot.col$BMChimera %in% "B6+SLE.yaa"]<-"Both"
row.names(annot.col) <- annot.col$cdr3
annot.col$cdr3 <- NULL
mat_breaks <- quantile_breaks(as.data.frame.matrix(Proportions), n = 101)
colfunc <- colorRampPalette(c("#EDF8E9", "#006D2C"))
p1<-pheatmap(as.data.frame.matrix(Proportions),annotation_col=annot.col,annotation_names_col=F,fontsize_col=5,
             # color=cividis(length(mat_breaks) - 1),breaks=mat_breaks,
             annotation_colors=list(BMChimera=c("WT"="grey","SLE.yaa"="red","Both"=alpha("red",0.5)))) #
save_pheatmap_png <- function(x, filename, width=1200, height=800, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.public.cluster.heatmap.png")

###CDR3 Network Analysis by condition
df<-clone.data.ab.ms_cdr3[clone.data.ab.ms_cdr3$cdr3.freq>1,]
# extract link between diag and patients ID for network
tmp1=data.frame(Diagnosis=df$condition,ID=df$mouse_ID,stringsAsFactors = F)
tmp1=tmp1[which(duplicated(paste(tmp1[,1],tmp1[,2],sep="_"))==F),]
colnames(tmp1)=c("x","y")
# extract link between clonotype_id.long and patients ID for network
tmp2=data.frame(df$mouse_ID,df$cdr3,stringsAsFactors = F)
tmp2=tmp2[which(duplicated(paste(tmp2[,1],tmp2[,2],sep="_"))==F),]
colnames(tmp2)=c("x","y")
# create table for qgraph
toQgraph=rbind(tmp1,tmp2)
# color of Dx
l=unique(c(toQgraph[,1],toQgraph[,2]))
col.tmp=rep(adjustcolor("grey85",.3),length(l))
col.tmp[which(l %in% c("B6"))]<-adjustcolor("black",.8)
col.tmp[which(l %in% c("SLE.yaa"))]<-adjustcolor("red",.8)
# color of subjects
col.tmp[which(l %in% unique(df$mouse_ID[which(df$condition=="B6")]))]<-adjustcolor("black",.3)
col.tmp[which(l %in% unique(df$mouse_ID[which(df$condition=="SLE.yaa")]))]<-adjustcolor("red",.3)
# size of nodes
size.tmp=rep(1,length(l))
size.tmp[which(l %in% c("B6","SLE.yaa"))]<-30
size.tmp[which(l %in% c(names(table(df$mouse_ID))))]<-15
for(i in which(size.tmp==1)){
  size.tmp[i]=gm_mean(df$cdr3.freq[which(df$cdr3==l[i])])
}
# color of edges
line.tmp=rep("grey",nrow(toQgraph))
line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("B6"))]<-"black"
line.tmp[which(c(toQgraph[,1],toQgraph[,2]) %in% c("SLE.yaa"))]<-"red"
labels.cex.tmp=l
labels.cex.tmp=rep(0.00000000001,length(l))
labels.cex.tmp[which(l %in% c("B6","SLE.yaa"))]<-1.5
labels.cex.tmp[which(l %in% c(names(table(df$mouse_ID))))]<-1
#graph
dim(toQgraph)
toQgraph$x[toQgraph$x %in% "B6"]<-"WT"
toQgraph$x[toQgraph$x %in% "SLE.yaa"]<-"SLE.yaa"
png("cdr3.network.by.mouse.png",width=5,height=5,units="in",res=400)
qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,edge.width=0.25,
       labels=TRUE,label.cex=labels.cex.tmp,directed=F,repulsion=1)
dev.off()

##between clusters
#calculate shared clones number (between clusters)
clust.overlap.df<-data.frame()
clust.overlap.list<-list() #list of lists (in case want to extract shared clone info)
for(i in unique(clone.data.ab.seurat$my.clusters2)){
  clust1.clones<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters2==i,]
  for (j in unique(clone.data.ab.seurat$my.clusters2)){
    if(i!=j){
      clust2.clones<-clone.data.ab.seurat[clone.data.ab.seurat$my.clusters2==j,]
      clust.overlap.list[[paste0(i,"+",j)]]<-intersect(clust1.clones$cdr3,clust2.clones$cdr3)
      clust.overlap.df[j,i]<-length(clust.overlap.list[[paste0(i,"+",j)]])
    }else{clust.overlap.df[j,i]<-NA}
  }
}
p1<-pheatmap(clust.overlap.df,angle_col=0)
cluster.names<-names(sort(table(clone.data.ab.seurat$my.clusters2),decreasing=T))
color.key<-data.frame(my.clusters2=cluster.names,color=hue_pal()(6))
annot.col<-data.frame(cluster=cluster.names)
rownames(annot.col)<-cluster.names
annot.colors<-list(cluster=hue_pal()(6))
names(annot.colors$cluster)<-cluster.names
#better coloring
mat_breaks <- quantile_breaks(clust.overlap.df, n = 101)
colfunc <- colorRampPalette(c("#EDF8E9", "#006D2C"))
p1<-pheatmap(clust.overlap.df,annotation_col=annot.col,annotation_names_col=F,show_colnames = F,
             annotation_colors=annot.colors,annotation_legend=F,
             color=colfunc(length(mat_breaks) - 1),breaks=mat_breaks)
save_pheatmap_png <- function(x, filename, width=700, height=500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "clone.overlap.clust.heatmap.png")
#circos plot
png("clone.overlap.clust.circos.png",width=10,height=10,units="in",res=200)
circos.par(start.degree=270,cell.padding=c(0,0,0,0))
chordDiagram(data.matrix(clust.overlap.df),annotationTrack = "grid",
             grid.col=hue_pal()(length(cluster.names)),#col=col_fun,
             preAllocateTracks = list(track.height = 0.6),
             order=cluster.names,symmetric = T)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) < 20) {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",niceFacing = TRUE, adj = c(0, 0.5))} 
  if(abs(xplot[2] - xplot[1]) > 20){
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",niceFacing = TRUE, adj = c(0.5, 0))}
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()
dev.off()

###CDR3 Network Analysis by cluster
df<-clone.data.ab.clust_cdr3[clone.data.ab.clust_cdr3$cdr3.freq>1,]
# extract link between cdr3 and cluster for network
tmp2=data.frame(df$my.clusters2,df$cdr3,stringsAsFactors = F)
tmp2=tmp2[which(duplicated(paste(tmp2[,1],tmp2[,2],sep="_"))==F),]
colnames(tmp2)=c("x","y")
toQgraph=tmp2
# color of clusters
l=unique(c(toQgraph[,1],toQgraph[,2]))
col.tmp=rep(adjustcolor("grey85",.3),length(l))
for(i in c(1:nrow(color.key))){
  col.tmp[which(l %in% color.key$my.clusters2[i])]<-adjustcolor(color.key$color[i],.8)
}
# size of nodes
size.tmp=rep(1,length(l))
size.tmp[which(l %in% cluster.names)]<-30
for(i in which(size.tmp==1)){
  size.tmp[i]=gm_mean(df$cdr3.freq[which(df$cdr3==l[i])])
}
# color of edges
line.tmp=rep("grey",nrow(toQgraph))
labels.cex.tmp=l
labels.cex.tmp=rep(0.00000000001,length(l))
labels.cex.tmp[which(l %in% cluster.names)]<-1
#graph
dim(toQgraph)
png("cdr3.network.by.cluster.png",width=5,height=5,units="in",res=400)
qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,edge.width=0.25,
       labels=TRUE,label.cex=labels.cex.tmp,directed=F,repulsion=1000)
dev.off()

###Stacked bar plots and pie charts
##making reference dfs
#clone frequency within condition
cdr3.list<-list()
v_gene.list<-list()
for(i in names(table(clone.data.ab$condition))){
  df<-clone.data.ab[clone.data.ab$condition==i,]
  cdr3.df<-as.data.frame(table(df$cdr3))
  colnames(cdr3.df)<-c("cdr3","cdr3.freq")
  cdr3.df$condition<-i
  cdr3.list[[i]]<-cdr3.df
  v_gene.df<-as.data.frame(table(df$v_gene))
  colnames(v_gene.df)<-c("v_gene","v_gene.freq")
  v_gene.df$condition<-i
  v_gene.df<-left_join(x = v_gene.df %>% group_by(v_gene) %>% mutate(id = row_number()),
                       y = clone.data.ab %>% select(v_gene,chain) %>% group_by(v_gene) %>% mutate(id = row_number()), 
                       by = c("v_gene","id"),keep=F)
  v_gene.list[[i]]<-v_gene.df
}
cdr3.data.bycondition<-bind_rows(cdr3.list)
cdr3.data.bycondition<-cdr3.data.bycondition%>% arrange(cdr3.freq)
v_gene.data.bycondition<-bind_rows(v_gene.list)
v_gene.data.bycondition<-v_gene.data.bycondition%>% arrange(v_gene.freq)
#clone frequency within mice
cdr3.list<-list()
v_gene.list<-list()
for(i in names(table(clone.data.ab$mouse_ID))){
  df<-clone.data.ab[clone.data.ab$mouse_ID==i,]
  cdr3.df<-as.data.frame(table(df$cdr3))
  colnames(cdr3.df)<-c("cdr3","cdr3.freq")
  cdr3.df$mouse_ID<-i
  cdr3.list[[i]]<-cdr3.df
  v_gene.df<-as.data.frame(table(df$v_gene))
  colnames(v_gene.df)<-c("v_gene","v_gene.freq")
  v_gene.df$mouse_ID<-i
  v_gene.df<-left_join(x = v_gene.df %>% group_by(v_gene) %>% mutate(id = row_number()),
                       y = clone.data.ab %>% select(v_gene,chain) %>% group_by(v_gene) %>% mutate(id = row_number()), 
                       by = c("v_gene","id"),keep=F)
  v_gene.list[[i]]<-v_gene.df
}
cdr3.data.bymouse<-bind_rows(cdr3.list)
cdr3.data.bymouse<-cdr3.data.bymouse%>% arrange(cdr3.freq)
cdr3.data.bymouse<-left_join(x = cdr3.data.bymouse, y = metadata, by = "mouse_ID",keep=F)
v_gene.data.bymouse<-bind_rows(v_gene.list)
v_gene.data.bymouse<-v_gene.data.bymouse%>% arrange(v_gene.freq)
v_gene.data.bymouse<-left_join(x = v_gene.data.bymouse, y = metadata, by = "mouse_ID",keep=F)
#clone frequency within condition + cluster
cdr3.list<-list()
v_gene.list<-list()
for(i in names(table(clone.data.ab.seurat$condition))){
  df2<-clone.data.ab.seurat[clone.data.ab.seurat$condition==i,]
  for(j in names(table(df2$my.clusters))){
    df3<-df2[df2$my.clusters==j,]
    cdr3.df<-as.data.frame(table(df3$cdr3))
    colnames(cdr3.df)<-c("cdr3","cdr3.freq")
    cdr3.df$condition<-i
    cdr3.df$my.clusters<-j
    cdr3.list[[paste0(i,".",j)]]<-cdr3.df
    v_gene.df<-as.data.frame(table(df3$v_gene))
    colnames(v_gene.df)<-c("v_gene","v_gene.freq")
    v_gene.df$condition<-i
    v_gene.df$my.clusters<-j
    v_gene.df<-left_join(x = v_gene.df %>% group_by(v_gene) %>% mutate(id = row_number()),
                         y = clone.data.ab %>% select(v_gene,chain) %>% group_by(v_gene) %>% mutate(id = row_number()), 
                         by = c("v_gene","id"),keep=F)
    v_gene.list[[paste0(i,".",j)]]<-v_gene.df
  }
}
cdr3.data.bycluster<-bind_rows(cdr3.list)
cdr3.data.bycluster<-cdr3.data.bycluster%>% arrange(cdr3.freq)
v_gene.data.bycluster<-bind_rows(v_gene.list)
v_gene.data.bycluster<-v_gene.data.bycluster%>% arrange(v_gene.freq)
save(list=c(ls(pattern="data")),file="clone.freq.prop.RData")

#by cluster
df<-cdr3.data.bycluster
png("cdr3.bycluster.pie.png",width=15,height=5.5,units="in",res=400)
# par(mfrow=c(2,length(names(table(df$my.clusters)))),mai=c(1,.1,.25,.1) )
layout(matrix(c(1:12,rep(13,6)),ncol=6,byrow=T),heights=c(4,4,1))
par(mai=c(.1,.1,.1,.1) )
for(i in names(table(df$condition))){
  for(k in names(table(df$my.clusters))){
    df2<-df[df$condition==i&df$my.clusters==k,]
    # df2<-df2[df2$cdr3%in% names(head(sort(table(df2$cdr3),decreasing=T),n=100)),] #select only top 100 clones
    df2$cdr3.prop2<- df2$cdr3.freq/sum(df2$cdr3.freq)
    df2.lbl<-df2[df2$cdr3.prop2>0.02,] #which clones to label
    if(length(df2.lbl$cdr3)!=0){
      df2.lbl$lbl<-paste0(df2.lbl$cdr3, "\n", round(100*df2.lbl$cdr3.prop2, 1),"%")
      df2<-left_join(df2,df2.lbl)
    } else {df2$lbl<-NA}
    # pie(df2$cdr3.freq,labels=df2$lbl,col=sample(brewer.pal(12,"Paired")),border=NA, xlab=j, cex=.4,line=-1)
    df2$color<-"darkgrey"
    for(j in 1:length(df2$cdr3.freq)){
      if(between(df2$cdr3.freq[j],2,2)){df2$color[j]<-"grey30" }
      if(between(df2$cdr3.freq[j],3,3)){df2$color[j]<-"tan1" }
      if(between(df2$cdr3.freq[j],4,4)){df2$color[j]<-"blue4"}
      if(between(df2$cdr3.freq[j],5,9)){df2$color[j]<-"red2"}
      if(df2$cdr3.freq[j]>=10){df2$color[j]<-"forestgreen"}
    }
    # pie(df2$cdr3.freq,labels=df2$lbl,col=sample(brewer.pal(12,"Paired")),border=NA, xlab=j, cex=0.4,line=-1)
    pie(df2$cdr3.freq,col=df2$color,border=NA,labels=NA)#,labels=df2$lbl,cex=0.3)
    # mtext(k,side=1,line=-3,cex=0.75)
    mtext(paste0("n = ",sum(df2$cdr3.freq)),side=1,line=-2,cex=1)
  }
  # mtext(i,side=3,cex=1.5,line=-1)
}
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=6,legend=c("Unique","\u2265 2","\u2265 3","\u2265 4","\u2265 5","\u2265 10"),
       fill=c("darkgrey","grey30","tan1","blue4","red2","forestgreen"),bty="n",cex=2)
dev.off()

save.image("clone.analysis.paired.RData")

