rm(list=ls())
library(pheatmap)
library(reticulate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(Seurat)

meta.list<-list()
for(i in c("B6","SLE.yaa")){
  file.names<-list.files(path=paste0("../scTfh.data/",i))
  if(i=="B6"){cond<-"B6"}else{cond<-"SLE.yaa"}
  meta.list[[i]]<-data.frame(condition2=cond,file.name=file.names)
}
meta<-do.call("rbind",meta.list)

#ML feature heatmap
df<-read.csv("unsupervised.rep.features.csv",header=T,row.names=1)
deep.feats<-t(df)
pheatmap(deep.feats)
annot.col<-data.frame(file.name=colnames(deep.feats))
annot.col<-merge(annot.col,meta,by="file.name")
names(annot.col)[names(annot.col) == "condition2"] <- "BMChimera"
row.names(annot.col) <- annot.col$file.name
annot.col$file.name <- NULL
p<-pheatmap(deep.feats,annotation_col=annot.col,fontsize_col=3,annotation_names_col=F,show_colnames=F,show_rownames=F,
         annotation_colors=list(BMChimera=c("B6"="grey","SLE.yaa"="red")),main="DeepTCR")
save_pheatmap_png <- function(x, filename, width=700, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p,"unsupervised.rep.sample.heatmap.png")

