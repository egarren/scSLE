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
library(yingtools2)
library(qgraph)
library(GGally)
library(network)
library(sna)
library(tools)
library(ggalluvial)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(circlize)
library(useful)
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
gm_mean = function(x, na.rm=TRUE){  #geometric mean functions for graphing
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 
gm_ymax<- function(x, na.rm = TRUE)
{
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))+exp(sd(log(x), na.rm = na.rm))
}

## export masterdb to GLIPH
load("masterdb3.RData")
df<-unique(masterdb)
df[df==""]<-NA
masterdb.GLIPH<-unique(data.frame(CDR3b=df$cdr3.b,TRBV=df$TRBV,TRBJ=df$TRBJ,CDR3a="NA",
                                  # subject.condition=paste0(df$disease),#,".",df$antigen,".",df$epitope),
                                  subject.condition=paste0(df$study2,":",df$disease2),
                                  count=df$freq))
write.table(masterdb.GLIPH,file="./GLIPH2/ref.GLIPH.txt",sep="\t",row.names=F,col.names=F,quote=F)
#submit file to: http://50.255.35.37:8080/

#load GLIPH ref data
GLIPH.ref.df<-read.csv("./GLIPH2/P829.csv",header=T)
GLIPH.ref.df$index<-paste0("GLIPH.ref_",GLIPH.ref.df$index)
GLIPH.ref.groups<-unique(GLIPH.ref.df[,1:12])
GLIPH.ref.df.filter<-GLIPH.ref.df#[GLIPH.ref.df$final_score<1e-7,]
GLIPH.ref.filter.groups<-unique(GLIPH.ref.df.filter[,1:12])
head(GLIPH.ref.filter.groups[order(GLIPH.ref.filter.groups$final_score),],n=9)

##add back masterdb metadata
df<-unique(masterdb) #P735
df[df==""]<-NA
df$Sample_cdr3<-paste0(df$study2,":",df$disease2,df$cdr3.b,df$TRBV)
GLIPH.ref.df.filter$Sample_cdr3<-paste0(GLIPH.ref.df.filter$Sample,GLIPH.ref.df.filter$TcRb,GLIPH.ref.df.filter$V)
GLIPH.ref.df.filter.meta<-left_join(GLIPH.ref.df.filter[,c(1:19,31)],df,by="Sample_cdr3")
GLIPH.ref<-unique(GLIPH.ref.df.filter.meta[,c("pattern","type","TcRb","disease","antigen",
                                              "epitope","disease2","study2","final_score","Species")])


load("GLIPH.data.RData")
##making public repertoire df
df<-GLIPH.data
length(unique(intersect(GLIPH.ref.df.filter$pattern,GLIPH.data$pattern)))
length(unique(GLIPH.data$pattern))
imm.list<-list()
for(i in names(table(df$mouse_ID))){ #creating immunarch list from collapsed data
  df2<-df[df$mouse_ID==i&!is.na(df$index),]
  df2<-transform(df2,Clones=ave(seq(nrow(df2)),index,FUN=length))
  df2<-unique(df2[,c("index","Clones")])
  df2$Proportion<-df2$Clones/sum(df2$Clones)
  imm.list[[paste0(i,"_")]]<-tibble(Clones=df2$Clones, Proportion=df2$Proportion, CDR3.aa=df2$index )
}
immdata = repLoad("./rep")
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
pr.B6 = pubRepFilter(pr, immdata$meta, c(condition = "B6"))
pr.SLE.yaa = pubRepFilter(pr, immdata$meta, c(condition = "SLE.yaa"))
pr.B6[is.na(pr.B6)]<-0
pr.SLE.yaa[is.na(pr.SLE.yaa)]<-0
pr.B6[["avgfreq.B6"]] = rowMeans(public_matrix(pr.B6), na.rm = T)
pr.SLE.yaa[["avgfreq.SLE.yaa"]] = rowMeans(public_matrix(pr.SLE.yaa), na.rm = T)
pr.B6[["sum.B6"]] = rowSums(public_matrix(pr.B6)[,1:2], na.rm = T)
pr.SLE.yaa[["sum.SLE.yaa"]] = rowSums(public_matrix(pr.SLE.yaa)[,1:2], na.rm = T)
pr.res.GLIPH = as.data.frame(dplyr::full_join(pr.B6, pr.SLE.yaa, by = c("CDR3.aa"))) #,"J.name"
pr.res.GLIPH[is.na(pr.res.GLIPH)]<-0
pr.res.GLIPH[["Samples.sum"]] = pr.res.GLIPH[["Samples.x"]] + pr.res.GLIPH[["Samples.y"]]
pr.res.GLIPH[["freq.ratio"]] = apply(pr.res.GLIPH[, c("avgfreq.B6", "avgfreq.SLE.yaa")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.GLIPH[["log2FC"]]= apply(pr.res.GLIPH[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log2((x[2])/(x[1])))
pr.res.GLIPH<-as.data.frame(left_join(pr.res.GLIPH,GLIPH.filter.groups[,c("index","final_score","number_subject","pattern",
                                                            "number_unique_cdr3","total.size","tfr.score")],
                        by=c("CDR3.aa"="index")))
labels<-pr.res.GLIPH[pr.res.GLIPH$sum.SLE.yaa>10|pr.res.GLIPH$sum.B6>10|(pr.res.GLIPH$sum.SLE.yaa>5&pr.res.GLIPH$sum.B6>5),]#,]
ggplot(pr.res.GLIPH,aes(x = sum.SLE.yaa, y =  sum.B6,size=number_subject))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(alpha=0.25,aes(colour=number_unique_cdr3))+ #tfr.score
  guides(fill = guide_legend(override.aes = list(size = 7)))+#,color=F)+ #,size=F
  labs(x = "SLE.yaa", y = "WT", size="Samples",color="Clones")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+
  geom_label_repel(data = labels,aes(label = pattern),size = 3,alpha=0.7,
                   color="royalblue",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)+
  scale_color_gradient(low="grey",high="darkgreen",trans="log",breaks=c(5,25,100,400),labels=format(c(5,25,100,400)))
ggsave2("overlap.B6.SLE.yaa.GLIPH.scatter.indexlabel.png",width=6, height=4,device="png")
pr.res.GLIPH.db<-pr.res.GLIPH

#highlighting GLIPHs corresponding to specific antigen
df.db<-unique(GLIPH.ref[GLIPH.ref$antigen!="unknown"&!is.na(GLIPH.ref$antigen),]) #&GLIPH.ref$disease2!="species"
df.db<-df.db[order(df.db$final_score),]
pr.res.db<-unique(join(pr.res.GLIPH, df.db[,!names(df.db)%in%c("final_score")], by="pattern", type="left", match="first"))
match.df<-unique(pr.res.db[,c("pattern","disease","antigen","epitope","disease2")])
match.freq<-1-sum(is.na(match.df$disease))/nrow(match.df)
pr.res.db$disease2[is.na(pr.res.db$disease2)]<-"unknown"
brewer.pal(7,"Dark2") #display.brewer.pal
ggplot(pr.res.db,aes(x = sum.SLE.yaa, y =  sum.B6,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(aes(colour=factor(disease2),fill = factor(disease2)), shape=21) + 
  scale_color_manual(values = c("#1B9E77","#7570B3",alpha("black",0.5),"#A6761D"))+
  scale_fill_manual(values = c(alpha("#1B9E77",0.5),alpha("#7570B3",0.5),"#1C00ff00",alpha("#A6761D",0.5)))+
  guides(fill = guide_legend(override.aes = list(size = 5)),color=F)+ #,size=F
  labs(x = "SLE.yaa", y = "B6", size="Samples",fill="Disease")+
  theme(legend.direction = "vertical", legend.box = "horizontal")#+#theme(legend.title = element_blank())+
ggsave2("annot.index.GLIPH.scatter.bycondition.png",width=6.5, height=4,device="png")

#UMAP for specific GLIPH
df<-T_all.combined
df@meta.data<-unique(plyr::join(df@meta.data, df.db[,!names(df.db)%in%c("final_score")], by="pattern", type="left", match="first"))
df@meta.data$disease2[is.na(df@meta.data$disease2)]<-"unknown"
for(j in c("disease","antigen","epitope","pattern")){
  Idents(df)<-j
  cell.list<-list()
  for(i in names(head(sort(table(df[[j]]), decreasing = T), n=9))){
    cell.list[[i]]<-WhichCells(df, idents = i)
  }
  DimPlot(df, reduction = "umap", split.by = "condition",order=T, #"pattern","type","TcRb","disease","antigen","epitope","disease2","study2","final_score","Species"
          cells.highlight=cell.list,cols.highlight=brewer.pal(n = length(cell.list), name = "Paired"),
          sizes.highlight=1)+ #NoLegend()+cols="grey",
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
  ggsave2(paste0("GLIPH.",j,".umap.png"),width=8, height=3.5,device="png")
}

