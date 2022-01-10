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
library(patchwork)
library(vegan)
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

load("clone.ab.data.RData") 

###GLIPH export
dir.create("./GLIPH2")
setwd("./GLIPH2")
df<-clone.data
df.TRA<-df[df$chain=="TRA",]
df.TRB<-df[df$chain=="TRB",]
cells<-inner_join(df.TRB[,c("barcode","cdr3","v_gene","j_gene")],
                  df.TRA[,c("barcode","cdr3","v_gene","j_gene")],by="barcode")
cells<-left_join(cells,unique(clone.data[,c("barcode","mouse_ID","condition","my.clusters")]),by="barcode")
cells$cdr3_ID<-paste0(cells$cdr3.x,"+",cells$cdr3.y)
cdr3s<-as.data.table(cells)[, lapply(.SD, data_concater), by=cdr3_ID]
cells$clone_ID<-paste0(cells$cdr3.x,cells$v_gene.x,cells$j_gene.x,
                       cells$cdr3.y,cells$v_gene.y,cells$j_gene.y,cells$mouse_ID)
cells.cloneID<-transform(cells,freq=ave(seq(nrow(cells)),clone_ID,FUN=length))
cells.cloneID<-as.data.table(cells.cloneID)[, lapply(.SD, data_concater), by=clone_ID]
cells.cloneID$freq<-as.numeric(as.character(cells.cloneID$freq))
cells.GLIPH2<-data.frame(CDR3b=cells.cloneID$cdr3.x,TRBV=cells.cloneID$v_gene.x,TRBJ=cells.cloneID$j_gene.x,
                         CDR3a=cells.cloneID$cdr3.y,subject.condition=paste0(cells.cloneID$mouse_ID,".",cells.cloneID$condition),
                         count=cells.cloneID$freq)
write.table(cells.GLIPH2,file="cdr3.paired.GLIPH2.txt",sep="\t",row.names=F,col.names=F,quote=F)
save(list=c("cells","cdr3s","cells.GLIPH2","cells.cloneID"),file="meta.for.GLIPH.RData")
setwd("../")
save.image("GLIPH.submit.RData")
#submit file to: http://50.255.35.37:8080/

#load GLIPH data
GLIPH.df<-read.csv("./GLIPH2/P2769.csv",header=T) # ./GLIPH2/P784.csv
GLIPH.df$index<-paste0("GLIPH_",GLIPH.df$index)
GLIPH.groups<-unique(GLIPH.df[,1:12])
#add metadata
GLIPH.df$cdr3_ID<-paste0(GLIPH.df$TcRb,"+",GLIPH.df$TcRa)
GLIPH.df.filter<-GLIPH.df[GLIPH.df$number_subject>1&
                            GLIPH.df$vb_score<0.05&GLIPH.df$final_score<1e-5,]
GLIPH.filter.groups<-unique(GLIPH.df.filter[,1:12])
head(GLIPH.filter.groups[order(GLIPH.filter.groups$final_score),],n=9)
save.image("GLIPH.results.RData")

#Add Seurat to GLIPH data
#adding barcode to GLIPH df
df<-clone.data
cells<-inner_join(df[df$chain=="TRB",c("barcode","cdr3")],df[df$chain=="TRA",c("barcode","cdr3")],by="barcode")
cells$cdr3_ID<-paste0(cells$cdr3.x,"+",cells$cdr3.y)
cells.meta<-left_join_replace(cells,clone.data.ab,by="barcode")
GLIPH.barcode.ref<-inner_join(cells[,c("barcode","cdr3_ID")],unique(GLIPH.df.filter[,c(1:19,31)]),by="cdr3_ID") #17
#adding Seurat using barcode
GLIPH.data<-left_join_replace(x = GLIPH.barcode.ref,y = cells.meta, by = "barcode")
GLIPH.data<-left_join_replace(GLIPH.data,metadata,by="mouse_ID")
for(i in 1:nrow(GLIPH.filter.groups)){
  index<-GLIPH.filter.groups$index[i]
  select.index<-GLIPH.data[GLIPH.data$index==index,]
  select.index.seurat<-select.index[!is.na(select.index$my.clusters),]
  GLIPH.filter.groups$total.size[i]<-nrow(select.index)
  GLIPH.filter.groups$tfr.score[i]<-100*nrow(select.index.seurat[select.index.seurat$my.clusters=="3",])/nrow(select.index.seurat)
  GLIPH.filter.groups$tcm.score[i]<-100*nrow(select.index.seurat[select.index.seurat$my.clusters=="0",])/nrow(select.index.seurat)
  GLIPH.filter.groups$top.cluster[i]<-names(sort(table(select.index.seurat$my.clusters2),decreasing=T))[1]
  GLIPH.filter.groups$gm.clone.size[i]<-gm_mean(table(select.index.seurat$cdr3_ID))
  GLIPH.filter.groups$clone.div.shannon[i]<-diversity(table(select.index.seurat$cdr3_ID),index="shannon")
  GLIPH.filter.groups$clone.div.simpson[i]<-diversity(table(select.index.seurat$cdr3_ID),index="simpson")
  GLIPH.filter.groups$clone.div.invsimp[i]<-diversity(table(select.index.seurat$cdr3_ID),index="invsimp")
}
GLIPH.df.filter<-left_join(GLIPH.df.filter,GLIPH.filter.groups[,c("index","total.size","tfr.score","tcm.score","top.cluster","gm.clone.size",
                                                                  "clone.div.shannon","clone.div.simpson","clone.div.invsimp")],by="index")
GLIPH.data<-left_join_replace(GLIPH.data,GLIPH.filter.groups[,c("index","total.size","tfr.score","tcm.score","top.cluster","gm.clone.size",
                                                                "clone.div.shannon","clone.div.simpson","clone.div.invsimp")],by="index")
GLIPH.data.seurat<-GLIPH.data[!is.na(GLIPH.data$my.clusters),]
GLIPH.barcode.ref<-left_join(GLIPH.barcode.ref,GLIPH.filter.groups[,c("index","total.size","tfr.score","tcm.score","top.cluster","gm.clone.size",
                                                                      "clone.div.shannon","clone.div.simpson","clone.div.invsimp")],by="index")
clusters<-names(table(GLIPH.data.seurat$my.clusters))
mice<-names(table(metadata$mouse_ID))
B6.mice<-names(table(metadata[metadata$condition=="B6",]$mouse_ID))
SLE.yaa.mice<-names(table(metadata[metadata$condition=="SLE.yaa",]$mouse_ID))
save(list=c(ls(pattern="GLIPH."),ls(pattern="mice"),ls(pattern="cells"),ls(pattern="clone"),"metadata","clusters"),
     file="GLIPH.data.RData")
save.image("GLIPH.data2.RData")

##Venn diagram
#between conditions
df<-GLIPH.data
dz.list<-list()
for(i in unique(df$condition)){
  dz.list[[i]]<-unique(df$index[df$condition==i])
}
venn.diagram(
  x = dz.list,category.names = c("B6","SLE.yaa"),filename = "condition.GLIPH.venn.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("black", 'red'),fill=c(alpha("black",0.3), alpha('red',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(dz.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.6,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("black", 'red'),cat.dist = c(0.1, 0.1),
  ext.pos=180,ext.line.lwd=0.25,ext.dist=-0.2
)

#NMDS
df<-GLIPH.data
community_matrix<-as.data.frame.matrix(table(df$mouse_ID,df$index))
example_NMDS=metaMDS(community_matrix, k=2,trymax=100)#"The stress, or the disagreement between 2-D configuration and predicted values from the regression"
#A good rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation
png("gliph.NMDS.stress.png",width=4,height=4,units="in",res=200)
stressplot(example_NMDS)
dev.off()
plot(example_NMDS)
treat=rep("B6",nrow(community_matrix))
treat[rownames(community_matrix) %in% SLE.yaa.mice]<-"SLE.yaa"
colors=rep("black",nrow(community_matrix))
colors[rownames(community_matrix) %in% SLE.yaa.mice]<-"red"
elevation=runif(10,0.5,1.5)
png("gliph.NMDS.png",width=6,height=6,units="in",res=200)
ordiplot(example_NMDS,type="n")
ordiellipse(example_NMDS$point[grep("B6",treat),],draw="polygon",groups=treat[treat=="B6"],col="grey",alpha=50)
ordiellipse(example_NMDS$point[grep("SLE.yaa",treat),],draw="polygon",groups=treat[treat=="SLE.yaa"],col="red",alpha=50)
orditorp(example_NMDS,display="species",col="black",pch=21,bg=adjustcolor("grey",0.5),pcex=1,air=50000,cex=0.1)
orditorp(example_NMDS,display="sites",col=colors,air=0.0001,cex=1.25)
dev.off()

#public GLIPH repeteroire, by condition
df<-GLIPH.data
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
pr.res.GLIPH<-as.data.frame(left_join(pr.res.GLIPH,GLIPH.filter.groups[,c("index","final_score","number_subject",
                                                            "number_unique_cdr3","total.size","tfr.score","tcm.score",
                                                            "top.cluster","gm.clone.size","clone.div.shannon",
                                                            "clone.div.simpson","clone.div.invsimp")],
                        by=c("CDR3.aa"="index")))
labels<-pr.res.GLIPH[abs(pr.res.GLIPH$log2FC)>1,]
ggplot(pr.res.GLIPH,aes(x = sum.SLE.yaa, y =  sum.B6,size=number_subject))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(alpha=0.25,aes(colour=number_unique_cdr3))+ #tfr.score
  guides(fill = guide_legend(override.aes = list(size = 7)))+#,color=F)+ #,size=F
  labs(x = "SLE.yaa", y = "B6", size="Samples",color="Clones")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+
  geom_label_repel(data = labels,aes(label = CDR3.aa),size = 3,alpha=0.7,
                   color="royalblue",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)+
  scale_color_gradient(low="grey",high="darkgreen",trans="log",breaks=c(1,5,10,15,25),labels=format(c(1,5,10,15,25)))
ggsave2("overlap.B6.SLE.yaa.GLIPH.scatter.all.labeled.png",width=12, height=8,device="png")

### B6 vs 564 specific expanded GLIPHs (log2FC, 2,3)
B6.public.expand.GLIPH<-pr.res.GLIPH[pr.res.GLIPH$log2FC< -3 &pr.res.GLIPH$sum.B6>2,]#[(pr.res.cdr3$avgfreq.B6-10)>10*pr.res.cdr3$avgfreq.SLE.yaa,]
SLE.yaa.public.expand.GLIPH<-pr.res.GLIPH[pr.res.GLIPH$log2FC>3&pr.res.GLIPH$sum.SLE.yaa>2,]#pr.res.cdr3[(pr.res.cdr3$avgfreq.SLE.yaa-10)>10*pr.res.cdr3$avgfreq.B6,]
pr.res.GLIPH$annotate<-ifelse(pr.res.GLIPH$CDR3.aa %in% B6.public.expand.GLIPH$CDR3.aa, "B6-expanded","Public")
pr.res.GLIPH$annotate[pr.res.GLIPH$CDR3.aa %in% SLE.yaa.public.expand.GLIPH$CDR3.aa]<-"SLE.yaa-expanded"
pr.res.GLIPH$annotate <- factor(pr.res.GLIPH$annotate, levels = c("Public","B6-expanded","SLE.yaa-expanded"))
ggplot(pr.res.GLIPH,aes(x = sum.SLE.yaa, y =  sum.B6,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(4, 10))+
  geom_point(aes(colour=factor(annotate),fill = factor(annotate)), shape=21) +
  scale_color_manual(values = c(alpha("black",0.5),"black","red"))+
  scale_fill_manual(values = c("#1C00ff00",alpha("black",0.2),alpha("red",0.2)))+
  guides(fill = guide_legend(override.aes = list(size = 7)),color=F,size=F)+ #,size=F
  labs(x = "SLE.yaa", y = "B6", size="Samples",fill="Expansion")+
  theme(legend.direction = "vertical", legend.box = "horizontal")#+#theme(legend.title = element_blank())+
ggsave2("GLIPH.expansion.scatter.png",width=6, height=4,device="png")


###GLIPH Network Analysis by condition
df<-GLIPH.data
# extract link between diag and patients ID for network
tmp1=data.frame(Diagnosis=df$condition,ID=df$mouse_ID,stringsAsFactors = F)
tmp1=tmp1[which(duplicated(paste(tmp1[,1],tmp1[,2],sep="_"))==F),]
colnames(tmp1)=c("x","y")
# extract link between clonotype_id.long and patients ID for network
tmp2=data.frame(df$mouse_ID,df$index,stringsAsFactors = F)
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
size.tmp[which(l %in% c("B6","SLE.yaa","AID","m564"))]<-30
size.tmp[which(l %in% c(names(table(df$mouse_ID))))]<-15
for(i in which(size.tmp==1)){
  size.tmp[i]=gm_mean(df$total.size[which(df$index==l[i])])/10+1
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
png("gliph.network.by.mouse.png",width=5,height=5,units="in",res=400)
qgraph(toQgraph,color=col.tmp,vsize=size.tmp/2,edge.color=line.tmp,edge.width=0.25,
       labels=TRUE,label.cex=labels.cex.tmp,directed=F,repulsion=70)
dev.off()



### GLIPH-tcr network
diff.GLIPH<-head(GLIPH.filter.groups$index[order(-GLIPH.filter.groups$total.size)],n=30)
df<-GLIPH.data[GLIPH.data$index %in% c(diff.GLIPH) ,]#,exp.GLIPH&GLIPH.data$cdr3.freq>1
tmp3=data.frame(df$index,df$cdr3,stringsAsFactors = F)
tmp3=tmp3[which(duplicated(paste(tmp3[,1],tmp3[,2],sep="_"))==F),]
colnames(tmp3)=c("x","y")
# create table for qgraph
toQgraph=rbind(tmp3)
# color of Dx
l=unique(c(toQgraph[,1],toQgraph[,2]))
col.tmp=rep(adjustcolor("blue",.4),length(l)) #adjustcolor("purple",.3)
col.tmp[which(l %in% unique(df$cdr3[which(df$condition=="B6")]))]<-adjustcolor("black",.2)
col.tmp[which(l %in% unique(df$cdr3[which(df$condition=="SLE.yaa")]))]<-adjustcolor("red",.2)
col.tmp[which(l %in% unique(intersect(df$cdr3[which(df$condition=="B6")],
                                      df$cdr3[which(df$condition=="SLE.yaa")])))]<-adjustcolor("forestgreen",.9)
# size of nodes
size.tmp=rep(1,length(l))
size.tmp[which(l %in% unique(df$index))]<-10
for(i in which(size.tmp==1)){
  size.tmp[i]=log(gm_mean(df$cdr3.freq[which(df$cdr3==l[i])]))*5+1
}
# color of edges
line.tmp=rep("grey",nrow(toQgraph))
#labels
labels.cex.tmp=rep(0.0001,length(l))
labels.cex.tmp[which(l %in% c(GLIPH.data$CDR3.aa))]<-1 #,exp.GLIPH
#graph
dim(toQgraph)
png("diff.GLIPH.network.by.cdr3s.png",width=10,height=10,units="in",res=400)
qgraph(toQgraph,color=col.tmp,vsize=size.tmp,edge.color=line.tmp,edge.width=0.15,
       directed=F,repulsion=1,borders=F,labels=F)
dev.off()


load("graphed.RData")
T_all.combined<-SLE.obj.combined
#Add GLIPH to Seurat data
GLIPH.barcode.ref.ordered<-GLIPH.barcode.ref[order(GLIPH.barcode.ref$final_score),]
meta<-GLIPH.barcode.ref.ordered[!duplicated(GLIPH.barcode.ref.ordered$barcode),]
rownames(meta)<-meta$barcode
T_all.combined@meta.data$orig.barcode<-rownames(T_all.combined@meta.data)
T_all.combined<-AddMetaData(T_all.combined,cbind(meta[,1:14],dplyr::select(meta, total.size)))
Idents(T_all.combined) <- "condition"
SLE.yaa.combined<-subset(T_all.combined, idents='SLE.yaa')
B6.combined<-subset(T_all.combined, idents='B6')

save.image("GLIPH.analyzed.RData")
