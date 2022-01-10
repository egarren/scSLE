rm(list=ls())
library(dplyr)
library(tidyverse)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(umap)
library(gridExtra)
library(data.table)
library(ggpubr)
library(clusterProfiler)
library(msigdbr)
library(Seurat)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(scales)
data_concater2 <- function(x){
  x<- names(sort(table(x),decreasing=T))
  paste(x[1])
}
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}
rxn.consistency<- function(df,min.range=1e-3) {
  df = -log(df+1)
  df<-df[apply(df, MARGIN =  1, FUN = max, na.rm = T) - apply(df, MARGIN =  1, FUN = min, na.rm = T) >=min.range,]
  df=df-min(apply(df, MARGIN =  1, FUN = min, na.rm = T))
  return(df)
}
cohens_d<- function(x,y) {
  pooled_std = sqrt(((length(x)-1) * var(x) + (length(y)-1) * var(y)) / (length(x) + length(y) - 2))
  return ((mean(x) -mean(y)) / pooled_std)
}
my.wilcox<-function(df,groupA,groupB){
  results<-matrix(ncol=3,nrow=nrow(df))
  for(i in 1:nrow(df)){
    groupA.val<-as.numeric(df[i,groupA])
    groupB.val<-as.numeric(df[i,groupB])
    obj<-try(wilcox.test(groupA.val,groupB.val,alternative="two.sided"))
    if (!is(obj, "try-error")) {
      results[i,1]=obj$statistic
      results[i,2]=obj$p.value
      results[i,3]=cohens_d(groupA.val,groupB.val)
    } 
  }
  results<-data.frame(results)
  rownames(results)<-rownames(df)
  colnames(results)<-c("wilcox_stat","wilcox_pval","cohens_d")
  results$adjusted_pval=p.adjust(results$wilcox_pval,method="BH")
  return(results)
}
amino_acid_metab = c("Alanine and aspartate metabolism",
                    "Arginine and Proline Metabolism",
                    "beta-Alanine metabolism",
                    "Cysteine Metabolism",
                    "D-alanine metabolism",
                    "Folate metabolism",
                    "Glutamate metabolism",
                    "Glycine, serine, alanine and threonine metabolism",
                    "Histidine metabolism",
                    "Lysine metabolism",
                    "Methionine and cysteine metabolism",
                    "Taurine and hypotaurine metabolism",
                    "Tryptophan metabolism",
                    "Tyrosine metabolism",
                    "Urea cycle",
                    "Valine, leucine, and isoleucine metabolism")


#micropool analysis
meta<-read.csv("../compass.meta.csv",row.names=1)
cell.meta<-meta #if no micropooling
rownames(cell.meta)<-gsub("-",".",cell.meta$cell_id)

#reaction info
a="metarxn" # or ""
reaction.penalties<-read.csv("reactions.tsv",sep="\t",row.names=1,na.strings=c("NA","N/A"))
reaction.meta<-read.csv("reaction_metadata.csv",na.strings=c("NA","N/A"))
rxn.names<-read_csv("rxn.names.csv")
rxn.names<-unique(rxn.names) #get rid duplicates
rxn.names <- rxn.names[match(unique(rxn.names$rxn), rxn.names$rxn),] #only keep first occurance of remaining duplicates
if(a=="metarxn"){
  metarxns<-read.csv("metareactions.csv",row.names=1,na.strings=c("NA","N/A")) #if using metareactions
  metarxn.map<-read.csv("metarxn.map.csv",row.names=1,na.strings=c("NA","N/A")) #if using metareactions
  metarxn.map$rxn<-rownames(metarxn.map)
  for(i in 1:nrow(metarxn.map)){
    r=metarxn.map$rxn[i]
    if(r %in% reaction.meta$reaction_no_direction){
      metarxn.map$metadata_r_id[i]<-r
    }else if(substr(r,1,nchar(r)-4) %in% reaction.meta$reaction_no_direction){
      metarxn.map$metadata_r_id[i]<-substr(r,1,nchar(r)-4) 
    }
  }
  metarxn.map<-left_join(metarxn.map,reaction.meta,by=c("metadata_r_id"="reaction_no_direction"))
  #normalization
  reaction.consistencies<-rxn.consistency(metarxns)
}else{
  reaction.consistencies<-rxn.consistency(reaction.penalties)
}


#calculate reaction consistency and wilcox test
for(j in c("all",unique(cell.meta$cell_type))){
  if(j=="all"){
    group1<-cell.meta$condition=="SLE.yaa" 
    group2<-cell.meta$condition=="B6"
  }else{
    group1<-cell.meta$condition=="SLE.yaa" & cell.meta$cell_type==j
    group2<-cell.meta$condition=="B6" & cell.meta$cell_type==j
  }
  wilcox.results<-my.wilcox(reaction.consistencies,rownames(cell.meta)[group1],rownames(cell.meta)[group2])
  #add in reaction metadata
  if(a=="metarxn"){
    wilcox.results$metarxn_id<-as.integer(rownames(wilcox.results))
    wilcox.results<-left_join(wilcox.results,metarxn.map,by="metarxn_id")
    rownames(wilcox.results)<-wilcox.results$rxn
    W<-wilcox.results
  }else{
    wilcox.results$rxn<-rownames(wilcox.results)
    for(i in 1:nrow(wilcox.results)){
      r=wilcox.results$rxn[i]
      if(r %in% reaction.meta$reaction_no_direction){
        wilcox.results$metadata_r_id[i]<-r
      }else if(substr(r,1,nchar(r)-4) %in% reaction.meta$reaction_no_direction){
        wilcox.results$metadata_r_id[i]<-substr(r,1,nchar(r)-4) 
      }
    }
    W=left_join(wilcox.results,reaction.meta,by=c("metadata_r_id"="reaction_no_direction"))
    rownames(W)<-W$rxn
  }
  #filtering reactions
  W<-W[(W$confidence %in% c(0,4)) , ]
  W<-W[!is.na(W$EC_number),]
  W$subsystem<-as.character(W$subsystem)
  W$subsystem[grepl("Citric acid cycle",W$subsystem) & !grepl("\\[m\\]",W$formula,perl=T)]<-"Other"
  assign(paste0(a,"W.",j),W)
}

#volcano plots
select.systems<-c("Glycolysis/gluconeogenesis","Citric acid cycle","Fatty acid oxidation","Nucleotide interconversion",
                  "Pyruvate metabolism","Purine synthesis","Pentose phosphate pathway","Pyrmidine catabolism")
for(j in ls(pattern="W.")){
  data<-get(j)
  data$logP = -log10(data$adjusted_pval)
  for(i in 1:nrow(data)){
    if(data$subsystem[i] %in% select.systems){
      data$subsystem.select[i] <-data$subsystem[i] 
    }else if (data$subsystem[i] %in% amino_acid_metab){
      data$subsystem.select[i] <-"Amino Acid Metabolism"
    }else{
      data$subsystem.select[i] <-"Other"
    }
  }
  data<-data[order(data$adjusted_pval),] #sort by significance
  for(i in 1:nrow(data)){
    if(data$rxn[i] %in% rxn.names$rxn & !is.na(data$adjusted_pval[i]) & 
       data$adjusted_pval[i]<0.01 & abs(data$cohens_d[i])>0.04){ #select significant values to label
      # data$label[i]<-data$rxn[i]
      if(!(rxn.names$label[rxn.names$rxn==data$rxn[i]] %in% data$label)){ #only label each rxn once
        data$label[i]<-rxn.names$label[rxn.names$rxn==data$rxn[i]]
      }else{
        data$label[i]<-""
      }
    }else{
      data$label[i]<-""
    }
  }
  # #set labels for publication
  if(j=="metarxnW.B"){data$label[!(data$label %in% c("mevalonate kinase","fucokinase","LCFA beta oxidation"))]<-""}
  if(j=="metarxnW.CD4"){data$label[!(data$label %in% c("guanosine aminohydrolase","ketohexokinase","phosphoglucomutase","lipase","carnitine O-palmitoyltransferase","diphosphoglyceromutase"))]<-""}
  if(j=="metarxnW.CD8"){data$label[!(data$label %in% c("nucleoside-diphosphate kinase","creatine kinase","adenylate kinase","transaldolase","phosphoglucomutase","nucleoside-triphosphatase"))]<-""}
  if(j=="metarxnW.Mac"){data$label[!(data$label %in% c("cytochrome P450","glucose-6-phosphate isomerase","triose-phosphate isomerase","aldose reductase","creatine kinase","purine-nucleoside phosphorylase"))]<-""}
  ggplot(data,aes(x=cohens_d,y=logP))+geom_point(size=0.2,color="grey20")+
    geom_point(data=data[data$subsystem.select!="Other",],aes(x=cohens_d,y=logP,color=subsystem.select),size=0.3)+
    geom_hline(yintercept=1,linetype="dashed",size=0.4)+geom_vline(xintercept=0,linetype="dashed",size=0.4)+
    xlim(-1,1)+labs(title=gsub("W.","",j),x="Cohen's d (SLE.yaa - B6)",y=expression(-log[10]("Wilcoxan-adjusted P")))+theme_bw()+
    theme(legend.title=element_blank())+
    guides(color=guide_legend(override.aes=list(size=1)))+
    scale_color_brewer(palette="Set1")+
    geom_text_repel(size=2,aes(label=label),box.padding = 0.5,segment.size=0.2,max.overlaps=Inf,force=50)#, ,box.padding = 0.005,
  ggsave2(paste0(gsub("W.","",j),".volcano.png"),width=6, height=4,device="png")
  write.csv(data, file=paste0(gsub("W.","",j),".csv"))
}

#Dot plot
for(j in ls(pattern="W.")){
  data<-get(j)
  data<-data[!data$subsystem %in% c("Miscellaneous","Unassigned","Other"),]
  data<-data[!grepl("Transport|Exchange",data$subsystem),]
  data<-data[ data$subsystem %in%  names(table(data$subsystem))[table(data$subsystem) >5] , ] # only keep subsystems with at least 5 cells
  data$sig.level<-"NS"
  data$sig.level[data$adjusted_pval<0.01]<-"BH-adjusted p<0.01"
  data$up.down<-"down"
  data$up.down[data$cohens_d>=0]<-"up"
  ggplot(data = data, aes(x = cohens_d, y = reorder(subsystem,cohens_d,FUN=median))) + 
    geom_point(mapping = aes_string(color = "up.down",alpha="sig.level")) +
    geom_vline(xintercept=0,linetype="dotted")+
    scale_alpha_manual(values=c(NS=0.2,`BH-adjusted p<0.01`=1))+
    scale_color_manual(values=c(up="red",down="black"))+
    theme_bw()+
    theme(legend.position="bottom")+
    labs(x="Cohen's d",y="",title=gsub("W.","",j))+
    guides(color=F,alpha = guide_legend(title=""))
  ggsave2(paste0(gsub("W.","",j),".dotplot.png"),width=6, height=8,device="png")
}


#Filter core reactions
core.rxns<-reaction.meta
core.rxns<-core.rxns[(core.rxns$confidence %in% c(0,4)) , ]
core.rxns<-core.rxns[!is.na(core.rxns$EC_number),]
if(a=="metarxn"){
  core.metarxns<-unique(metarxn.map$metarxn_id[metarxn.map$rxn %in% paste0(core.rxns$reaction_no_direction,"_pos") |
                                                 metarxn.map$rxn %in% paste0(core.rxns$reaction_no_direction,"_neg")])
  core.reaction.consistencies<-subset(reaction.consistencies,rownames(reaction.consistencies) %in% core.metarxns)
}else{
  core.reaction.consistencies<-subset(reaction.consistencies, rownames(reaction.consistencies) %in% paste0(core.rxns$reaction_no_direction,"_pos") |
                                        rownames(reaction.consistencies) %in% paste0(core.rxns$reaction_no_direction,"_neg") )
}

#collapsed metarxn.map
metarxn.map.collapse<-left_join(metarxn.map,rxn.names,by="rxn")
metarxn.map.collapse<- as.data.table(metarxn.map.collapse)[, lapply(.SD, data_concater), by=metarxn_id]
metarxn.map.collapse2<-left_join(metarxn.map,rxn.names,by="rxn")
metarxn.map.collapse2<- as.data.table(metarxn.map.collapse2)[, lapply(.SD, data_concater2), by=metarxn_id]
select.systems<-c("Glycolysis/gluconeogenesis","Citric acid cycle","Fatty acid oxidation","Nucleotide interconversion",
                  "Pyruvate metabolism","Purine synthesis","Pentose phosphate pathway","Pyrmidine catabolism")
for(i in 1:nrow(metarxn.map.collapse2)){
  if(metarxn.map.collapse2$subsystem[i] %in% select.systems){
    metarxn.map.collapse2$subsystem.select[i] <-metarxn.map.collapse2$subsystem[i] 
  }else if (metarxn.map.collapse2$subsystem[i] %in% amino_acid_metab){
    metarxn.map.collapse2$subsystem.select[i] <-"Amino Acid Metabolism"
  }else{
    metarxn.map.collapse2$subsystem.select[i] <-"Other"
  }
}
  
  
###run PCA
df2<-t(core.reaction.consistencies)
res.pca <- PCA(df2,  graph = FALSE)
pca.plot<-as.data.frame(res.pca$ind$coord)[1:5]
pca.plot<-cbind(pca.plot,cell.meta)
ggplot(data=pca.plot,aes(x=Dim.1,y=Dim.3,color=condition))+geom_point(size=0.5)
var.plot<-facto_summarize(res.pca,element="var")
var.plot$metarxn_id<-as.integer(rownames(var.plot))
pca.plot$cell_type<-factor(pca.plot$cell_type,levels=names(sort(table(pca.plot$cell_type),decreasing=T)))
var.plot.select<-head(var.plot[order(-var.plot$contrib),],n=10)
scale=40
ggplot(data=pca.plot,aes(x=Dim.1,y=Dim.2,color=cell_type))+
  geom_point(size=0.5)+
  geom_segment(data=var.plot.select,aes(x=0,y=0,xend=(Dim.1*scale),yend=(Dim.2*scale)),
               arrow = arrow(length = unit(0.4, "picas"),type="closed"),color="black")+
  annotate("text", x = (var.plot.select$Dim.1*scale*1.1), y = (var.plot.select$Dim.2*scale*1.1),label = var.plot.select$metarxn_id)
#manual labelling
x=1247
metarxn.map[metarxn.map$metarxn_id==x,]
manual.labels=data.frame(metarxn_id=c(1414, 1258, 1447, 1433, 1247),
                         label2=c("oxalate exchange", "FA oxidation","ER transport","carnitine shuttle","cholesterol synthesis"))
var.plot<-left_join(var.plot,manual.labels,by="metarxn_id")
var.plot.select<-var.plot[var.plot$metarxn_id %in% manual.labels$metarxn_id,]
scale=40
ggplot(data=pca.plot,aes(x=Dim.1,y=Dim.2,color=cell_type))+
  geom_point(size=0.5)+
  geom_segment(data=var.plot.select,aes(x=0,y=0,xend=(Dim.1*scale),yend=(Dim.2*scale)),
               arrow = arrow(length = unit(0.4, "picas"),type="closed"),color="black")+
  annotate("text", x = (var.plot.select$Dim.1*scale*1.3), y = (var.plot.select$Dim.2*scale*1.3),label = var.plot.select$label2)+
  theme_classic()+labs(x="PC1",y="PC2")+theme(legend.title=element_blank())#+ 
ggsave2("PCA.cluster.png",width=6, height=5,device="png")
#PCA plot with histograms
scat<-ggplot(data=pca.plot,aes(x=Dim.1,y=Dim.2,color=condition))+
  geom_point(size=0.5)+scale_color_manual(values=c("black","red"))+
  labs(x="PC1",y="PC2")+theme_classic() +
  theme(legend.position="none")
xden<-ggplot(pca.plot, aes(Dim.1, fill=condition)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c("black","red","grey")) + theme_classic()+
  theme(legend.position = "none",axis.title = element_blank())
yden<-ggplot(pca.plot, aes(Dim.2, fill=condition)) + 
  geom_density(alpha=.5) + coord_flip()+
  scale_fill_manual(values = c("black","red","grey")) + theme_classic()+
  theme(legend.position = "none",axis.title = element_blank())
blankPlot <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
grid.arrange(xden, blankPlot, scat, yden, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
ggsave2("PCA.condition.png",plot=arrangeGrob(xden, blankPlot, scat, yden, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4)),
        width=8, height=8,device="png")

###run UMAP
custom.config = umap.defaults
umap<-umap(df2, config=custom.config)
umap.plot<-cbind(as.data.frame(umap$layout),cell.meta)
umap.plot$cell_type<-factor(umap.plot$cell_type,levels=names(sort(table(umap.plot$cell_type),decreasing=T)))
ggplot(data=umap.plot,aes(x=V1,y=V2,color=cell_type))+
  geom_point(size=0.5)+facet_grid(. ~ condition)+
  theme(strip.background = element_rect(colour="black", fill="white",linetype="solid"),
        panel.border = element_rect(colour = "black", size=1,fill=NA),panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.key = element_rect(fill = "white", colour = NA),legend.title=element_blank())+
  labs(x="UMAP_1",y="UMAP_2")
ggsave2("compass.umap.png",width=12, height=6,device="png")

load("../graphed.RData")


#correlate PCA with gene expression
genes<-as.data.frame(t(as.matrix(SLE.obj.combined[["RNA"]]@data)))
mat<-pca.plot
select.genes<-read_csv("select.genes.csv")
cor.mat<-setNames(data.frame(matrix(NA,ncol=5,nrow=nrow(select.genes)),row.names=select.genes$gene),colnames(mat)[1:5])
for(i in select.genes$gene){
  for(j in 1:5){
    cor<-cor.test(mat[[j]],genes[[i]],method="spearman")
    if(cor$p.value<0.05){
      cor.mat[i,j]<-cor$estimate
    }else{
      cor.mat[i,j]<-NA
    }
  }
}
cor.mat<-cor.mat[rowSums(is.na(cor.mat)) != ncol(cor.mat), ] #get rid of empty rows
pca.genes.cor.mat<-cor.mat

#correlate metarxns with gene expression
mat<-setNames(data.frame(t(core.reaction.consistencies)),rownames(core.reaction.consistencies))
cor.mat<-setNames(data.frame(matrix(NA,ncol=ncol(mat),nrow=nrow(select.genes)),row.names=select.genes$gene),colnames(mat))
for(i in select.genes$gene){
  for(j in 1:ncol(mat)){
    cor<-cor.test(mat[[j]],genes[[i]],method="spearman")
    if(cor$p.value<0.05){
      cor.mat[i,j]<-cor$estimate
    }else{
      cor.mat[i,j]<-NA
    }
  }
}
cor.mat<-cor.mat[rowSums(is.na(cor.mat)) != ncol(cor.mat), ] #get rid of empty rows
metarxns.genes.cor.mat<-cor.mat

# compute module scores
m_df = msigdbr(species = "Mus musculus")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
gs.id<-m_df %>% split(x = .$gene_symbol, f = .$gs_id) 
modules<-read.csv("gene.sets.csv",header=T)
for(h in 110:nrow(modules)){#
  i<-modules$name[h]
  if(modules$type[h]=="kegg"){genes<-bitr(gsub("kegg.","",i), fromType="PATH", toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL} #see: keytypes(org.Mm.eg.db)
  if(modules$type[h]=="GO"){genes<-bitr(i, fromType="GO", toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL} # see: keytypes(org.Mm.eg.db)
  if(modules$type[h]=="gsea"){genes<-gs.id[[i]]} #egmt, egmt2, msgidb
  mod.genes<-genes[genes %in% rownames(SLE.obj.combined)]
  if(length(mod.genes)>5){
    SLE.obj.combined <- AddModuleScore(object = SLE.obj.combined,features = list(mod.genes),name = paste0("mod.",i,modules$info[h],".score"))
    #UMAP
    mod.name<-gsub(":|-",".",paste0("mod.",i,modules$info[h],".score1"))
    FeaturePlot(SLE.obj.combined, features =mod.name, split.by = "condition",
                cols=viridis(100, begin = 0),order=T,pt.size=0.001,combine=T)
    ggsave2(paste0(i,".",modules$info[h],".umap.png"),width=6.5, height=3,device="png")
    #UMAP compare
    p<-FeaturePlot(SLE.obj.combined, features = mod.name, split.by = "condition",cols=viridis(100, begin = 0),
                   order=T,pt.size=0.001,combine=F)
    for(j in 1:length(p)) {
      p[[j]] <- p[[j]] + NoLegend()+NoAxes()+
        theme(panel.border = element_rect(colour = "black", size=1),
              plot.title=element_blank(),axis.title.y.right=element_blank(),
              axis.line=element_blank())
    }
    cowplot::plot_grid(p[[1]],p[[2]],ncol=2)
    ggsave2(paste0(i,".",modules$info[h],".umap2.png"),width=6, height=2.75,device="png")
    #VLN plot
    VlnPlot(SLE.obj.combined, features = mod.name, split.by = "condition",
            group.by = "my.clusters2",cols=c("red","grey"),pt.size = 0,combine=F)
    ggsave2(paste0(i,".",modules$info[h],".vln.png"),width=5, height=4,device="png")
  }
}
colnames(SLE.obj.combined@meta.data)

#correlation mod with PCA
mod.scores<-select(SLE.obj.combined@meta.data,contains("mod."))
mat<-pca.plot
cor.mat<-setNames(data.frame(matrix(NA,ncol=5,nrow=ncol(mod.scores)),row.names=colnames(mod.scores)),colnames(mat)[1:5])
for(i in colnames(mod.scores)){
  for(j in 1:5){
    cor<-cor.test(mat[[j]],mod.scores[[i]],method="spearman")
    if(cor$p.value<0.05){
      cor.mat[i,j]<-cor$estimate
    }else{
      cor.mat[i,j]<-NA
    }
  }
}
cor.mat<-cor.mat[rowSums(is.na(cor.mat)) != ncol(cor.mat), ] #get rid of empty rows
pca.mods.cor.mat<-cor.mat

#mod x mod scatter
#metarxn x mods heatmap
scatter<-as.data.frame(t(metarxns.mods.cor.mat.all))
modules$mod.id<-gsub(":|-",".",paste0("mod.",modules$name,modules$info,".score1"))
#add metarxn metadata
scatter$metarxn_id<-rownames(scatter)
scatter<-merge(scatter,metarxn.map.collapse2,by="metarxn_id")
x.mod=c("KEGG_BCR_signaling","KEGG_TCR_signaling","KEGG_antigen_presentation","KEGG_complement_cascade")
y.mod="KEGG_oxidative_phosphorylation"
plot.list<-list()
for(i in 1:length(x.mod)){
  plot.list[[i]]<-ggplot(data=scatter,aes_string(x=modules$mod.id[modules$label==x.mod[i]],y=modules$mod.id[modules$label==y.mod],color="subsystem.select"))+
    geom_point(size=0.2)+theme_bw()+labs(x=x.mod[i],y=y.mod)+theme(legend.title=element_blank())+
    guides(color=guide_legend(override.aes=list(size=1)))+
    geom_vline(xintercept=0,linetype="dotted")+
    geom_hline(yintercept=0,linetype="dotted")+
    scale_color_brewer(palette="Paired")
}
ggarrange(plotlist=plot.list, common.legend=TRUE,legend="right")
ggsave2("mods.mods.metarxn.scatter.png",width=14, height=7,device="png")
ggarrange(plotlist=plot.list, common.legend=TRUE,legend="right")
ggsave2("mods.mods.metarxn.scatter2.png",width=8, height=6,device="png")

#pca x mods heatmap
mat<-pca.mods.cor.mat
colnames(mat)<-c("PC1","PC2","PC3","PC4","PC5")
pheatmap(mat)
#relabel rows
modules$mod.id<-gsub(":|-",".",paste0("mod.",modules$name,modules$info,".score1"))
annot.row<-data.frame(mod.id=rownames(mat),order=1:nrow(mat))
annot.row<-merge(annot.row,unique(modules[,c("mod.id","label")]),by="mod.id")
annot.row<-annot.row[order(annot.row$order),]
rownames(mat)<-annot.row$label
#plot
range <- max(abs(mat),na.rm=T)
p1<-pheatmap(mat,fontsize_col=7,fontsize_row=5,cluster_cols = F,angle_col=0,
             color = colorRampPalette(c("blue", "white", "red"))(100),breaks = seq(-range, range, length.out = 100))
save_pheatmap_png <- function(x, filename, width=1500, height=1500, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "pca.mods.heatmap.png")