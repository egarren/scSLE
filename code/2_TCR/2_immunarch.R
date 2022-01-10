rm(list=ls())
# Load the package
library(immunarch)
library(cowplot)
library(ggrepel)
library(Seurat)
library(dplyr)
library(ggpubr)
library(EnhancedVolcano)
library(ggseqlogo)
library(ggplot2)
library(ggsci)
library(viridis)
library(RColorBrewer)
library(plyr)
library(patchwork)
library(gridExtra)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}

load("clone.ab.data.RData")

##VDJtools
#aa strength scores
meta<-read.table("./vdjtools.cluster/clust.out/annot/metadata.txt", header=T, fill=T, stringsAsFactors=F, sep="\t")
meta$aa.score<-0
for (i in 1:nrow(meta)) {
  df<-read.table(paste0("./vdjtools.cluster/clust.out/annot/",meta$file_name[i]), header=T, fill=T, stringsAsFactors=F, sep="\t")
  meta$aa.score[i]<-sum(df$freq*df$aaprop.strength,na.rm=T)
}
meta$Cluster<-as.factor(meta$Cluster) 
ggbarplot(meta, x = "Cluster", y = "aa.score",add = c("mean_se", "jitter"),color = "condition",
          position = position_dodge(0.8), legend="right",ylim=c(4,5.7),palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", method="t.test",size=2.5,label.y=5.5)+
  # stat_compare_means(method="anova",label.y=6)+
  labs( y = "Weighted strong aa score")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = c("CD4","CD8","CD8-eff","Treg","Cd74","CD4-eff"))
ggsave("aa.score.cluster.png",width=5.5, height=3.5,device="png")

##Immunarch
#by condition
#create repertoire folder
dir.create("./rep")
# dir.create("./immunarch.results")
rep.tabs<-list.files("./vdjtools/out/annot",".txt$")
file.copy(file.path("./vdjtools/out/annot",rep.tabs),"./rep")
df<-read.table("./vdjtools/out/annot/metadata.txt",header=T)
df<-cbind(Sample=gsub(".txt$","",df$file_name),df)
write.table(df, file = "./rep/metadata.txt", sep = "\t",quote=F,row.names = FALSE)

# Load the data to the package
# setwd("./immunarch.results")
immdata = repLoad("./rep")
save(immdata, file="immdata.Rdata")
#Basic anlaysis (explore)
exp_vol = repExplore(immdata$data, .method = "volume")
vis(exp_vol, .by = c("condition"), .meta = immdata$meta) #.by=c("condition',"cluster")
ggsave2("number.unique.clonotypes.png",width=4, height=6,device="png") #for publication use fixVis
exp_len = repExplore(immdata$data, .method = "len", .col = "aa")
vis(exp_len,.by=c("condition"), .meta = immdata$meta)
ggsave2("cdr3.length.png",width=20, height=6,device="png") #for publication use fixVis
exp_cnt = repExplore(immdata$data, .method = "count")
vis(exp_cnt)#.by=c("condition)"), .meta = immdata$meta)
ggsave2("clonotype.abundance.png",width=6, height=6,device="png") #for publication use fixVis
#clonality
for(i in c("clonal.prop","homeo","top","rare")){
  imm_pr = repClonality(immdata$data, .method = i)
  png(paste0(i,".clonality.png"),width=10,height=6,units="in",res=200)
  grid.arrange(vis(imm_pr),vis(imm_pr,.by=c("condition"),.meta=immdata$meta),ncol=2)
  dev.off()
}

#public repertoire (identify shared clones)
pr = pubRep(immdata$data, "aa", .coding = T, .verbose = F)
# vis(pr,.by = c("condition"), .meta = immdata$meta)
# ggsave2("public.clonotypes.png",width=7, height=7,device="png")
pr.B6 = pubRepFilter(pr, immdata$meta, c(condition = "B6"))
pr.SLE.yaa = pubRepFilter(pr, immdata$meta, c(condition = "SLE.yaa"))
pr.B6[is.na(pr.B6)]<-0
pr.SLE.yaa[is.na(pr.SLE.yaa)]<-0
pr.B6[["avgfreq.B6"]] = rowMeans(public_matrix(pr.B6), na.rm = T)
pr.SLE.yaa[["avgfreq.SLE.yaa"]] = rowMeans(public_matrix(pr.SLE.yaa), na.rm = T)
pr.B6[["sum.B6"]] = rowSums(public_matrix(pr.B6)[,1:2], na.rm = T)#+1
pr.SLE.yaa[["sum.SLE.yaa"]] = rowSums(public_matrix(pr.SLE.yaa)[,1:2], na.rm = T)#+1
pr.res.cdr3 = dplyr::full_join(pr.B6, pr.SLE.yaa, by = "CDR3.aa")
pr.res.cdr3.graph<-as.data.table(pr.res.cdr3)
pr.res.cdr3.graph[is.na(pr.res.cdr3.graph)]<-0
# pr.res.cdr3.graph$sum.B6[pr.res.cdr3.graph$sum.B6==0]<-1
# pr.res.cdr3.graph$sum.SLE.yaa[pr.res.cdr3.graph$sum.SLE.yaa==0]<-1
pr.res.cdr3.graph[["Samples.sum"]] = pr.res.cdr3.graph[["Samples.x"]] + pr.res.cdr3.graph[["Samples.y"]]
pr.res.cdr3.graph[["freq.ratio"]] = apply(pr.res.cdr3.graph[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.cdr3.graph[["log2FC"]]= apply(pr.res.cdr3.graph[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log2((x[2])/(x[1])))
rownames(pr.res.cdr3.graph)<-pr.res.cdr3.graph$CDR3.aa
pr.res.cdr3.graph$aa_length<-nchar(as.character(pr.res.cdr3.graph$CDR3.aa))
pr.res.cdr3.graph<-left_join(pr.res.cdr3.graph,clone.meta[,c("cdr3","chain")],by=c("CDR3.aa"="cdr3"),keep=F)
pr.res.cdr3.graph<-as.data.frame(pr.res.cdr3.graph)
plot.list<-list()
plot.list2<-list()
for(h in c("TRA","TRB")){
  pr.df<-pr.res.cdr3.graph[pr.res.cdr3.graph$chain==h&!is.na(pr.res.cdr3.graph$chain),]
  cdr3.public.FC<-pr.df[,c("CDR3.aa","Samples.sum","log2FC")]
  colnames(cdr3.public.FC)<-c("CDR3.aa","Samples.sum.mice","log2FC.mice")
  most.shared.cdr3<-pr.df[abs(pr.df$log2FC)>6&
                                        pr.df$sum.SLE.yaa>20&pr.df$sum.B6>20,]
  ggplot(pr.df,aes(x = sum.B6, y = sum.SLE.yaa,size=Samples.sum,label=CDR3.aa))+geom_point(alpha=0.25)+ 
    scale_x_continuous(trans='log2') +scale_y_continuous(trans='log2')+
    geom_abline(intercept = 0, slope = 1) +
    geom_text_repel(data = most.shared.cdr3,aes(label = CDR3.aa),size = 3,
                    color="red",min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
  ggsave2(paste0(h,".overlap.B6.SLE.yaa.cdr3.scatter.png"),width=7, height=6,device="png")
  #Volcano
  # pr.res.cdr3.volcano = dplyr::inner_join(pr.B6, pr.SLE.yaa, by = "CDR3.aa")
  # pr.res.cdr3.volcano[["log2FC"]]= apply(pr.res.cdr3.volcano[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log2((x[2])/(x[1])))
  pr.df$p_val<-NA
  for(i in 1:length(pr.df$p_val)){
    pr.df$p_val[i]<-my.ttest(pr.df[i,3:7],pr.df[i,11:15])
  }
  pr.df$p_val_adj<-p.adjust(pr.df$p_val,method="BH")
  EnhancedVolcano(pr.df,lab = pr.df$CDR3.aa,
                  x = 'log2FC',y = 'p_val',title="Public CDR3",col=c("black","black","black","red3"),
                  pCutoff = 0.05,FCcutoff = 2,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                  legendPosition="none",drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                  subtitle="", caption="",border="full",cutoffLineWidth=0,
                  gridlines.major=F,gridlines.minor=F,titleLabSize=10
  )
  ggsave2(paste0(h,".overlap.B6.SLE.yaa.cdr3.volcano.png"),width=5, height=5,device="png")
  #WebLogo on public repertoire
  B6.xp<-pr.df[pr.df$log2FC< -5,]#&pr.df$sum.B6>0.1,]
  SLE.yaa.xp<-pr.df[pr.df$log2FC>5,]#&pr.df$sum.SLE.yaa>0.1,]
  df.list<-list()
  for(j in c("B6","SLE.yaa")){
    df<-get(paste0(j,".xp"))
    if(j=="B6"){name="WT"}else{name="SLE.yaa"}
    for(k in 14){ #10:18
      df2<-df[df$aa_length==k,]
      plot.list[[paste0(h,j,k)]]<-ggplot()+geom_logo(as.character(df2$CDR3.aa))+
        theme_logo()+ggtitle(paste0(h,"_",name))+
        theme(axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5)) #,"_",k
    }
    df$condition<-name
    df.list[[j]]<-df
  }
  #cdr3 length
  aa.df<-rbindlist(df.list)
  aa.df$condition <- factor(aa.df$condition, levels = c("WT","SLE.yaa"))
  mu <- ddply(aa.df, "condition", summarise, grp.mean=mean(aa_length))
  plot.list2[[h]]<-ggplot(aa.df, aes(x=aa_length,color=condition,fill=condition)) + 
    geom_histogram(aes(y=..density..),binwidth=1,position="identity",alpha=0.2) +#geom_density(fill=NA)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+ggtitle(h)+theme_classic()+
    labs(x = "CDR3 Length", y = "") +theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())+
    scale_color_manual(values = c("black", "red"))+
    scale_fill_manual(values = c("black", "red"))
}
CombinePlots(plots=list(plot.list[[1]],plot.list[[3]],plot.list[[2]],plot.list[[4]]),ncol=2,legend="bottom")
ggsave2("cdr3.pr.weblogo.bycondition.png",width=6,height=4,device="png")
CombinePlots(plot.list2,ncol=2,legend="right")
ggsave2("cdr3.pr.aa_length.bycondition.histo.png",width=7,height=4,device="png")

#Annotating using antigen prediction databases
load("masterdb3.RData")
pr.df<-pr.res.cdr3.graph[pr.res.cdr3.graph$chain=="TRB"&!is.na(pr.res.cdr3.graph$chain),]
df.db<-unique(masterdb[masterdb$antigen!="unknown"&!is.na(masterdb$antigen)&masterdb$disease2!="species",])
# df.db<-unique(masterdb[!is.na(masterdb$chain)&masterdb$chain=="TRB"&masterdb$TRBV!=""&masterdb$Species=="Mouse"&masterdb$Tcell=="CD4",]) #
pr.res.db<-unique(left_join(x = pr.df, y = df.db[,c("cdr3.b","antigen","epitope","study2","Species","disease2","disease")], 
                        by = c("CDR3.aa"="cdr3.b"),keep=F))
pr.res.db$disease2[is.na(pr.res.db$disease2)]<-"unknown"
labels<-unique(pr.res.db[abs(pr.res.db$log2FC)>3&pr.res.db$Samples.sum>1&!is.na(pr.res.db$antigen)&
                           (pr.res.db$sum.B6>50|pr.res.db$sum.SLE.yaa>50),])
# labels2<-labels[labels$antigen %in% c("NP177","M45","IE1","Tat","EBNA4","p65"),]
ggplot(pr.res.db,aes(x = sum.SLE.yaa, y =  sum.B6,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(aes(colour=factor(disease2),fill = factor(disease2)), shape=21) + 
  scale_color_manual(values = c(brewer.pal(n = 4, name = "Dark2"),alpha("black",0.5),"#A6761D"))+
  scale_fill_manual(values = c(alpha(brewer.pal(n = 4, name = "Dark2"),0.5),"#1C00ff00",alpha("#A6761D",0.5)))+
  guides(fill = guide_legend(override.aes = list(size = 5)),color=F)+ #,size=F
  labs(x = "SLE.yaa", y = "WT", size="Samples",fill="Disease")+
  theme(legend.direction = "vertical", legend.box = "horizontal")+#theme(legend.title = element_blank())+
  geom_label_repel(data =labels,aes(label = antigen),size = 3,min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
ggsave2("annot.db.scatter.bycondition.png",width=6, height=3.5,device="png")

#pubrep using a/b pairing
df<-clone.data.ab
imm.list<-list()
for(i in unique(df$mouse_ID)){ #creating immunarch list from collapsed data
  df2<-df[df$mouse_ID==i,]
  # df2<-df2[!is.na(df2$chain),]
  # df2$clone_ID<-paste0(df2$cdr3,"_",df2$v_gene)
  df2<-transform(df2,Clones=ave(seq(nrow(df2)),cdr3,FUN=length))
  df2<-as.data.table(df2)[, lapply(.SD, data_concater), by=cdr3]
  df2$Clones<-as.numeric(as.character(df2$Clones))
  df2$Proportion<-df2$Clones/sum(df2$Clones)
  imm.list[[i]]<-tibble(
    Clones=df2$Clones, Proportion=df2$Proportion, CDR3.nt=df2$cdr3_nt,CDR3.aa=df2$cdr3,
    V.name=df2$v_gene,D.name=df2$d_gene,J.name=df2$j_gene
  )
}
meta2<-immdata$meta
meta2$Sample<-meta2$mouse_ID
pr = pubRep(imm.list, "aa", .coding = T, .verbose = F)
# vis(pr,.by = c("condition"), .meta = immdata$meta)
# ggsave2("public.clonotypes.png",width=7, height=7,device="png")
pr.B6 = pubRepFilter(pr, meta2, c(condition = "B6"))
pr.SLE.yaa = pubRepFilter(pr, meta2, c(condition = "SLE.yaa"))
pr.B6[is.na(pr.B6)]<-0
pr.SLE.yaa[is.na(pr.SLE.yaa)]<-0
pr.B6[["avgfreq.B6"]] = rowMeans(public_matrix(pr.B6), na.rm = T)
pr.SLE.yaa[["avgfreq.SLE.yaa"]] = rowMeans(public_matrix(pr.SLE.yaa), na.rm = T)
pr.B6[["sum.B6"]] = rowSums(public_matrix(pr.B6)[,1:2], na.rm = T)
pr.SLE.yaa[["sum.SLE.yaa"]] = rowSums(public_matrix(pr.SLE.yaa)[,1:2], na.rm = T)
pr.res.cdr3 = dplyr::full_join(pr.B6, pr.SLE.yaa, by = "CDR3.aa")
pr.res.cdr3.graph<-as.data.table(pr.res.cdr3)
pr.res.cdr3.graph[is.na(pr.res.cdr3.graph)]<-0
# pr.res.cdr3.graph$sum.B6[pr.res.cdr3.graph$sum.B6==0]<-1
# pr.res.cdr3.graph$sum.SLE.yaa[pr.res.cdr3.graph$sum.SLE.yaa==0]<-1
pr.res.cdr3.graph[["Samples.sum"]] = pr.res.cdr3.graph[["Samples.x"]] + pr.res.cdr3.graph[["Samples.y"]]
pr.res.cdr3.graph[["freq.ratio"]] = apply(pr.res.cdr3.graph[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log10(x[1])/log10(x[2]))
pr.res.cdr3.graph[["log2FC"]]= apply(pr.res.cdr3.graph[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log2((x[2])/(x[1])))
pr.res.cdr3.graph[["log2FC.avg"]]= apply(pr.res.cdr3.graph[, c("avgfreq.B6", "avgfreq.SLE.yaa")],1, function(x) log2((x[2])/(x[1])))
rownames(pr.res.cdr3.graph)<-pr.res.cdr3.graph$CDR3.aa
pr.df<-pr.res.cdr3.graph
cdr3.public.FC<-pr.df[,c("CDR3.aa","Samples.sum","log2FC")]
colnames(cdr3.public.FC)<-c("CDR3.aa","Samples.sum.mice","log2FC.mice")
labels<-pr.df[abs(pr.df$log2FC)>2&(pr.df$sum.B6>20|pr.df$sum.SLE.yaa>50)&pr.df$Samples.sum>1,]
ggplot(pr.df,aes(x = sum.SLE.yaa, y =  sum.B6,size=Samples.sum))+theme_classic()+
  scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10')+ scale_size(range = c(1, 10))+
  geom_point(color=alpha("black",0.5),fill="#1C00ff00", shape=21) + 
  labs(x = "SLE.yaa", y = "WT", size="Samples")+
  geom_text_repel(data =labels,aes(label = CDR3.aa),size = 3,color="forestgreen",
                   min.segment.length=unit(0,'lines'),nudge_x=0.1,nudge_y=0.1,segment.size=0.25)
ggsave2("overlap.TCRab.bycondition.png",width=5, height=4,device="png")
#Volcano
# pr.res.cdr3.volcano = dplyr::inner_join(pr.B6, pr.SLE.yaa, by = "CDR3.aa")
# pr.res.cdr3.volcano[["log2FC"]]= apply(pr.res.cdr3.volcano[, c("sum.B6", "sum.SLE.yaa")],1, function(x) log2((x[2])/(x[1])))
pr.df$p_val<-NA
for(i in 1:length(pr.df$p_val)){
  pr.df$p_val[i]<-my.ttest(pr.df[i,3:7],pr.df[i,11:15])
}
pr.df$p_val_adj<-p.adjust(pr.df$p_val,method="BH")
EnhancedVolcano(pr.df,lab = pr.df$CDR3.aa,
                x = 'log2FC.avg',y = 'p_val',title="Public CDR3",col=c("black","black","black","red3"),
                pCutoff = 0.05,FCcutoff = 1,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                legendPosition="none",drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                subtitle="", caption="",border="full",cutoffLineWidth=0,
                gridlines.major=F,gridlines.minor=F,titleLabSize=10
)
ggsave2("overlap.TCRab.bycondition.volcano.png",width=5, height=5,device="png")


# Compute and visualise gene usage with samples, grouped by their disease status
gu = geneUsage(immdata$data, .type="family", .norm=T)#, .gene="musmus.trbv") #see gene_stats()
vis(gu, .by="condition", .meta=immdata$meta)#, .plot="box")
ggsave2("vgene.usage.bar.png",width=24, height=8,device="png")
# png(paste0("vgene.usage.tree.png"),width=12,height=12,units="in",res=200)
# vis(gu, .plot="tree")
# dev.off()

# Compute Jensen-Shannon divergence among gene distributions of samples,
# cluster samples using the hierarchical clustering and visualise the results
plot.list<-list()
plot.list2<-list()
for(i in c("js","cor","cosine")){
  for(l in c("heatmap", "heatmap2", "circos")){ #, "radar"
    gu.clust<-geneUsageAnalysis(gu, .method = i)
    png(paste0("vgene.usage.overlap.",i,".",l,".png"),width=6,height=6,units="in",res=200)
    vis(gu.clust,.plot=l,.title=i)
    dev.off()
  }
  gu.clust<-geneUsageAnalysis(gu, .method = paste0(i,"+dbscan"))
  plot.list[[i]]<-vis(gu.clust)+ggtitle(paste0(i,".dbscan"))
  for(k in c("pca","mds","tsne")){
    gu.clust = geneUsageAnalysis(gu, .method = paste0(i,"+",k,"+dbscan"))
    plot.list[[paste0(i,"+",k)]]<-vis(gu.clust)+ggtitle(paste0(i,".",k,".dbscan"))
    for(j in c("kmeans","dbscan")){
      df<-as.data.frame(geneUsageAnalysis(gu, .method = paste0(i,"+",k,"+",j))$data)
      df$Sample<-rownames(df)
      df<-left_join(x = df, y = immdata$meta, by = "Sample",keep=F)
      plot.list2[[paste0(i,"+",k,"+",j)]]<-ggplot(df,aes(x = DimI, y =DimII,label=Sample,color=condition))+
        geom_text()+ggtitle(paste0(i,"+",k,"+",j))
    }
  }
  for(j in c("hclust","kmeans")){
    gu.clust<-geneUsageAnalysis(gu, .method = paste0(i,"+",j))
    png(paste0("vgene.usage.overlap.",i,".",j,".png"),width=6,height=3,units="in",res=200)
    vis(gu.clust, .title=paste0(i,"+",j))
    dev.off()
    for(k in c("pca","mds","tsne")){
      gu.clust = geneUsageAnalysis(gu, .method = paste0(i,"+",k,"+",j))
      png(paste0("vgene.usage.overlap.",i,".",j,".",k,".png"),width=6,height=3,units="in",res=200)
      vis(gu.clust,.title=paste0(i,"+",k,"+",j))
      dev.off()
    }
  }
}
CombinePlots(plot.list,ncol=4)
ggsave2("vgene.usage.overlap.dbscan.png",width=12,height=9,device="png")
CombinePlots(plot.list2,ncol=6, legend="right")
ggsave2("vgene.usage.overlap.dimred.labeled.png",width=18,height=9,device="png")

#Spectratyping
p1 = vis(spectratype(immdata$data[[1]], .quant = "id", .col = "nt"))
p2 = vis(spectratype(immdata$data[[1]], .quant = "count", .col = "aa+v"))
png("spectratype.png",width=12,height=4,units="in",res=200)
grid.arrange(p1, p2, ncol = 2)   
dev.off()

# Compare diversity of repertoires and visualise samples, grouped by two parameters
plot.list<-list()
for(j in c("nt","aa","v","aa+v")){
  for (i in c("chao1","hill","div","gini.simp","inv.simp","d50")){#"gini","raref"
    div = repDiversity(immdata$data, .method = i,.col=j)
    plot.list[[paste0(i,"+",j)]]<-vis(div, .by="condition", .meta=immdata$meta)
  }
}
CombinePlots(plot.list,ncol=6)
ggsave2("diversity.png",width=15, height=15,device="png")
plot.list<-list()
for(j in c("nt","aa","v","aa+v")){
  div = repDiversity(immdata$data, "gini",.col=j)
  temp<-data.frame(Sample=rownames(div),Value=div[,1])
  div<-add_class(temp,"immunr_ginisimp")
  plot.list[[j]]<-vis(div, .by="condition", .meta=immdata$meta)
}
CombinePlots(plot.list,ncol=1,legend="right")
ggsave2("diversity.gini.png",width=3, height=15,device="png")
plot.list<-list()
for(j in c("nt","aa","aa+v")){
  div = repDiversity(immdata$data, "raref",.col=j)
  plot.list[[j]]<-vis(div, .by="condition", .meta=immdata$meta)
  plot.list[[paste0(j,"samples")]]<-vis(div)
}
CombinePlots(plot.list,ncol=2,legend="right")
ggsave2("diversity.raref.png",width=8, height=15,device="png")


#clonotype tracking
pr = pubRep(immdata$data, "aa", .coding = T, .verbose = F)
target = head(pr$CDR3.aa,n=10)
tc = trackClonotypes(immdata$data, target, .col = "aa")
vis(tc, .order=order(immdata$meta$condition)) #.plot = "area"
ggsave2("clonotype.track.png",width=8, height=6,device="png")
#kmer stats
kmers = getKmers(immdata$data, 5)
vis(kmers,.head=20,.position="stack",.log=T)
ggsave2("kmer.png",width=6, height=6,device="png")
# for(i in immdata$meta$Sample){
#   kmers = getKmers(immdata$data[[i]], 5)
#   for(j in c("freq","prob","wei","self")){
#     kp<-kmer_profile(kmers,j) #compute sequence motif matrix
#     p1 = vis(kp)
#     p2 = vis(kp, .plot = "seq")
#     png(paste0("kmer.profile.",i,".",j,".png"),width=12,height=6,units="in",res=200)
#     grid.arrange(p1, p2, ncol=2)
#     dev.off()
#   }
# }
save.image("immunarch.bycondition.RData")


###by clusters
#create repertoire folders
dir.create("./clust.rep")
rep.tabs<-list.files("./vdjtools.cluster/clust.out/annot",".txt$")
file.copy(file.path("./vdjtools.cluster/clust.out/annot",rep.tabs),"./clust.rep")
df<-read.table("./vdjtools.cluster/clust.out/annot/metadata.txt",header=T)
df<-cbind(Sample=gsub(".txt$","",df$file_name),df)
write.table(df, file = "./clust.rep/metadata.txt", sep = "\t",quote=F,row.names = FALSE)

# Load the data to the package
# setwd("./bycluster.immunarch.results")
immdata = repLoad("./clust.rep")
save(immdata, file="immdata.clust.Rdata")
#Basic anlaysis (explore)
exp_vol = repExplore(immdata$data, .method = "volume")
vis(exp_vol, .by = c("Cluster","condition"), .meta = immdata$meta,.test=F) #.by=c("condition',"cluster")
ggsave2("bycluster.number.unique.clonotypes.png",width=6, height=6,device="png") #for publication use fixVis
exp_cnt = repExplore(immdata$data, .method = "count")
vis(exp_cnt,.by=c("Cluster","condition"), .meta = immdata$meta)
ggsave2("bycluster.clonotype.abundance.png",width=6, height=6,device="png") #for publication use fixVis
#clonality
imm_pr = repClonality(immdata$data, .method = "clonal.prop")
vis(imm_pr,.by=c("Cluster","condition"),.meta=immdata$meta,.test=F)
ggsave2("bycluster.clone.prop.clonality.png",width=6, height=6,device="png") #for publication use fixVis
#nicer graph
imm_pr2<-as.data.frame(imm_pr)
imm_pr2$Sample<-rownames(imm_pr2)
imm_pr2<-left_join(imm_pr2,immdata$meta,by="Sample")
ggbarplot(imm_pr2, x = "Cluster", y = "Clones", add = c("mean_se","jitter"),color = "condition",
          position = position_dodge(0.8),legend="right")+
  stat_compare_means(aes(group = condition),label = "p.format")+
  labs(x="Cluster", y = "Number of Clones in top 10%")
ggsave2("bycluster.clone.prop2.png",width=10, height=7,device="png")

#using a/b pairs:
#create new imm db
df<-clone.data.ab.seurat
imm.list<-list()
for(i in unique(df$mouse_ID)){ #creating immunarch list from collapsed data
  for(j in unique(df$my.clusters)){
    df2<-df[df$mouse_ID==i&df$my.clusters==j,]
    # df2<-df2[!is.na(df2$chain),]
    # df2$clone_ID<-paste0(df2$cdr3,"_",df2$v_gene)
    df2<-transform(df2,Clones=ave(seq(nrow(df2)),cdr3,FUN=length))
    df2<-as.data.table(df2)[, lapply(.SD, data_concater), by=cdr3]
    df2$Clones<-as.numeric(as.character(df2$Clones))
    df2$Proportion<-df2$Clones/sum(df2$Clones)
    imm.list[[paste0(i,"_.clust",j)]]<-tibble(
      Clones=df2$Clones, Proportion=df2$Proportion, CDR3.nt=df2$cdr3_nt,CDR3.aa=df2$cdr3,
      V.name=df2$v_gene,D.name=df2$d_gene,J.name=df2$j_gene
    )
  }
}

# Compare diversity of repertoires and visualise samples, grouped by two parameters
plot.list<-list()
for(j in c("aa","v","aa+v")){#"nt",
  for (i in c("chao1","hill","div","gini.simp","inv.simp","d50")){#"gini","raref"
    div = repDiversity(imm.list, .method = i,.col=j)
    plot.list[[paste0(i,"+",j)]]<-vis(div, .by=c("Cluster","condition"), .meta=immdata$meta,.test=F)
  }
}
CombinePlots(plot.list,ncol=6)
ggsave2("bycluster.diversity.png",width=50, height=15,device="png",limitsize=F)
#nicer graph
div = as.data.frame(repDiversity(imm.list, .method = "chao1",.col="aa"))
div$Sample=rownames(div)
div2<-left_join(div,immdata$meta,by="Sample")
ggbarplot(div2, x = "Cluster", y = "Estimator", add = c("mean_se","jitter"),color = "condition",
          position = position_dodge(0.8),legend="right")+
  stat_compare_means(aes(group = condition),label = "p.format")+
  labs(x="Cluster", y = "Diversity")
# div2$condition2<-"WT"
# div2$condition2[div2$condition=="SLE.yaa"]<-"SLE.yaa"
# div2$condition2 <- factor(div2$condition2, levels = c("WT","SLE.yaa"))
ggbarplot(div2, x = "Cluster", y = "Estimator", add = c("mean_se","jitter"),color = "condition",
          position = position_dodge(0.8),legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.format", size=2)+#,label.y=200
  labs(y = "Chao1 Diversity")+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  # scale_color_manual(labels=c("WT","SLE.yaa"))+
  scale_x_discrete(labels=c("CD4","CD8","CD8-eff","Treg","Cd74","CD4-eff")) #+scale_y_log10()
ggsave2("bycluster.chao1.png",width=5, height=4,device="png")


save.image("immunarch.bycluster.RData")






