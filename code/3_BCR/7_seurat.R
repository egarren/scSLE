rm(list=ls())
library(Seurat)
library(ggpubr)
gm_mean = function(x, na.rm=TRUE){  #geometric mean functions for graphing
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 
barcoder2<-function(df, prefix){
  df$ms.barcoder<-prefix
  df$orig.barcode<-df$cell_id
  # df$barcode <- gsub("-1", "", df$cell_id)
  df$barcode <- paste0(prefix, df$cell_id)
  df$ms_clone <- paste0(prefix, df$clone_id) 
  df$ms_v <- paste0(prefix, df$v_gene) 
  df$ms_j<- paste0(prefix, df$j_gene) 
  df$ms_cdr3 <- paste0(prefix, df$cdr3) 
  df
}


#load cdr3 data
metadata<-read.csv("SLE.metadata.csv", header=T) #define path
metadata[] <- lapply(metadata, as.character)
mice<-metadata$mouse_ID
for (i in mice){
  df<-read.table(file=paste0("./cellranger/",i,"/filtered_contig_heavy_germ-pass.tsv"),header=T,sep="\t",stringsAsFactors = F)
  assign(paste0(i,".b.imm"),barcoder2(df,prefix=paste0(i,"_")))
}
ExampleDb<-rbind(m232.b.imm,m233.b.imm,m234.b.imm,m235.b.imm)
ExampleDb<-left_join(ExampleDb,metadata,by="ms.barcoder")

#add to ExampleDb
store<-ExampleDb
for(i in c("db_props","clones","db_obs","db_obs_v","db_obs_ch","expected")){
  df<-get(i)
  keep<-setdiff(colnames(df),colnames(store))
  new_df<-df[,keep]
  colnames(new_df)<-paste0(i,".",colnames(new_df))
  key<-"barcode"
  if(i=="clones"){key="ms_clone"}
  new_df[[key]]<-df[[key]]
  ExampleDb<-left_join(ExampleDb,new_df,by=key)
}
rownames(ExampleDb)<-ExampleDb$barcode

##clone size
# all clones w/n all mice, comparing conditions
df.list<-list()
for(j in unique(ExampleDb$condition)){
  df<-ExampleDb[ExampleDb$condition ==j,]
  freq.tab<-as.data.frame(table(df$clone_id))  #ms_cdr3 frequency table for given cluster
  freq.tab$condition<-j
  df.list[[j]]<-freq.tab #add to list for each cluster
}
df<-rbindlist(df.list)
ggbarplot(df, x = "condition", y = "Freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="none",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),method="wilcox.test")+
  labs(x="", y = "Clone Size") #+scale_y_log10()
ggsave2("bar.dot.condition.clone.freq.png",width=3, height=4,device="png")
mu <- ddply(df, "condition", summarise, grp.mean=gm_mean(Freq))
# df$condition <- factor(df$condition, levels = c("WT","SLE.yaa"))
ggplot(df, aes(x=Freq,color=condition,fill=condition)) + 
  geom_histogram(aes(y=..density..),bins=10,position="identity",alpha=0.2) +#geom_density(fill=NA)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=condition),linetype="dashed")+theme_classic()+
  labs(x = "Clone Size", y = "Density")+scale_x_log10()+theme(legend.title = element_blank())+
  scale_color_manual(values = c("black", "red"))+
  scale_fill_manual(values = c("black", "red"))
ggsave2("histogram.condition.clone.freq.png",width=4, height=3,device="png")
#mut freq
ggbarplot(ExampleDb, x = "condition", y = "db_obs.mu_freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="none",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),method="wilcox.test")+
  labs(x="", y = "Mutation frequency") #+scale_y_log10()
ggsave2("bar.dot.condition.mu.freq.png",width=3.5, height=4,device="png")

#GEO export
df<-ExampleDb
write.csv(df,"scBCRseq_barcodes.csv")

#add to Seurat
load("graphed.RData")
DATA<-SLE.obj.combined
DimPlot(DATA, reduction = "umap", label=T)
keep<-setdiff(colnames(ExampleDb),colnames(DATA@meta.data))
new_df<-ExampleDb[,keep]
DATA<-AddMetaData(DATA, metadata=new_df)
save.image("immcantation.seurat.RData")



#all clones w/n all mice w/n cluster
df.list <- list()
plot.list<-list()
ExampleDb.seurat<-xDATA@meta.data
mycluster.names<-names(sort(table(ExampleDb.seurat$my.clusters2),decreasing=T))
for (i in mycluster.names){
  for(j in c("B6","SLE.yaa")){
    df<-ExampleDb.seurat[ExampleDb.seurat$my.clusters2==i&ExampleDb.seurat$condition==j,] #gating on given cluster
    freq.tab<-as.data.frame(table(df$clone_id))  #ms_cdr3 frequency table for given cluster
    freq.tab$condition<-j
    freq.tab$my.clusters2<-i
    df.list[[paste0(i,j)]]<-freq.tab #add to list for each cluster
  }
  df2<-rbindlist(df.list)
}
df<-rbindlist(df.list)
ggbarplot(df, x = "my.clusters2", y = "Freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.signif", method="wilcox.test",label.y=20,size=2)+
  labs(y = "Clone Size")+theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(limits = mycluster.names) #+scale_y_log10()
ggsave2("bar.dot.cluster.clone.freq.png",width=5, height=4,device="png")
#mut freq
ggbarplot(ExampleDb.seurat, x = "my.clusters2", y = "db_obs.mu_freq", color="condition", add = c("mean","jitter"),
          position = position_dodge(0.8), legend="right",palette=c("black","red"),xlab=F)+
  stat_compare_means(aes(group = condition),label = "p.signif", method="wilcox.test",size=2)+
  labs(y = "Mutation Frequency")+theme(legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(limits = mycluster.names) #+scale_y_log10()
ggsave2("bar.dot.cluster.mu.freq.png",width=5, height=4,device="png")




#my own plots
dir.create("plots")
feature.list<-names(DATA@meta.data)[sapply(DATA@meta.data, is.numeric)]
for(j in feature.list){
  p<-FeaturePlot(DATA, features = j, split.by = "condition",pt.size=0.001,combine=F)
  for(i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoAxes()+
      theme(panel.border = element_rect(colour = "black", size=1),
            axis.title.y.right=element_blank(),plot.title=element_blank(),
            axis.line=element_blank())
  }
  cowplot::plot_grid(plotlist=p,ncol=2)
  ggsave2(paste0("./plots/umap.",j,".bycondition.png"),width=10, height=4,device="png")
}
for(j in feature.list){
  p<-try(VlnPlot(DATA2, features = j, split.by = "condition",  
                 group.by = "my.clusters2",cols=c("grey","red"),pt.size = 0.0001))
  if(is(p, "try-error")){
    cells.keep<-rownames(DATA@meta.data)[!is.na(DATA@meta.data[[j]])]
    DATA2<-subset(DATA,cells=cells.keep)
    p<-VlnPlot(DATA2, features = j, split.by = "condition",  
               group.by = "my.clusters2",cols=c("grey","red"),pt.size = 0.0001)+
      stat_compare_means(label = "p.signif", method="t.test")
  }
  p
  ggsave2(paste0("./plots/vln.",j,".bycondition.png"),width=12, height=8,device="png")
}
colnames(DATA@meta.data)
feature.list<-names(DATA@meta.data)[!sapply(DATA@meta.data, is.numeric)]
for(j in feature.list){
  if(length(unique(DATA@meta.data[[j]]))<15 & length(unique(DATA@meta.data[[j]]))>1){
    Idents(DATA)<-j
    p<-DimPlot(DATA, split.by = "condition",pt.size=0.001,combine=F)
    for(i in 1:length(p)) {
      p[[i]] <- p[[i]] + NoAxes()+
        theme(panel.border = element_rect(colour = "black", size=1),
              axis.title.y.right=element_blank(),plot.title=element_blank(),
              axis.line=element_blank())
    }
    cowplot::plot_grid(plotlist=p,ncol=2)
    ggsave2(paste0("./plots/umap.",j,".bycondition.png"),width=10, height=4,device="png")
  }
}

