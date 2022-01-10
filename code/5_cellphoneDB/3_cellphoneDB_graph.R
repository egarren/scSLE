rm(list=ls())
library(ggplot2)
library(readr)
library(tidyverse)
library(viridis)
library(data.table)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(dplyr)

#import cellphoneDB output
b6.means<-read_tsv("./B6_cellphoneDB/out/means.txt")
b6.p<-read_tsv("./B6_cellphoneDB/out/pvalues.txt")
SLE.yaa.means<-read_tsv("./SLE.yaa_cellphoneDB/out/means.txt")
SLE.yaa.p<-read_tsv("./SLE.yaa_cellphoneDB/out/pvalues.txt")
# 
test<-b6.means
test<-test[test$interacting_pair =="CSF1_SIRPA",]
test<-dplyr::select(test,contains("PC"))

#create df
b6.df<-gather(b6.means, colnames(b6.means)[12:ncol(b6.means)], key="Pair", value = "Score")
b6.df2<-gather(b6.p, colnames(b6.p)[12:ncol(b6.p)], key="Pair", value = "p-val")
df<-full_join(b6.df,b6.df2)
df$condition<-"B6"
SLE.yaa.df<-gather(SLE.yaa.means, colnames(SLE.yaa.means)[12:ncol(SLE.yaa.means)], key="Pair", value = "Score")
SLE.yaa.df2<-gather(SLE.yaa.p, colnames(SLE.yaa.p)[12:ncol(SLE.yaa.p)], key="Pair", value = "p-val")
df2<-full_join(SLE.yaa.df,SLE.yaa.df2)
df2$condition<-"SLE.yaa"
keep<-intersect(df$interacting_pair,df2$interacting_pair)
df3<-rbind(df[df$interacting_pair %in% keep,],df2[df2$interacting_pair %in% keep,])
df3$logP= -log(df3$`p-val`)
df3$logP[is.infinite(df3$logP)]<-6
df3$receptor.cell<-gsub("\\|.*","",df3$Pair)
df3$ligand.cell<-gsub(".*\\|","",df3$Pair)
df3<-df3[order(-df3$Score,df3$`p-val`),]


#scatter plot
ggplot(df3, aes(x=Score,y=logP,color=condition))+geom_point()
sig<-df3[df3$`p-val`<0.05,]
sig$label=paste0(sig$interacting_pair," in ",sig$Pair)
sig2<-sig[,c("label","condition","Score")]
scores<-spread(sig2,key=condition,value=Score)
scores$delta=abs(scores$B6-scores$SLE.yaa)
df4<-scores[order(-scores$delta),]
ggplot(scores, aes(x=SLE.yaa,y=B6,color=delta))+#
  geom_point()+
  scale_color_viridis(option = "D",direction=-1)+
  theme_classic()+
  theme(legend.position="none")+
  geom_text_repel(size=3,aes(label=ifelse(delta>0.7,as.character(label),"")),
                  box.padding = 0.1, segment.size=0.2,max.overlaps=50)
ggsave2("scores.png",width=4, height=4,device="png")

#dotplot
head(unique(df3$interacting_pair),n=50)
unique(grep("VSIR",df3$interacting_pair,value=T))
ssig.interact=c("CD40_CD40LG","FCER2_CR2","CSF1_SIRPA","SIRPA_CD47","CEACAM1_CD209","VCAM1_a4b1 complex","SELL_SELPLG")
df4<-df3[df3$interacting_pair %in% sig.interact,]
df4<-df4[df4$`p-val`<0.05,]
ggplot(data = df4, mapping = aes_string(x = "receptor.cell", y = "interacting_pair")) + 
  geom_point(mapping = aes_string(size = "logP", fill = "Score"),pch=21,color="black") +
  facet_grid(ligand.cell ~ condition)+
  scale_fill_viridis(option="H",guide = guide_colorbar(frame.colour = "black",frame.linewidth=1,
                                                       ticks.colour="black",ticks.linewidth=1,
                                                       title.position="top",title.hjust=0.5))+
  theme(legend.position="bottom",strip.background = element_rect(colour="black", fill="white",linetype="solid"),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1),
        panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA,   colour = "grey20"),
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5),linetype="dashed"),
        panel.grid.major=element_line(linetype="dashed"),
        legend.key = element_rect(fill = "white", colour = NA))+
  labs(x="Receptor-expressing cell",y="Receptor-Ligand pair")+
  guides(size = guide_legend(title="-logP",title.position="top", title.hjust = 0.5,override.aes=list(fill="black")))+
  scale_size_continuous(limits=c(2,6),breaks=c(2,4,6),labels = c("2","4","6"))
ggsave2("dotplot2.png",width=8, height=16,device="png")

#dotplot for single cell type
for(i in unique(df3$ligand.cell)){
  df5<-df3[df3$ligand.cell==i,]
  df5<-df5[df5$`p-val`<0.05,]
  ggplot(data = df5, mapping = aes_string(x = "condition", y = "interacting_pair")) + 
    geom_point(mapping = aes_string(size = "logP", fill = "Score"),pch=21,color="black") +
    facet_grid(. ~receptor.cell)+
    scale_fill_viridis(option="H",guide = guide_colorbar(frame.colour = "black",frame.linewidth=1,
                                                         ticks.colour="black",ticks.linewidth=1,
                                                         title.position="top",title.hjust=0.5))+
    theme(legend.position="bottom",strip.background = element_rect(colour="black", fill="white",linetype="solid"),
          axis.text.x=element_text(angle=45,vjust=1,hjust=1),
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,   colour = "grey20"),
          panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5),linetype="dashed"),
          panel.grid.major=element_line(linetype="dashed"),
          legend.key = element_rect(fill = "white", colour = NA))+
    labs(x="",y="Ligand-Receptor pair")+
    guides(size = guide_legend(title="-logP",title.position="top", title.hjust = 0.5,override.aes=list(fill="black")))+
    scale_size_continuous(limits=c(2,6),breaks=c(2,4,6),labels = c("2","4","6"))
  ggsave2(paste0(i,".ligand.dot.png"),width=8, height=length(unique(df5$interacting_pair))/4,device="png")
}

#heatmap
sig<-df3[df3$`p-val`<0.05,]
sig$interaction_condition<-paste0(sig$interacting_pair,".",sig$condition)
sig2<-sig[,c("Pair","Score","interaction_condition")]
matrix<-as.data.frame(spread(sig2,key=interaction_condition,value=Score))
rownames(matrix)<-matrix$Pair
matrix$Pair<-NULL
matrix[is.na(matrix)]<-0
pheatmap(matrix)

#add heatmap metadata
annot.col<-data.frame(interaction_condition=colnames(matrix))
annot.col<-merge(annot.col,unique(sig[,c("interaction_condition","condition","secreted","is_integrin")]),by="interaction_condition")
names(annot.col)[names(annot.col) == "condition"] <- "Genotype"
names(annot.col)[names(annot.col) == "secreted"] <- "Secreted"
names(annot.col)[names(annot.col) == "is_integrin"] <- "Integrin"
annot.col$Secreted<-ifelse(annot.col$Secreted==T,"Yes","No")
annot.col$Integrin<-ifelse(annot.col$Integrin==T,"Yes","No")
row.names(annot.col) <- annot.col$interaction_condition
annot.col$interaction_condition <- NULL
p1<-pheatmap(matrix,annotation_col=annot.col,fontsize_col=3,fontsize_row=5,annotation_names_col=T,
             annotation_colors=list(Genotype=c("B6"="grey","SLE.yaa"="red"),
                                    Integrin=c("No"="blue4","Yes"="forestgreen"),
                                    Secreted=c("No"="blue4","Yes"="forestgreen")))
save_pheatmap_png <- function(x, filename, width=5000, height=3400, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(p1, "scores.heatmap.png")
 
##interaction heatmap
for(i in c("B6","SLE.yaa")){
  df<-read_tsv(paste0("./",i,"_cellphoneDB/out/count_network.txt"))
  df<-as.data.frame(spread(df, key=TARGET, value=count))
  rownames(df)<-df$SOURCE
  df$SOURCE<-NULL
  p1<-pheatmap(df,cluster_cols=F,cluster_rows=F, color=hcl.colors(50,"Oslo"))
  save_pheatmap_png <- function(x, filename, width=1100, height=1000, res = 300) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  save_pheatmap_png(p1, paste0(i,".interaction.heatmap.png"))
}
