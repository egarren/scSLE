rm(list=ls())
###Pathway analysis
library(org.Mm.eg.db)
library(tidyverse)
library(RDAVIDWebService)
library(Seurat)
library(cowplot)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(topGO)
library(scde)
library(biomaRt)
library(GO.db)
library(DBI)
# library(VISION)
library(msigdbr)
# library(msigdb)
# library(goseq)
# library(nicethings)
library(DOSE)
library(ggpubr)

#load DEs
dir.create("./gsea")
setwd("./gsea")
load("../graphed.RData")

m_df = msigdbr(species = "Mus musculus")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
pathways<-m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
msig.df<-m_df %>% dplyr::select(gs_name,entrez_gene)

#load DE df
DE.df<-SLE.obj.combined.autoimmune.response
head(DE.df,n=20)

if(!("p_val_adj" %in% colnames(DE.df))){
  if("ttest" %in% colnames(DE.df) ){DE.df$p_val_adj<-p.adjust(DE.df$ttest,method="BH")}
  if("p_val" %in% colnames(DE.df) ){DE.df$p_val_adj<-p.adjust(DE.df$p_val,method="BH")}
  if("pval" %in% colnames(DE.df) ){DE.df$p_val_adj<-DE.df$pval}
}
if(!("log2FC" %in% colnames(DE.df))){
  if("delta" %in% colnames(DE.df)){DE.df$log2FC<-DE.df$delta}
  if("rho" %in% colnames(DE.df)){DE.df$log2FC<-DE.df$rho}
  if("avg_log2FC" %in% colnames(DE.df)){DE.df$log2FC<-DE.df$avg_log2FC}
  if("avg_logFC" %in% colnames(DE.df)){DE.df$log2FC<-log2(exp(DE.df$avg_logFC))}
}
if(sum(DE.df$log2FC)!=0){
  ## FSGEA
  df<-DE.df
  df<-df[df$p_val_adj<0.05,] #select for sig genes, abs(df$log2FC)>0.2&
  if(nrow(df)>50){
    if(!("gene" %in% colnames(df))){df$gene<-rownames(df)}
    ranks <- df$log2FC
    names(ranks) <- df$gene
    png(paste0(i,".ranks.png"),width=4,height=4,units="in",res=200)
    barplot(sort(ranks, decreasing = T))
    dev.off()
    for(h in c("pathways","CTD")){
      path<-get(h)
      fgseaRes <- fgsea(path, ranks, minSize=15, maxSize = 500)#, nperm=1000
      if(!is.null(fgseaRes)){if(nrow(fgseaRes)>1){
        #enrichment plot
        head(fgseaRes[order(padj, -abs(NES)), ], n=15)
        plotEnrichment(path[[fgseaRes[order(padj, -abs(NES)), ]$pathway[132]]], ranks)+
          ggtitle(fgseaRes[order(padj, -abs(NES)), ]$pathway[132])+
          theme(plot.title = element_text(size=5),axis.title=element_text(size=5),axis.text=element_text(size=5))#plot top pathway enrichment
        ggsave2(paste0(i,".",h,".topenrichment.png"),width=2, height=2,device="png")
        #gsea table
        topUp <- as.data.frame(head(arrange(fgseaRes %>% filter(ES > 0), padj),n=5))
        topDown <-  as.data.frame(head(arrange(fgseaRes %>% filter(ES < 0), padj),n=5))
        topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
        png(paste0(i,".",h,".gsea.table.png"),width=10,height=4,units="in",res=200)
        plotGseaTable(path[topPathways$pathway], ranks, fgseaRes, gseaParam = 0.5)
        dev.off()
        write.csv(fgseaRes[order(padj, -abs(NES)), ][,1:7],file=paste0(h,".",i,".gsea.results.csv"))
        
      }}
    }
  }
  
  
  ##clusterProfiler (https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)
  df<-DE.df
  if(!(i %in% c(ls(pattern=".size.correl"),ls(pattern="gene.correl")))){df$gene<-rownames(df)}
  genes.meta<-bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db") #pull in entrezID, could also use bitr_kegg, see: keytypes(org.Mm.eg.db)
  df<-cbind(df,genes.meta[match(df$gene, genes.meta$SYMBOL),])
  if(i %in% c(ls(pattern=".size.correl"),ls(pattern="gene.correl"))){
    sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.001])
    if(length(sigGenes)<50){sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.01])} 
    if(length(sigGenes)<50){sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.05])} 
  }else{
    sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.00001 &abs(df$log2FC) > 0.3 ]) #select sig genes, abs(df$log2FC) 
    if(length(sigGenes)<50){sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.0001 &abs(df$log2FC) > 0.25 ])} 
    if(length(sigGenes)<50){sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.001 &abs(df$log2FC) > 0.2 ])} 
    if(length(sigGenes)<50){sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.01 &abs(df$log2FC) > 0.1 ])} 
  }
  df2<-df[df$p_val_adj<0.01 & !is.na(df$ENTREZID),]
  ranks <- df2$log2FC
  names(ranks) <- df2$ENTREZID
  geneList <- sort(ranks, decreasing = TRUE)
  #topGO
  ggo.table<-groupGO(gene=sigGenes,OrgDb=org.Mm.eg.db,ont="BP",level=4,readable=T) #BP (biological process), MF (molecular function), CC (cellular componnent)
  if(!is.null(ggo.table)){if(dim(ggo.table)[1]!=0 & sum(ggo.table$Count) != 0){
    barplot(ggo.table, drop=TRUE, showCategory=6)+ggtitle(i)
    ggsave2(paste0(i,".gGO.bar.png"),width=8, height=2,device="png")
  }}
  for(j in c("eGO","gsea","kk","kk2","mkk","mkk2","davidKEGG","davidBP","egmt","egmt2","ctd","ctd2")){ #
    if(j=="eGO"){tab<-try(enrichGO(gene=sigGenes,OrgDb=org.Mm.eg.db,ont="BP", readable=T))}
    if(j=="gsea"){tab<-try(gseGO(geneList,OrgDb=org.Mm.eg.db,ont="BP",pvalueCutoff = 0.05))}
    if(j=="kk"){tab <- try(enrichKEGG(gene = sigGenes,organism = 'mmu'))}#search_kegg_organism('mmu', by='kegg_code')}
    if(j=="kk2"){tab<-try(gseKEGG(geneList,organism="mmu",pvalueCutoff=0.05))} #KEGG gsea}
    if(j=="mkk"){tab<-try(enrichMKEGG(sigGenes,organism="mmu"))}
    if(j=="mkk2"){tab<-try(gseMKEGG(geneList,organism="mmu",pvalueCutoff=0.05)) }
    if(j=="davidKEGG"){tab<-try(enrichDAVID(sigGenes,annotation="KEGG_PATHWAY",david.user="XX")) }
    if(j=="davidBP"){tab<-try(enrichDAVID(sigGenes,annotation="GOTERM_BP_FAT",david.user="XX"))}
    if(j=="egmt"){tab<-try(enricher(sigGenes,TERM2GENE = msig.df))}
    if(j=="egmt2"){tab<-try(GSEA(geneList,TERM2GENE = msig.df,pvalueCutoff=0.05))}
    if(j=="ctd"){tab<-try(enricher(sigGenes,TERM2GENE = CTD.entrez))}
    if(j=="ctd2"){tab<-try(GSEA(geneList,TERM2GENE = CTD.entrez,pvalueCutoff=0.05))}
    if(!is.null(tab)&!is(tab,"try-error")){if(nrow(tab)>1){ 
      if(j %in% c("eGO","gsea")){tab<-setReadable(tab,OrgDb=org.Mm.eg.db)}else{tab<-setReadable(tab,OrgDb=org.Mm.eg.db,keyType="ENTREZID")}
      write.csv(tab,file=paste0(i,".",j,".tab.csv"))
      if(j %in% c("eGO","kk","mkk","davidKEGG")){
        barplot(tab, drop=TRUE, showCategory=6)+ggtitle(i)
        ggsave2(paste0(i,".",j,".bar.png"),width=8, height=2,device="png")
      }
      clusterProfiler::dotplot(tab, showCategory=10)+ggtitle(i)+  
        theme(legend.position="bottom",legend.direction = "horizontal", legend.box = "vertical",
              axis.text.x=element_text(angle=45,hjust=1))
      ggsave2(paste0(i,".",j,".dot.png"),width=5, height=5,device="png")
      clusterProfiler::dotplot(tab, showCategory=12,font.size=7)+ggtitle(i)+  
        guides(size=F)+scale_color_continuous(name="P-val",guide=guide_colorbar(reverse=TRUE))
      ggsave2(paste0(i,".",j,".dot2.png"),width=4.25, height=2,device="png")
      if (length(tab@termsim) == 0) {
        emapplot(pairwise_termsim(tab),showCategory=15)+ggtitle(i)
      }else{emapplot(tab,showCategory=15)+ggtitle(i)}
      ggsave2(paste0(i,".",j,".emap.png"),width=10, height=9,device="png")
      clusterProfiler::cnetplot(tab, categorySize="pvalue", foldChange=geneList)+
        ggtitle(i)+scale_color_gradient2(name="LogFC")
      ggsave2(paste0(i,".",j,".cnet.png"),width=12, height=10,device="png")
      clusterProfiler::cnetplot(tab, categorySize="pvalue", foldChange=geneList,node_label="gene")+
        guides(size=F)+scale_color_gradient2(name="LogFC")
      ggsave2(paste0(i,".",j,".cnet2.png"),width=12, height=10,device="png")
      
      if(j %in% c("eGO","gsea")){
        png(paste0(i,".",j,".graph.png"),width=20,height=15,units="in",res=200)
        clusterProfiler::plotGOgraph(tab)
        dev.off()
      }
      if(j %in% c("gsea","kk2","egmt2","ctd2")){
        clusterProfiler::gseaplot(tab,geneSetID=tab$ID[1],title=tab$Description[1])
        ggsave2(paste0(i,".",j,".plot.png"),width=4, height=5,device="png")
      }
      if(j %in% c("kk","kk2","davidKEGG")){
        pathview(gene.data = geneList, pathway.id = tab$ID[2], species = "mmu", out.suffix=paste0(i,".",j,".pathview"))
      }
    }}
  }
}

for(i in c("size.correl","gene.correl")){
  dir.create(paste0("./",i))
  file.names<-list.files(pattern=i)
  file.rename(from=paste0("./",file.names),to=paste0("./",i,"/",file.names))
}


#Cluster GO dotplot
for(i in c("condition","my.clusters")){ #"GLIPH.exp",
  Idents(SLE.obj.combined)<-i
  clust.mark<- FindAllMarkers(object = SLE.obj.combined, test.use = "MAST",only.pos = TRUE,
                              min.pct=0,logfc.threshold = 0.1,return.thresh=0.05)
  genes.meta<-bitr(clust.mark$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db") #pull in entrezID, could also use bitr_kegg, see: keytypes(org.Mm.eg.db)
  clust.mark<-cbind(clust.mark,genes.meta[match(clust.mark$gene, genes.meta$SYMBOL),])
  clust.mark<-clust.mark[!is.na(clust.mark$ENTREZID),]
  clust.mark[] <- lapply(clust.mark, as.character)
  clust.mark.genes<-clust.mark %>% split(x = .$ENTREZID, f = .$cluster) 
  ck<-compareCluster(geneClusters=clust.mark.genes,fun="enrichKEGG",organism="mmu") 
  clusterProfiler::dotplot(ck,showCategory=6)+ggtitle(i)+  theme(legend.direction = "vertical", legend.box = "horizontal")
  ggsave2(paste0(i,".cluster.KEGG.dot.png"),width=9, height=4,device="png")
  c.ggo<-compareCluster(geneCluster=clust.mark.genes,fun="groupGO",OrgDb=org.Mm.eg.db,ont="BP",level=4,readable=T) 
  clusterProfiler::dotplot(c.ggo,showCategory=6)+ggtitle(i)+  theme(legend.direction = "vertical", legend.box = "horizontal")
  ggsave2(paste0(i,".cluster.gGO.dot.png"),width=9, height=4,device="png")
  c.ego<-compareCluster(geneCluster=clust.mark.genes,fun="enrichGO",OrgDb=org.Mm.eg.db,ont="BP", readable=T)
  clusterProfiler::dotplot(c.ego,showCategory=10,font.size=7)+ggtitle(i) +
    guides(size=F)+scale_color_continuous(name="P-val",trans="log",guide=guide_colorbar(reverse=TRUE),
                                          breaks=c(1e-5,1e-10,1e-15,1e-20,1e-25,1e-30),
                                          labels=format(c(1e-5,1e-10,1e-15,1e-20,1e-25,1e-30)))+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave2(paste0(i,".cluster.eGO.dot.png"),width=4.5, height=5,device="png")
  c.msigdb<-compareCluster(geneCluster=clust.mark.genes,fun="enricher",TERM2GENE = msig.df)
  clusterProfiler::dotplot(c.msigdb,showCategory=10,font.size=7)+ggtitle(i) +
    guides(size=F)+scale_color_continuous(name="P-val",trans="log",guide=guide_colorbar(reverse=TRUE),
                                          breaks=c(1e-5,1e-10,1e-15,1e-20,1e-25,1e-30),
                                          labels=format(c(1e-5,1e-10,1e-15,1e-20,1e-25,1e-30)))+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave2(paste0(i,".cluster.msigdb.dot.png"),width=6.5, height=6,device="png")
  c.ctd<-compareCluster(geneCluster=clust.mark.genes,fun="enricher",TERM2GENE = CTD.entrez)
  clusterProfiler::dotplot(c.ctd,showCategory=6)+ggtitle(i)+  theme(legend.direction = "vertical", legend.box = "horizontal")
  ggsave2(paste0(i,".cluster.ctd.dot.png"),width=14, height=4,device="png")
  save.image(paste0("temp.",i,".clustdots.RData"))
}

