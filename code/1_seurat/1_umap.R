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

#load data
metadata<-read.csv("SLE.metadata.csv", header=T) 
metadata[] <- lapply(metadata, as.character)
mice<-metadata$mouse_ID
for (i in mice){
  path<-paste0(i,"/outs/count/filtered_feature_bc_matrix") #define path
  seurat.object<-CreateSeuratObject(counts = Read10X(data.dir = path,gene.column=1) , 
                                    min.cells = 3, min.features  = 200, project = i, assay = "RNA")
  assign(i,seurat.object)
}

##Add Count Matrices 
seurat.list<-lapply(mice,get)
SLE.obj<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=mice,project="SLE")
# add sample metadata
SLE.obj<-AddMetaData(SLE.obj, metadata=rownames(SLE.obj@meta.data),col.name = "cellID")
SLE.obj@meta.data<-left_join(x = SLE.obj@meta.data, y = metadata, by = c("orig.ident"="mouse_ID"),keep=T)
rownames(SLE.obj@meta.data)<-SLE.obj@meta.data$cellID
rm(list=setdiff(ls(), c("SLE.obj", "metadata")))

##QC and Filter
#adding feature (gene) metadata
SLE.obj@assays[["RNA"]]@meta.features$original_ensembl<-rownames(SLE.obj@assays[["RNA"]]@meta.features)
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position", "chromosome_name", 
                      "percentage_gene_gc_content", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "ensembl_gene_id",values=list(rownames(SLE.obj@assays[["RNA"]]@meta.features)),
                  mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="ensembl.org"),useCache=F) 
m <- match(SLE.obj@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
unique(genes.meta$mgi_symbol[genes.meta$gene_biotype=="IG_C_gene"])
SLE.obj@assays[["RNA"]]@meta.features<-cbind(SLE.obj@assays[["RNA"]]@meta.features,genes.meta[m,])
#renaming features to gene symbols
rownames(SLE.obj@assays[["RNA"]]@meta.features) <- make.names(SLE.obj@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(SLE.obj@assays[["RNA"]]@data) <- make.names(SLE.obj@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(SLE.obj@assays[["RNA"]]@counts) <- make.names(SLE.obj@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
# #TCR and Ig identification and filtering
ig_list <- c( "IG_D_gene", "IG_D_pseudogene", "IG_J_gene", "IG_LV_gene", 
             "IG_pseudogene", "IG_V_gene", "IG_V_pseudogene")#"IG_C_gene", "IG_C_pseudogene",
tr_list <-c("TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_C_gene")
SLE.obj <- subset(SLE.obj, features=rownames(SLE.obj[!(SLE.obj@assays[["RNA"]]@meta.features$gene_biotype %in% ig_list),]))
SLE.obj <- subset(SLE.obj, features=rownames(SLE.obj[!(SLE.obj@assays[["RNA"]]@meta.features$gene_biotype %in% tr_list),]))
SLE.obj <- subset(SLE.obj, features=rownames(SLE.obj[!(is.na(SLE.obj@assays[["RNA"]]@meta.features$mgi_symbol)),]))
SLE.obj <- PercentageFeatureSet(SLE.obj, pattern = "^mt", col.name = "percent.mt")
SLE.obj <- PercentageFeatureSet(SLE.obj, pattern = "^Rpl", col.name = "percent.Rpl")
SLE.obj <- PercentageFeatureSet(SLE.obj, pattern = "^Rps", col.name = "percent.Rps")
#cell cycle regression
m_cc<-readRDS("m_cc.rds") #define path
SLE.obj <- CellCycleScoring(SLE.obj, s.features = m_cc$s.genes, g2m.features = m_cc$g2m.genes, set.ident = TRUE)
SLE.obj$CC.Difference <- SLE.obj$S.Score - SLE.obj$G2M.Score
#HSP regression
HSP_genes<-genes.meta[genes.meta$go_id=="GO:0009408",]$mgi_symbol #0034605
SLE.obj <- AddModuleScore(object = SLE.obj,features = list(HSP_genes),name = 'HSP.score')
Idents(SLE.obj)<-"condition"
VlnPlot(object = SLE.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                     "percent.Rpl","percent.Rps","HSP.score1",
                                     "S.Score","G2M.Score"),pt.size=0, ncol = 3)
ggsave2("vln.QC.png",device="png")

dim(SLE.obj)
SLE.obj<- subset(x = SLE.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA>1000 & 
                  percent.mt >  -Inf & percent.mt < 5 & percent.Rpl < 25 & percent.Rps < 25 &
                 S.Score <0.15 & G2M.Score<0.15) 
dim(SLE.obj)

##Integrating by condition
#split by condition 
SLE.obj.list <- SplitObject(SLE.obj, split.by = "condition")
# SCTransform
for (i in 1:length(SLE.obj.list)) {
  SLE.obj.list[[i]] <- SCTransform(SLE.obj.list[[i]], vars.to.regress = c("nCount_RNA", "percent.mt","HSP.score1",
                                                                      "percent.Rpl","percent.Rps")) #", "CC.Difference",
}
#integrate and anchor
SLE.obj.features <- SelectIntegrationFeatures(object.list = SLE.obj.list, nfeatures = 3000)
SLE.obj.list <- PrepSCTIntegration(object.list = SLE.obj.list, anchor.features = SLE.obj.features)
SLE.obj.anchors <- FindIntegrationAnchors(object.list = SLE.obj.list, normalization.method = "SCT",anchor.features = SLE.obj.features)
SLE.obj.combined <- IntegrateData(anchorset = SLE.obj.anchors, normalization.method = "SCT")
SLE.obj.combined <- RunPCA(SLE.obj.combined)
rm(list=setdiff(ls(), c("SLE.obj.combined", "metadata","genes.meta")))

## t-SNE and Clustering
ElbowPlot(object = SLE.obj.combined)
ggsave2("elbow.png",device="png")
SLE.obj.combined <- RunUMAP(SLE.obj.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
SLE.obj.combined <- FindNeighbors(SLE.obj.combined, reduction = "pca", dims = 1:25) #change based on elbow plot
SLE.obj.combined <- FindClusters(SLE.obj.combined, resolution = 0.15) #adjust resolution (bigger=more clusters), initially used 0.2
SLE.obj.combined <- RunTSNE(object = SLE.obj.combined, dims.use = 1:25, do.fast = TRUE) #change number of PCs to use, change perplexity (https://distill.pub/2016/misread-tsne/)
DimPlot(SLE.obj.combined, reduction = "umap",label=T)
for (i in c("condition","Phase","mouse_ID","gender","batch","seurat_clusters")){
  Idents(SLE.obj.combined)<-i
  DimPlot(SLE.obj.combined, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(SLE.obj.combined, features= i,split.by = "condition",pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap.png"),width=10, height=5,device="png")
}
Idents(SLE.obj.combined)<-"seurat_clusters"
save.image("merged.RData")

clusters<-as.data.frame(table(SLE.obj.combined@meta.data$seurat_clusters))
letters<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
for(i in 1:nrow(clusters)){
  SLE.obj.combined@meta.data$cluster_letter[SLE.obj.combined@meta.data$seurat_clusters%in%clusters$Var1[i]]<-letters[i]
}
SLE.obj.combined.sce <- as.SingleCellExperiment(SLE.obj.combined)
save(SLE.obj.combined.sce,file="sce.clusters.RData")
save.image("sce.RData")
