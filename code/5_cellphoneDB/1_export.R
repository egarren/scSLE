require(EWCE)
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
library(Seurat)
library(plyr)

load("graphed.RData")

# Basic function to convert mouse to human gene names
alldata<-SLE.obj.combined 
alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000)
allgenes <- rownames(alldata)
matrix1 <- as.data.frame(alldata@assays$RNA@data)
matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]

### If you are using a mouse data, then its needed to convert the gene names to human orthologs
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = rownames(alldata@assays$RNA@data) , mart = mouse, 
                 attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))
matrix2 <- matrix1[match(genesV2$MGI.symbol,rownames(alldata),nomatch=F),] #only keep human orthologs
matrix2$gene <- genesV2$Gene.stable.ID
#dealing with duplicates
sum(duplicated(matrix2$gene)) #count number of duplicate orthologs
matrix3<-ddply(matrix2,"gene",numcolwise(sum)) #sum counts for multiple mapped orthologs 
rownames(matrix3) <-  matrix3$gene
matrix3$gene<-NULL
save.image("cellphoneDB.RData")

### Subseting the matrix
B6 <- grepl('B6',alldata@meta.data$condition)
SLE.yaa <- grepl('SLE.yaa',alldata@meta.data$condition)
write.table(cbind(Gene=rownames(matrix3),matrix3[,B6]), 'B6_filtered_hcount.txt',row.names=F,sep='\t',quote=F)
write.table(cbind(Gene=rownames(matrix3),matrix3[,SLE.yaa]), 'SLE.yaa_filtered_hcount.txt',row.names=F,sep='\t',quote=F)


##Write metadata
metadata_B6 <- data.frame(Cell=rownames(alldata@meta.data[grepl('B6',alldata@meta.data$condition),]),cell_type=alldata@meta.data$my.clusters2[grepl('B6',alldata@meta.data$condition)])
metadata_SLE.yaa <- data.frame(Cell=rownames(alldata@meta.data[grepl('SLE.yaa',alldata@meta.data$condition),]),cell_type=alldata@meta.data$my.clusters2[grepl('SLE.yaa',alldata@meta.data$condition)]) ## Just negate grepl('state1',alldata@meta.data$stim),]
write.table(metadata_B6, 'B6_filtered_meta.txt',row.names=F,sep='\t',quote=F)
write.table(metadata_SLE.yaa, 'SLE.yaa_filtered_meta.txt',row.names=F,sep='\t',quote=F)





