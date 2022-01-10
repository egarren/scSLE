rm(list=ls())
library(tigger)
library(dplyr)
barcoder2<-function(df, prefix){
  df$ms.barcoder<-prefix
  df$orig.barcode<-df$cell_id
  # df$barcode <- gsub("-1", "", df$cell_id)
  df$barcode <- paste0(prefix, df$df$cell_id)
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
AIRRDb<-left_join(ExampleDb,metadata,by="ms.barcoder")

gl<-readIgFasta("./ref_genomes/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta")
# Detect novel alleles
novel <- findNovelAlleles(AIRRDb, gl, nproc=1)
# Extract and view the rows that contain successful novel allele calls
novel_rows <- selectNovel(novel)
# Plot evidence of the first (and only) novel allele from the example data
novel_row <- which(!is.na(novel$polymorphism_call))[1]
png("novel.alleles.png")
plotNovel(AIRRDb, novel[novel_row, ])
dev.off()

# Infer the individual's genotype, using only unmutated sequences and checking
# for the use of the novel alleles inferred in the earlier step.
geno <- inferGenotype(AIRRDb, germline_db=gl, novel=novel,
                      find_unmutated=TRUE)
# Save the genotype sequences to a vector
genotype_db <- genotypeFasta(geno, gl, novel)
# Visualize the genotype and sequence counts
print(geno)
# Make a colorful visualization. Bars indicate presence, not proportion.
png("genotypes.png",height=8,width=4,units="in",res=200)
plotGenotype(geno, text_size = 10)
dev.off()

# Infer the individual's genotype, using the bayesian method
geno_bayesian <- inferGenotypeBayesian(AIRRDb, germline_db=gl, 
                                       novel=novel, find_unmutated=TRUE)
# Visualize the genotype and sequence counts
print(geno_bayesian)
# Make a colorful visualization. Bars indicate presence, not proportion.
png("genotypes.bayesian.png",height=8,width=4,units="in",res=200)
plotGenotype(geno_bayesian, text_size=10)
dev.off()





