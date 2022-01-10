rm(list=ls())
# Load required packages
library(alakazam)
library(dplyr)
library(scales)
library(cowplot)
library(patchwork)
library(ggpubr)
library(igraph)

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

# Quantify usage at the gene level
gene <- countGenes(ExampleDb, gene="v_call", groups="condition", mode="gene")
head(gene, n=4)

# Assign sorted levels and subset to IGHV1
ighv1 <- gene %>%
  mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
  filter(getFamily(gene) == "IGHV1")
# Plot V gene usage in the IGHV1 family by sample
ggplot(ighv1, aes(x=gene, y=seq_freq)) +
  theme_bw() +
  ggtitle("IGHV1 Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous(labels=percent) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=condition), size=5, alpha=0.8)
ggsave2("IGHV1.usage.png",width=9, height=4,device="png")

# Quantify V family usage by sample
family <- countGenes(ExampleDb, gene="v_call", groups="condition", mode="family")
# Plot V family usage by sample
ggplot(family, aes(x=gene, y=seq_freq)) +
  theme_bw() +
  ggtitle("Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous(labels=percent) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=condition), size=5, alpha=0.8)
ggsave2("IGHV.usage.png",width=9, height=4,device="png")

#gene abundance by clone count
family <- countGenes(ExampleDb, gene="v_call", groups=c("condition", "c_call"), 
                     clone="ms_clone", mode="family")
head(family, n=4)
# Subset to IGHM and IGHG for plotting
family <- filter(family, c_call %in% c("IGHM", "IGHG2B"))
# Plot V family clonal usage by sample and isotype
ggplot(family, aes(x=gene, y=clone_freq)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous(labels=percent) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=condition), size=5, alpha=0.8) +
  facet_grid(. ~ c_call)
ggsave2("IGHV.C.usage.png",width=9, height=4,device="png")

#gene abundance by copy number
# Calculate V family copy numbers by sample and isotype
family <- countGenes(ExampleDb, gene="v_call", groups=c("condition", "c_call"), 
                     mode="family", copy="consensus_count")
head(family, n=4)
# Subset to IGHM and IGHG for plotting
family <- filter(family, c_call %in% c("IGHM", "IGHG2B"))
# Plot V family copy abundance by sample and isotype
ggplot(family, aes(x=gene, y=copy_freq)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous(labels=percent) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=condition), size=5, alpha=0.8) +
  facet_grid(. ~ c_call)
ggsave2("IGHV.C.copy.png",width=9, height=4,device="png")

##Amino acid anlaysis
df.list<-list()
for(i in unique(ExampleDb$condition)){
  db<-ExampleDb[ExampleDb$condition==i,]
  #calculate aa properties
  df.list[[i]]<- aminoAcidProperties(db, seq="junction", nt=TRUE, trim=TRUE,label="cdr3")
}
db_props<-do.call(rbind,df.list)

# The full set of properties are calculated by default
dplyr::select(db_props[1:3, ], starts_with("cdr3"))

# Generate plots for all four of the properties
plot.list<-list()
for (i in c("cdr3_aa_length","cdr3_aa_gravy","cdr3_aa_bulk","cdr3_aa_aliphatic","cdr3_aa_polarity","cdr3_aa_charge")){
  plot.list[[i]] <- ggplot(db_props, aes_string(x="c_call", y=i,fill="condition")) +
    xlab("Isotype") + ylab(i) + theme_bw() + 
    stat_compare_means(aes(group = condition),label = "p.signif", method="t.test")+
    # scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot()#aes(fill=c_call)
}
for (i in c("cdr3_aa_basic","cdr3_aa_acidic","cdr3_aa_aromatic")){
  plot.list[[i]] <- ggplot(db_props, aes_string(x="c_call", y=i,fill="condition")) +
    xlab("Isotype") + ylab(i) +theme_bw() + 
    scale_y_continuous(labels=scales::percent) +
    stat_compare_means(aes(group = condition),label = "p.signif", method="t.test")+
    # scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot()#aes(fill=c_call)
}

# Plot in a 2x2 grid
wrap_plots(plot.list,guides="collect")& theme(legend.position = 'bottom')
ggsave2("aa.analysis.png",width=15, height=10,device="png")

##Diversity
# Partitions the data based on the sample column
clones <- countClones(ExampleDb, clone="ms_clone",group="condition")
head(clones, 5)
clones <- countClones(ExampleDb, clone="ms_clone",group=c("condition","c_call"),copy="consensus_count")
head(clones, 5)
curve <- estimateAbundance(ExampleDb, group="condition", ci=0.95, nboot=200, clone="ms_clone")
plot(curve)
ggsave2("clonal.abundance.png",width=4, height=3,device="png")

#diversity curve
# Compare diversity curve across values in the "sample" column
# q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 increments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 200 resampling realizations are performed (nboot=200)
sample_curve <- alphaDiversity(ExampleDb, group="condition", clone="ms_clone",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=200)
# Compare diversity curve across values in the c_call column
# Analyse is restricted to c_call values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
df<-ExampleDb[ExampleDb$c_call %in% c("IGHA","IGHD","IGHG1","IGHG2B","IGHG2C","IGH3","IGHM"),]
isotype_curve <- alphaDiversity(df, group="c_call", clone="ms_clone",
                                min_q=0, max_q=4, step_q=0.1,
                                ci=0.95, nboot=200)
# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
plot(sample_curve, colors=c("B6"="black", "SLE.yaa"="red"),legend_title="")
ggsave2("sample.diversity.curve.png",width=4, height=3,device="png")
plot(isotype_curve)
ggsave2("isotype.diversity.curve.png",width=4, height=3,device="png")


#diversity test
# Test diversity at q=0, q=1 and q=2 (equivalent to species richness, Shannon entropy, 
# Simpson's index) across values in the sample_id column
# 200 bootstrap realizations are performed (nboot=200)
sample_test <- alphaDiversity(ExampleDb, group="condition", min_q=0, max_q=2, step_q=1, nboot=200, clone="ms_clone")
isotype_test <- alphaDiversity(df, group="c_call", min_q=0, max_q=2, step_q=1, nboot=200, clone="ms_clone")
# Print P-value table
print(sample_test@tests)
print(isotype_test@tests)
# Plot results at q=0 and q=2
# Plot the mean and standard deviations at q=0 and q=2
plot(isotype_test, 0, colors=IG_COLORS, legend_title="Isotype")
plot(isotype_test, 2, colors=IG_COLORS, legend_title="Isotype")
plot(sample_test,2)
ggsave2("sample.diversity.test.png",width=3, height=2.5,device="png")

#Ig lineage tree
head(sort(table(ExampleDb$ms_clone), decreasing=T))
sub_db <- subset(ExampleDb, ms_clone == "m233_2228_51")
# This example data set does not have ragged ends
# Preprocess clone without ragged end masking (default)
clone <- makeChangeoClone(sub_db, text_fields=c("condition", "c_call"),num_fields="consensus_count",clone="ms_clone")
# Show combined annotations
clone@data[, c("condition", "c_call", "consensus_count")]
# Run PHYLIP and parse output
phylip_exec <- "./phylip-3.697/exe/dnapars"
graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
# The graph has shared annotations for the clone
data.frame(clone_id=graph$clone,
           junction_length=graph$junc_len,
           v_gene=graph$v_gene,
           j_gene=graph$j_gene)
# The vertices have sequence specific annotations
data.frame(sequence_id=V(graph)$name, 
           c_call=V(graph)$c_call)
# Plot graph with defaults
plot(graph)
# Modify graph and plot attributes
V(graph)$color <- "steelblue"
V(graph)$color[V(graph)$name == "Germline"] <- "black"
V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
V(graph)$label <- V(graph)$c_call
E(graph)$label <- ""
# Remove large default margins
png("imm.tree.png")
par(mar=c(0, 0, 0, 0) + 0.1)
# Plot graph
plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=40)
# Add legend
legend("topleft", c("Germline", "Inferred", "Sample"), 
       fill=c("black", "white", "steelblue"), cex=0.75)
dev.off()
# convert to phylo
phylo <- graphToPhylo(graph)
#plot using ape
png("imm.phylo.png")
plot(phylo, show.node.label=TRUE)
dev.off()

#batch process clones
# Preprocess clones
clones2 <- ExampleDb %>%
  group_by(ms_clone) %>%
  do(CHANGEO=makeChangeoClone(., text_fields=c("condition", "c_call"), 
                              num_fields="consensus_count",clone="ms_clone"))
# Build lineages
fun<- function (...) {
  return(tryCatch(buildPhylipLineage(...), error=function(e) NULL))
}
graphs <- lapply(clones2$CHANGEO,fun,phylip_exec=phylip_exec, rm_temp=TRUE)
# Note, clones with only a single sequence will not be processed.
# A warning will be generated and NULL will be returned by buildPhylipLineage
# These entries may be removed for clarity
graphs[sapply(graphs, is.null)] <- NULL
# The set of tree may then be subset by node count for further 
# analysis, if desired.
graphs2 <- graphs[sapply(graphs, vcount) >= 5]
#plot
dir.create("trees")
for(i in 1:length(graphs2)){
  graph<-graphs2[[i]]
  # Modify graph and plot attributes
  V(graph)$color <- "steelblue"
  V(graph)$color[V(graph)$name == "Germline"] <- "black"
  V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
  V(graph)$label <- V(graph)$c_call
  E(graph)$label <- ""
  # Remove large default margins
  png(paste0("./trees/imm.tree",i,".png"))
  par(mar=c(0, 0, 0, 0) + 0.1)
  # Plot graph
  plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
       vertex.label.color="black", vertex.size=40)
  # Add legend
  legend("topleft", c("Germline", "Inferred", "Sample"), 
         fill=c("black", "white", "steelblue"), cex=0.75)
  dev.off()
  # convert to phylo
  phylo <- graphToPhylo(graph)
  #plot using ape
  png(paste0("./trees/imm.phylo",i,".png"))
  plot(phylo, show.node.label=TRUE)
  dev.off()
}