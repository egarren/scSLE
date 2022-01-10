rm(list=ls())
# Import required packages
library(alakazam)
library(shazam)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
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

db <- subset(ExampleDb, c_call %in% c("IGHA", "IGHG1","IGHG2B","IGHG2C","IGHG3"))# & sample_id == "+7d")

#mutation counts
# Calculate R and S mutation counts
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
# Show new mutation count columns
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)
# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)
# Calculate combined R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            combine=TRUE,
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)
#box plot
ggplot(db_obs, aes(x=c_call, y=mu_freq, fill=condition)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  stat_compare_means(aes(group = condition),label = "p.signif", method="t.test")+
  # scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot()
ggsave2("mutfreq.png",width=5, height=4,device="png")

#mutations within subregions
# Calculate R and S mutation counts for individual CDRs and FWRs
db_obs_v <- observedMutations(db, sequenceColumn="sequence_alignment",
                              germlineColumn="germline_alignment_d_mask",
                              regionDefinition=IMGT_V_BY_REGIONS,
                              frequency=FALSE, 
                              nproc=1)
# Show new FWR mutation columns
db_obs_v %>% 
  select(sequence_id, starts_with("mu_count_fwr")) %>%
  head(n=4)
# Calculate aggregate CDR and FWR V-segment R and S mutation frequencies
db_obs_v <- observedMutations(db_obs_v, sequenceColumn="sequence_alignment",
                              germlineColumn="germline_alignment_d_mask",
                              regionDefinition=IMGT_V,
                              frequency=TRUE, 
                              nproc=1)
# Show new CDR and FWR mutation frequency columns
db_obs_v %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)
#compare silent and replacement mutations
g2 <- ggplot(db_obs_v, aes(x=c_call, y=mu_freq_cdr_s, fill=condition)) +
  theme_bw() + ggtitle("CDR silent mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  stat_compare_means(aes(group = condition),label = "p.signif", method="t.test")+
  # scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot()
g3 <- ggplot(db_obs_v, aes(x=c_call, y=mu_freq_cdr_r, fill=condition)) +
  theme_bw() + ggtitle("CDR replacement mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  stat_compare_means(aes(group = condition),label = "p.signif", method="t.test")+
  # scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot()
png("silent.mut.freq.png",width=10, height=4,units="in",res=200)
wrap_plots(g2, g3, guides="collect")
dev.off()

#define mutations by physiochemical properties
# Calculate charge mutation frequency for the full sequence
db_obs_ch <- observedMutations(db, sequenceColumn="sequence_alignment",
                               germlineColumn="germline_alignment_d_mask",
                               regionDefinition=NULL,
                               mutationDefinition=CHARGE_MUTATIONS,
                               frequency=TRUE, 
                               nproc=1)
# Show new charge mutation frequency columns
db_obs_ch %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)
ggplot(db_obs_ch, aes(x=c_call, y=mu_freq_seq_r, fill=condition)) +
  theme_bw() + ggtitle("Charge replacement mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  stat_compare_means(aes(group = condition),label = "p.signif", method="t.test")+
  geom_boxplot()
ggsave2("charge.mutfreq.png",width=5, height=4,device="png")


###Selection pressure
# Collapse clonal groups into single sequences
clones3 <- collapseClones(ExampleDb, cloneColumn="ms_clone", 
                         sequenceColumn="sequence_alignment", 
                         germlineColumn="germline_alignment_d_mask", 
                         regionDefinition=IMGT_V, 
                         method="thresholdedFreq", minimumFrequency=0.6,
                         includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                         nproc=1)
# Count observed mutations and append mu_count columns to the output
observed <- observedMutations(clones3, 
                              sequenceColumn="clonal_sequence",
                              germlineColumn="clonal_germline",
                              regionDefinition=IMGT_V, nproc=1)
# Count expected mutations and append mu_exptected columns to the output
expected <- expectedMutations(observed, 
                              sequenceColumn="clonal_sequence",
                              germlineColumn="clonal_germline",
                              targetingModel=HH_S5F,
                              regionDefinition=IMGT_V, nproc=1)

# Calculate selection scores using the output from expectedMutations
baseline <- calcBaseline(expected, testStatistic="focused", 
                         regionDefinition=IMGT_V, nproc=1)

# Calculate selection on charge class with the mouse 5-mer model
baseline2 <- calcBaseline(clones3, testStatistic="focused", 
                         regionDefinition=IMGT_V, 
                         targetingModel=MK_RS5NF,
                         mutationDefinition=CHARGE_MUTATIONS,
                         nproc=1)

#convolve individual distributions
# Combine selection scores by time-point
grouped_1 <- groupBaseline(baseline, groupBy="condition")
# Subset the original data to switched isotypes
db_sub <- subset(ExampleDb, c_call %in% c("IGHM", "IGHG1","IGHG2B","IGHG2C","IGHG3"))
# Collapse clonal groups into single sequence
clones_sub <- collapseClones(db_sub, cloneColumn="ms_clone",
                             sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=IMGT_V, 
                             method="thresholdedFreq", minimumFrequency=0.6,
                             includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                             nproc=1)
# Calculate selection scores from scratch
baseline_sub <- calcBaseline(clones_sub, testStatistic="focused", 
                             regionDefinition=IMGT_V, nproc=1)
# Combine selection scores by time-point and isotype
grouped_2 <- groupBaseline(baseline_sub, groupBy=c("condition", "c_call"))

#convolve at multiple levels
# First group by subject and status
subject_grouped <- groupBaseline(baseline, groupBy=c("condition", "mouse_ID"))
# Then group the output by status
status_grouped <- groupBaseline(subject_grouped, groupBy="condition")

#comparing probability distributions
testBaseline(grouped_1, groupBy="condition")
# # Set sample and isotype colors
# sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
# isotype_colors <- c("IGHM"="darkorchid", "IGHD"="firebrick", 
#                     "IGHG"="seagreen", "IGHA"="steelblue")
# Plot mean and confidence interval by time-point
plotBaselineSummary(grouped_1, "condition")
ggsave2("selection.score.png",width=2, height=4,device="png")
plotBaselineSummary(status_grouped, "condition")
# Plot selection scores by time-point and isotype for only CDR
plotBaselineSummary(grouped_2, "condition", "c_call", subsetRegions="cdr")#groupColors=isotype_colors,
ggsave2("selection.score.isotype.png",width=5, height=4,device="png")
# Group by CDR/FWR and facet by isotype
plotBaselineSummary(grouped_2, "condition", "c_call", facetBy="group")
ggsave2("selection.score.isotype.fwr.png",width=4, height=9,device="png")
# Plot selection PDFs for a subset of the data
plotBaselineDensity(grouped_2, "c_call", groupColumn="condition", colorElement="group", sigmaLimits=c(-1, 1))#colorValues=sample_colors, 
ggsave2("selection.pdf.isotype.png",width=9, height=7,device="png")


#tuning clonal assignment
# Use nucleotide Hamming distance and normalize by junction length
dist_ham <- distToNearest(ExampleDb, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",#"v_call_genotyped"
                          model="ham", normalize="len", nproc=1)
# Use genotyped V assignments, a 5-mer model and no normalization
dist_aa <- distToNearest(ExampleDb, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",#"v_call_genotyped"
                          model="aa", normalize="none", nproc=1)
# Single-cell mode 
# Group cells in a one-stage process (VJthenLen=FALSE) and using
# both heavy and light chain sequences (onlyHeavy=FALSE)
dist_sc <- distToNearest(db, cellIdColumn="cell", locusColumn="locus", 
                         VJthenLen=FALSE, onlyHeavy=FALSE)

#threshold determination
#manual
# Generate Hamming distance histogram
ggplot(subset(dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
# Generate HH_S5F distance histogram
ggplot(subset(dist_aa, !is.na(dist_nearest)),aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("aa distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 50, 5)) +
  geom_histogram(color="white", binwidth=1) +
  geom_vline(xintercept=7, color="firebrick", linetype=2)

#automatic
# Find threshold using density method
output <- findThreshold(dist_ham$dist_nearest, method="density")
threshold <- output@threshold
# Plot distance histogram, density estimate and optimum threshold
png("hamming.png")
plot(output, title="Density Method")
dev.off()
print(output)

#automatic mixture model
# Find threshold using gmm method
output <- findThreshold(dist_ham$dist_nearest, method="gmm", model="gamma-gamma")
# Plot distance histogram, Gaussian fits, and optimum threshold
png("gmm.png")
plot(output, binwidth=0.02, title="GMM Method: gamma-gamma")
dev.off()
print(output)

#nearest neighor independent for subsets
dist_fields <- distToNearest(ExampleDb, model="ham", normalize="len", 
                             fields="condition", nproc=1)
ggplot(subset(dist_fields, !is.na(dist_nearest)), 
       aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Grouped Hamming distance") + 
  ylab("Count") +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2) +
  facet_grid(condition ~ ., scales="free_y")
ggsave2("grouped.hamming.png",width=5, height=4,device="png")

#nearest neighbor across groups
dist_cross <- distToNearest(ExampleDb, sequenceColumn="junction", 
                            vCallColumn="v_call", jCallColumn="j_call",#_genotyped
                            model="ham", first=FALSE, 
                            normalize="len", cross="condition", nproc=1)
ggplot(subset(dist_cross, !is.na(cross_dist_nearest)), 
       aes(x=cross_dist_nearest)) + 
  theme_bw() + 
  xlab("Cross-sample Hamming distance") + 
  ylab("Count") +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2) +
  facet_grid(condition ~ ., scales="free_y")
ggsave2("cross.hamming.png",width=5, height=4,device="png")


for(i in unique(ExampleDb$condition)){
  db2 <- subset(ExampleDb, condition==i)
  # Collapse sequences into clonal consensus
  clone_db <- collapseClones(db2, cloneColumn="ms_clone", 
                             sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             nproc=1)
  # Create targeting model in one step using only silent mutations
  # Use consensus sequence input and germline columns
  model <- createTargetingModel(clone_db, model="s", sequenceColumn="clonal_sequence", 
                                germlineColumn="clonal_germline", vCallColumn="v_call")
  # Generate hedgehog plot of mutability model
  png(paste0(i,".A.hedgehog.png"))
  plotMutability(model, nucleotides="A", style="hedgehog")
  dev.off()
  png(paste0(i,".C.hedgehog.png"))
  plotMutability(model, nucleotides="C", style="hedgehog")
  dev.off()
}







