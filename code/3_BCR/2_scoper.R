# Imports
library(scoper)
library(dplyr)


ExampleDb<-read.table(file="filtered_contig_heavy_productive-T.tsv",header=T,sep="\t",stringsAsFactors = F)
  
# Clonal assignment using chosen method
results <- hierarchicalClones(ExampleDb, threshold=0.15)
# Get results data.frame
results_db <- as.data.frame(results)
glimpse(results_db)
# Plot a histogram of inter clonal distances
png("plot.png")
plot(results, binwidth=0.02)
dev.off()
# Get summary data.frame
glimpse(summary(results))


write.table(results_db,file="filtered_contig_heavy_clone-pass.tsv",col.names=T,sep="\t")
