library(alakazam)
library(ape)
library(dplyr)

# read in the data
db <- readIgphyml("filtered_contig_heavy_germ-pass_igphyml-pass.tab", format="phylo",
                  branches="mutations")

png("graph.png",width=8,height=6,unit="in",res=300)
plot(db$trees[[1]],show.node.label=TRUE)
add.scale.bar(length=5)
dev.off()