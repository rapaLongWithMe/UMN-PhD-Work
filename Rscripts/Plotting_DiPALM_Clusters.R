### General script for plotting DiPALM clusters
 
## Set up the environment
setwd("~/Desktop/Git/JGI_Cold")

require(DiPALM)
require(WGCNA)
require(tidyverse)


##### IF ALREADY HAVE A PLOTTING DF AND PATTERN CLUSTERS OBJECT, SKIP DOWN TO LINE 43 #####
## Bring in the SigkMEs and LimmaModskMEs from the DiPALM output
input = "DiPALM/L58"

load(file.path(input, file = "L58_LimmaModskMEs.RData"))
load(file.path(input, file = "L58_SigkMEs.RData"))

## Double check that the correct geno was brought in 
rownames(LimmaModskMEs)
  # Clip off the variant descriptors to match with the SigkMEs
rownames(LimmaModskMEs) = gsub("\\.v[0-9].[0-9]", "", rownames(LimmaModskMEs))
  # Rename to make script more fluid between accessions
SigkMEs = L58_kMEs

## Cluster contrasts
LimmaModskMEsSig<-LimmaModskMEs[names(SigkMEs),]
patternCor<-cor(t(LimmaModskMEsSig))
patternTree<-hclust(as.dist(1-patternCor),method = "complete")

patternClusters<-cutreeDynamic(dendro = patternTree, minClusterSize = 100, distM = 1-patternCor, deepSplit = 1)
names(patternClusters)<-patternTree$labels
table(patternClusters)

## Convert clusters into a list of gene groups
patternClusters<-tapply(X = names(patternClusters), INDEX = patternClusters, function(x) x)
## Save the patternClusters so you don't have to rerun each time
save(patternClusters, file = "DiPALM/L58/L58_patternClusters.RData")
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####



## Bring in the MasterDF for plotting
load("ExpressedGenes/PlottingDFs/20231120_Pcglu_ExprsDF_forPlotting.RData")
  # Load the patternClusters if you're starting here
load("DiPALM/PcGlu/Pcglu_patternClusters.RData")


## Pull out a cluster to plot
  # Subset the exprsDf to contain genes only in cluster 1
#clust1 = exprsDF[which(exprsDF$Gene %in% patternClusters[[1]]),]
#clust2 = exprsDF[which(exprsDF$Gene %in% patternClusters[[2]]),]
#clust3 = exprsDF[which(exprsDF$Gene %in% patternClusters[[3]]),]
#clust4 = exprsDF[which(exprsDF$Gene %in% patternClusters[[4]]),]
#clust7 = exprsDF[which(exprsDF$Gene %in% patternClusters[[7]]),]
#clust11 = exprsDF[which(exprsDF$Gene %in% patternClusters[[11]]),]
#clust15 = exprsDF[which(exprsDF$Gene %in% patternClusters[[15]]),]
#clust36 = exprsDF[which(exprsDF$Gene %in% patternClusters[[36]]),]
#clust23 = exprsDF[which(exprsDF$Gene %in% patternClusters[[23]]),]
#clust = rbind(clust54, clust42, clust47)
clust = exprsDF[which(exprsDF$Gene %in% patternClusters[[3]]),]

## Calculate the average zScore
avgZ = clust %>% group_by(Treatment, ZT) %>%
  mutate(AvgZ = mean(zScore), 
         sdZ = sd(zScore)) %>%
  ungroup()

# Quick plot
avgZ %>% ggplot(aes(x = ZT, y = AvgZ, color = Treatment, group = Treatment)) +
  geom_line() +
  scale_color_manual(values = c("blue", "red"))


## Plot pretty normalized expression
  # Count the number of genes in this cluster of interest
length(avgZ$Gene %>% unique())

pdf("DiPALM/CompressedClusters_L58.pdf", width = 5, height = 4)
avgZ %>% ggplot(aes(x = ZT, y = AvgZ, 
                      color = Treatment, fill = Treatment, group = Treatment)) + 
  geom_line() +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  xlab("ZT") +
  ylab("Average normalized expression") +
  theme(legend.position = "bottom") +
  ggtitle(paste0("L58: ", length(avgZ$Gene %>% unique()), " genes"))
dev.off()


## Write out the gene annotations to match with the ATG to look at GO enrichment of these genes
clustGenes = avgZ$Gene %>% unique()
write.csv(clustGenes, file = "DiPALM/CompiledClusters_Pcglugenes.csv", row.names = FALSE)


## Bring in the syntenic table with ATGs
syntGenes = read.csv("Synteny/20231120_Genespace_SyntenicSheet.csv")
  # Simplify the DF
clustGenes_ATGs = syntGenes %>% filter(genome == "A03" | genome == "Athaliana")
  # Paralogs can be parsed form the og column; each og is orthologous to Arabidopsis
test = clustGenes_ATGs %>% pivot_wider(names_from = og, values_from = c(id, genome))

write.csv(clustGenes, file = "DiPALM/CompiledClusters_AO3genes.csv", row.names = FALSE)


