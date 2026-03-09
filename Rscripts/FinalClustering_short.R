### FINAL CLULSTERING

#Set Up
require(DiPALM)
require(tidyverse)
require(pheatmap)
require(WGCNA)

### CHANGE BASED ON GENOTYPE ###
setwd("~/Desktop/Git/RNAseq/R500")
#Define some path variables
inPath = "~/Desktop/Git/RNAseq/R500/R_input"
FigPath = "~/Desktop/Git/RNAseq/R500/Figures"
outPath = "~/Desktop/Git/RNAseq/Bnapus/Yu_Qu_DH20_DH12/R_output/DH12"

load(file.path(outPath, "DH12_LimmaModskMEs.RData"))
load(file.path(outPath, "DH12_SigkMEs.RData"))
load(file.path(outPath, "TCsR500.RData"))
load(file.path(outPath, "DH12_BlockModsAll18_Dec18.RData"))
load(file.path(outPath, "DH12_TestSumskMEs.RData"))
#################################


#cluster contrasts for differentially patterned genes into modules; will be used in heat map later
LimmaModskMEsSig<-LimmaModskMEs[names(SigkMEs),]
patternCor<-cor(t(LimmaModskMEsSig))
patternTree<-hclust(as.dist(1-patternCor),method = "complete")

#cluster contrasts for differential median expression genes into modules
#LimmaModsMedSig<-LimmaModsMed[names(SigMed),]
#patternCorMed<-cor(t(LimmaModsMedSig))
#patternTreeMed<-hclust(as.dist(1-patternCorMed),method = "complete")


#average the replicates to simplify the plot
#similar to TCsR500 but bound by columns (samples) not rows (genes)
expressionMat<-do.call(cbind,TCs_List)
colnames(expressionMat)

#rearrange so the colnames are "Treatment_Zt"
tmp<-expressionMat
spNms<-strsplit(x = colnames(expressionMat), split = "_")
tnms<-sapply(spNms,function(x) paste(x[c(3, 2)],collapse = "_"))
TreatByZt<-tapply(1:length(tnms),INDEX = tnms, function(x) tmp[,x])
names(TreatByZt)
  # reorder by zt
TreatByZt = TreatByZt[c(7,11, 12, 8, 9, 10, 1, 5, 6, 2, 3, 4)]


#want to average the 3 cols of each element into a single col
expressionAvg = list()
tc_means = list()
for (i in 1:length(TreatByZt)) {
  tc_means[[i]] = data.frame(AvgExprs = rowMeans(TreatByZt[[i]]))
  expressionAvg = do.call(cbind, tc_means)
}
colnames(expressionAvg) = paste("AvgExprs", names(TreatByZt), sep = "_")


### CHANGE BASED ON GENOTYPE ###
save(expressionAvg, file = file.path(outPath, "TSU_AvgExprs.RData"))

## Make a heatmap
  # set up a color palette
colFunc<-colorRampPalette(colors = c("darkblue","blue","lightblue","white","orange"))

#make sure the subsetting is happening correctly
nrow(expressionAvg[names(SigkMEs),])
#plot
pheatmap1_kME = pheatmap(mat = expressionAvg[names(SigkMEs),], cluster_rows = patternTree, cluster_cols = F,scale = "row", color = colFunc(25), gaps_col = 4, show_rownames = F, main = "Expression Patterns")
#pheatmap1_Med = pheatmap(mat = expressionAvg[names(SigMed),], cluster_rows = patternTree, cluster_cols = F,scale = "row", color = colFunc(25), gaps_col = 4, show_rownames = F, main = "Expression Patterns")


##Re-cluster genes based on similarity of expression pattern (kME)

### CHANGE BASED ON THE NUMBER OF DIFFERENTIAL PATTERN CHANGE GENES ###
  # FOR ANY SET OF GENES LESS <= 100 I AM REQUIRING ONLY 25 GENES PER CLUSTER
patternClusters<-cutreeDynamic(dendro =patternTree, minClusterSize = 100, distM = 1-patternCor, deepSplit = 1)
###########################################################################

#adds the gene names back as labels; now have a number for each gene that corresponds to what cluster it belongs in
names(patternClusters)<-patternTree$labels
#number of genes per cluster
table(patternClusters)

#put clusters into a list
patternClusters<-tapply(X = names(patternClusters), INDEX = patternClusters, function(x) x)

##### CHANGE BASED ON GENOTYPE #####
getwd()
save(patternClusters, file = "DH12_kME_PatternClusters_Jan18.RData")
##### ##### ##### ##### ##### #####

clustScores<-sapply(patternClusters,function(x) mean(TestSumskMEs[x]))


#print all clusters into a single pdf file
### CHANGE BASED ON GENOTYPE ### 
pdf(file.path(FigPath,"DH12_Ribbons_Jan18.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ### 
for(i in 1:length(patternClusters)) {
  PlotTCsRibbon(TClst = TCs_List, 
                tgenes = patternClusters[[i]],
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1", "darkgoldenrod1", "chartreuse4","chartreuse4","chartreuse4"), 
                tltys = c(1, 2, 3, 1, 2, 3), 
                main = paste("Pattern Change Cluster:", names(patternClusters)[i],"#Genes:",length(patternClusters[[i]]),sep=" "))
}
dev.off()


