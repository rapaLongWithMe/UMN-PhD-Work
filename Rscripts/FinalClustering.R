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
load(file.path(outPath, "DH12_TCs_List.RData"))
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

lapply(TreatByZt, names)
dim(TreatByZt[[1]])
#want to average the 3 cols of each element into a single col
expressionAvg = list()
tc_means = list()
for (i in 1:length(TreatByZt)) {
  tc_means[[i]] = rowMeans(TreatByZt[[i]])
  expressionAvg = do.call(cbind, tc_means)
}

#pull colnames 
eMatCols<-(unique(tnms))
#this seems to work but looks weird in the environment...fine in View
colnames(expressionAvg) = eMatCols
#save averaged expression
setwd(outPath)
getwd()

### CHANGE BASED ON GENOTYPE ###
save(expressionAvg, file = "DH12_AveragedExpression.RData")
################################

#make a color palette
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

### CHANGE BASED ON GENOTYPE ### 
pdf(file.path(FigPath,"DH12_Ribbons_Jan18.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ### 
lapply(rapply(TCs_List, function(x) gsub("Cold_", "", rownames(x)), how = "list"), matrix)
TCs_List2 = lapply(TCs_List2, function(x) gsub("Warm_", "", rownames(x)))


TCs_List2 = TCs_List
for (i in 1:length(TCs_List)) {
  rownames(TCs_List2[[i]]) = gsub("Cold_", "", rownames(TCs_List2[[i]])) 
  rownames(TCs_List2[[i]]) = gsub("Warm_", "", rownames(TCs_List2[[i]])) 
}

patternClusters2 = patternClusters
for (i in 1:length(patternClusters2)) {
patternClusters2[[i]] = gsub("Cold_", "", patternClusters2[[i]])
patternClusters2[[i]] = gsub("Warm_", "", patternClusters2[[i]])
}


for(i in 1:length(patternClusters2)) {
  PlotTCsRibbon(TClst = TCs_List2, 
                tgenes = patternClusters2[[i]],
                scale = T, 
                tcols = c("blue","blue","blue","blue",
                          "red","red","red","red"), 
                tltys = c(1, 2, 3, 4, 1, 2, 3,4), 
                main = paste("Pattern Change Cluster:", names(patternClusters)[i],"#Genes:",length(patternClusters[[i]]),sep=" "))
}
dev.off()


## ORGANIZING BY PHASE
#find the timepoint of max expression per sample
#first pull cold/warm samples, then add a column containing the max which is then used to re-order the average expression. Once that's done remove the Peak_tp column
#can do either cold or warm, starting with drought 
colnames(expressionAvg)
treat_expressionAvg = data.frame(expressionAvg[,1:4])
treat_expressionAvg = treat_expressionAvg %>%
  mutate(Peak_tp =  max.col(.))

test = treat_expressionAvg %>%
  arrange(desc(.))
#check
levels(as.factor(treat_expressionAvg$Peak_tp))
#reorder
treat_expressionAvg = arrange(treat_expressionAvg, Peak_tp)%>%
  select(-Peak_tp)

#maybe this is better
#treat_expressionAvg = arrange(treat_expressionAvg, desc(across(contains("Drought"))))

#this will re-order the (matched) gene names from SigkMEs and order them by treat_expressionAvg, which we have already organized by phase
phaseOrder = treat_expressionAvg[order(match(rownames(treat_expressionAvg), rownames(SigkMEs))),]

phaseOrder = test[order(match(rownames(test), rownames(SigkMEs))),]


#doubleCheck
head(phaseOrder)
head(treat_expressionAvg)
#patternTree should be different; that was the original order 
head(patternTree$labels)

#now want to subset all expressed genes and keep only the significant ones
#check the number of genes
nrow(phaseOrder)
length(patternTree$labels)

phaseOrder_SigkMEs = phaseOrder[names(SigkMEs),]

phaseOrder_SigkMEs = phaseOrder_SigkMEs %>% arrange(desc(.))


#check the number of genes
nrow(phaseOrder_SigkMEs)
length(patternTree$labels)

## DOUBLE CHECK HERE, AM I SUBSETTING BY THE RIGHT THING???
#subset the phaseOrder vector to only include sig. genes
phaseOrder_SigkMEs = phaseOrder_SigkMEs[order(match(rownames(phaseOrder_SigkMEs), rownames(treat_expressionAvg))),]
#final check 
head(phaseOrder_SigkMEs)
head(treat_expressionAvg)
#patternTree should be different; that was the original order 
head(patternTree$labels)

#reset the patternTree labels 
patternTree$labels = rownames(phaseOrder_SigkMEs)
#re-plot the heatmaps which should now be ordered by phase
#pheatmap2 = pheatmap(mat = treat_expressionAvg[names(phaseOrder_SigkMEs),], cluster_rows = patternTree, cluster_cols = F,scale = "row", color = colFunc(25), gaps_col = 4, show_rownames = F, main = "Expression Patterns")

## Phase order and plot coexpression modules
#genes are already clustered into modules, extract, average and order
#this prints the time point of max avg expression for each cluster
colnames(BlockModsAll$MEs)
rownames(BlockModsAll$MEs)
PeakPhase_MEs = data.frame(t(BlockModsAll$MEs)) %>% 
  rownames_to_column('module') %>%  # creates an ID column
  
  ### CHANGE BASED ON GENOTYPE ### 
  gather(timepoint, expression, DH12_1_Drought_R1:DH12_19_Drought_R1) %>% 
  ###############################
  group_by(module) %>% 
  slice(which.max(expression))

PeakPhase_MEs = arrange(PeakPhase_MEs, desc(expression))
#make a vector that contains the phase-ordered colors (MEs)
PeakPhase_Colors = PeakPhase_MEs$module
#so want to use this to re-order the eigengenes
MEs = BlockModsAll$MEs
PeakPhase_Eigengene = MEs[, PeakPhase_Colors]
PeakPhase_Eigengene = arrange(PeakPhase_Eigengene, desc(across(starts_with("ME"))))
#check
colnames(PeakPhase_Eigengene)
colnames(MEs)
PeakPhase_Colors

### FIX COLOR PALETTE!!
#pdf(file.path(outputPath,"Figure2_C.pdf"),width = 6,height = 6)
pheatmap3 = pheatmap(t(PeakPhase_Eigengene),cluster_rows = F, cluster_cols = F, scale = "row", main = "Module expression")


