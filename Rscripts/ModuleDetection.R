### MODULE DETECTION

## Set Up
require(tidyverse)
require(WGCNA)

### CHANGE BASED ON GENOTYPE ###
setwd("~/Desktop/Git/RNAseq/Brapa/ACC28")
outPath = "~/Desktop/Git/RNAseq/Brapa/ACC28"

load(file.path(outPath, "ACC28_AvgTCs_Mat.RData"))
colnames(TCs_Mat)
#the next line is only if the network has already been constructed
#load(file.path(outPath, "Da_BlockModsAll18_Jan4.RData"))
################################


#Estimate parameters
#Allow parallel processing for multiple cores
allowWGCNAThreads(8)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2))
cex1 = 0.9

#Rows correspond to samples and columns to genes
threshold = pickSoftThreshold(t(DF), powerVector = powers, verbose = 5)


#Plot to select parameters
# Scale free topology fit
plot(threshold$fitIndices[,1],
     -sign(threshold$fitIndices[,3])*threshold$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(threshold$fitIndices[,1],
     -sign(threshold$fitIndices[,3])*threshold$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power 
plot(threshold$fitIndices[,1], threshold$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(threshold$fitIndices[,1], threshold$fitIndices[,5], labels=powers, cex=cex1,col="red")

#construct network
### BLOCKING THIS OUT FOR NOW BECAUSE PROBABLY HAS TO BE RUN IN MSI
#BlockModsAll <- blockwiseModules(datExpr = t(TCs_Mat), power = 18, networkType = "signed", corType="bicor", TOMType="signed", minModuleSize=100, mergeCutHeight=0.30, deepSplit=1, pamRespectsDendro = F, nThreads = 4, verbose=3)

#save to avoid re-running in the future
## CHANGE FILE NAME BASED ON GENOTYPE ###
#save(BlockModsAll, file="Av_BlockModsAll.RData")

#Bring in the Modules if ran in MSI
load(file.path(outPath, "Av_BlockModsAll.RData"))

#Convert labels to colors for plotting
##the $colors dim has annotations as rownames; aka which genes are in which module
mergedColors = labels2colors(BlockModsAll$colors)
plotDendroAndColors(BlockModsAll$dendrograms[[1]], mergedColors[BlockModsAll$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                    main = "Module Relationships")

#identify how many modules were detected
dim(BlockModsAll[[3]])
MEs = (BlockModsAll[[3]])

#plot the module eigengenes to inspect patterns
#Add Zts to plot
MEs_Zts = rownames_to_column(MEs, var = "Zts")
colnames(MEs_Zts)
MEs_Zts$Zts

###Remove everything but the time point number and convert to zts###
MEs_Zts$Zts = as.numeric(c(17, 21, 1, 5, 9, 13))
#################################
MEs_Zts$Zts

#lengthen
plottingMEs = pivot_longer(MEs_Zts, cols = !Zts, names_to = "Modules", values_to = "Expression")

#plot
plottingMEs %>%
        group_by(Modules) %>%
        ggplot(aes(x = as.numeric(Zts), 
                   y = Expression, 
                   fill = "black")) +
        geom_line() +
        facet_wrap(~Modules) +
        xlab("Zt") +
        ylab("Time series response") +
        ggtitle("Co-expression Eigengene Patterns") +
  theme_bw()

