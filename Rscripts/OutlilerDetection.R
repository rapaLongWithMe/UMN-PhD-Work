### CHECKING FOR OUTLIERS

## Set Up
require(tidyverse)
require(WGCNA)

setwd("~/Desktop/Git/RNAseq/Brapa/ACC28")
#Define some path variables
inPath = "~/Desktop/Git/RNAseq/Brapa/ACC28/R_input/"
outPath = "~/Desktop/Git/RNAseq/Brapa/ACC28/R_output"

#load the pre-processed data
load(file.path(outPath, "ACC28_TCs_List.RData"))


#clean the data
DF = data.frame(TCs_List)
dim(DF)
names(DF)
#fix those crazy column names
names(DF) = gsub(x = names(DF), pattern = "Drought_R..", replacement = "")
names(DF) = gsub(x = names(DF), pattern = "Watered_R..", replacement = "")
names(DF)


#Cluster samples
sampleTree= hclust(dist(t(DF)), method = "average")

#plot
par(cex = 0.6)
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Using hclust to detect ACC28 outliers; all data", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

#Just control
colnames(DF)
ControlDF = DF[, 33:48]
colnames(ControlDF)
sampleTree2 = hclust(dist(t(ControlDF)), method = "average")

# Plot
plot(sampleTree2, main = "Using hclust to detect ACC28 outliers; control data", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

#Just cold
colnames(DF)
ColdDF = DF[, 1:32]
colnames(ColdDF)
sampleTree3 = hclust(dist(t(ColdDF)), method = "average")

# Plot
plot(sampleTree3, main = "Using hclust to detect ACC28 outliers; cold data", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)


#save the DFs to use downstream
setwd(outPath)
### CHANGE TO MATCH GENOTYPE ###
save(DF, file = "Av_DF.RData")
save(DF2, file = "Av_DroughtDF.RData")
save(DF3, file = "Av_WateredDF.RData")

