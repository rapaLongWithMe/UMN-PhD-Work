##Getting started with GENIE3
  #Trying to follow the Markdown from Katie's circadian drought paper

# Set up the environment
require(tidyverse)
require(GENIE3)

setwd("~/Desktop/Git/RNAseq/Brapa/GRNs/2022_YellowSarsons")
inputPath = "~/Desktop/Git/RNAseq/Brapa/GRNs/2022_YellowSarsons/R_input"
outputPath= "~/Desktop/Git/RNAseq/Brapa/GRNs/2022_YellowSarsons/R_output"
source(file.path(inputPath,"Rfunctions.R"))

# Define the set of transcription factors (TFs) for each of the networks.
  # Selected TFs using the MapMan;MERCATOR (https://mapman.gabipd.org/app/mercator), and taking all genes with the annotation "Transcriptional Control: 27.3"
BrTFs<-BrTFs$Accession

# Load clock genes; added as potential regulators
clock<-read.csv(file.path(inputPath,"Brapa_At_ClockIDs.txt"),stringsAsFactors = F)
BrClock<-unique(clock$BRA)
BrTFs<-union(BrTFs,BrClock)

# Only consider the TFs that are expressed
  ##### CHANGE BASED ON GENOTYPE #####
load(file.path(inputPath, file = "R500_TCs_List.RData"))
##### ##### ##### ##### ##### ##### ##### ##### ##### #####

#make a matrix
exprsMat = do.call(cbind, TCs_List)
#clip off the extra genotype and treatment information we added
rownames(exprsMat) = str_split(rownames(exprsMat), '_', simplify = TRUE)[,3]
is.matrix(exprsMat)
#pull out just the expressed TFs
BrTFs<-intersect(BrTFs,rownames(exprsMat))


# Next, we construct our Gene Regulator Networks (GRNs) using the GENIE3 framework (https://github.com/aertslab/GENIE3).
  # First, define the TF expression and target matrices
BrTFsMat<-exprsMat[BrTFs,]
  #choose the appropriate columns based on the type of network you are interested in
colnames(exprsMat)
BrTargetMat = exprsMat[, 25:48]
colnames(BrTFsMat)
BrTFsMat = BrTFsMat[, 25:48]


# Permute Target matrices to give 10,000 permuted targets
  # set the dimensions to be 10000*sample# and 10000 rows
BrTargetMatPerm<-matrix(data = sample(as.numeric(BrTargetMat),240000,replace = F),nrow = 10000, dimnames = list(row=paste("Decoy",1:10000,sep="_"),col=colnames(BrTargetMat)))


# Build the weight matrices
# ****If the R 'parallel' package is installed, I recommend increasing the 'nThr' parameter below ****
  # rows are targets and columns are TFs 
BrWeights<-RS.Get.Weigth.Matrix(target.matrix = t(BrTargetMat), input.matrix = t(BrTFsMat),nThr = 1)
  #write out the matrix to a csv to make sure these make sense
write.csv(BrWeights, file = "r500_WarmBrWeights.csv")
  #permute the weights to have something to test against
BrWeightsPerm<-RS.Get.Weigth.Matrix(target.matrix = t(BrTargetMatPerm), input.matrix = t(BrTFsMat), nThr=1)

##### CHANGE BASED ON GENOTYPE #####
save(BrWeights,BrWeightsPerm,file=file.path(outputPath,"r500_WarmWeightMatrix.RData"))
##### ##### ##### ##### ##### #####

# Impose an FDR edge weight cutoff to determine significant edge scores
BrWeightsCut.05<-CutGRNMat(tMat = BrWeights, pMat = BrWeightsPerm, cutPv = 0.05)

# Organize target groups
BrTargetGroups.05<-apply(BrWeightsCut.05,2,function(x) rownames(BrWeightsCut.05)[which(x>0)])

#save the target list
##### CHANGE BASED ON GENOTYPE #####
save(BrTargetGroups.05, file = file.path(outputPath, "r500_WarmWeightsCut.05.RData"))
lapply(BrTargetGroups.05, function(x) write.table(data.frame(x), 'r500_WarmWeightsCut.05.csv', append= T, sep=',' ))
##### ##### ##### ##### ##### #####



