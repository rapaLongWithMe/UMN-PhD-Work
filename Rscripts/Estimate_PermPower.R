#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=4
#SBATCH --mem=32g
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ricon001@umn.edu

module load R/4.3.0-openblas-rocky8
R

### MSI script to estimate SoftThresholdingPower across multiple iterations
require(WGCNA)
require(tidyverse)

## Load the expression matrix
load("ACC28MappedToR500_DF.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2))
cex1 = 0.9

#Rows correspond to samples and columns to genes
threshold = pickSoftThreshold(t(DF), powerVector = powers, verbose = 5)


#Plot to select parameters
pdf("ACC28_PowerEstimate1.pdf")
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
dev.off()