## SetUp
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("edgeR")

require(edgeR)
require(readr)
require(tidyverse)

##Bring in the data
  #compile all files at once into a single dataframe
setwd("~/Desktop/Git/RNAseq/Brapa/ACC28/RawCounts")
  #create some paths to save 
inPath = "~/Desktop/Git/RNAseq/Brapa/ACC28/R_input"

#create a vector that contains the name of files to loop over
filenames <- list.files(pattern="*.txt")
#make sure these names look correct and you have the right number
  #8tps*4reps*2treatments = 64
print(filenames)

#compile all files into a list with one element for each file
rapa = lapply(filenames, function(i) {
  read.table(i, header = TRUE) 
})
names(rapa) = filenames


#initialize some empty vectors to fill in
NormCounts = ""
cnts = ""
## ADD ROWNAMES FROM GENE ID ####
for (i in 1:length(rapa)) {
  print(i)
  rownames(rapa[[i]]) = rapa[[i]]$Geneid
  print(rownames(rapa[[i]]))
  geneLen = setNames(abs(rapa[[i]]$Length), nm = rapa[[i]]$Geneid)
  cnts[[i]] = rapa[[i]][7]
  cnts = do.call(cbind, cnts)
  cntsDge<- DGEList(counts = cnts)
  cntsDge<- calcNormFactors(cntsDge)
  NormCounts<- rpkm(cntsDge, log=T, gene.length = geneLen, prior.count=0)
}

#save the normalized counts to an appropriate location
  ##### CHANGE FILENAME BASED ON GENOTYPE #####
setwd(inPath)
save(NormCounts, file = "ACC28_NormCounts.RData") 
##### ##### ##### ##### ##### ##### ##### #####



