## Calculating response scores for all genotypes

#Set up
require(tidyverse)
require(WGCNA)

inPath = "~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes/AveragedExpression"
outPath = "~/Desktop/Git/RNAseq/Bnapus/SharedResponse/"

#load expressionAvg's from FinalClustering steps for all genotypes you want to re-cluster
#since will all show up in the environment as the same name, be sure to rename immediately before pulling in another genotype

### CHANGE BASED ON GENOTYPE ### 
load(file.path(inPath, "Av_AveragedExpression.RData"))
expressionAvg_G1 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G1)

load(file.path(inPath, "St_AveragedExpression.RData"))
expressionAvg_G2 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G2)

load(file.path(inPath, "Se_AveragedExpression.RData"))
expressionAvg_G3 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G3)

load(file.path(inPath, "Al_AveragedExpression.RData"))
expressionAvg_G4 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G4)

load(file.path(inPath, "Br_AveragedExpression.RData"))
expressionAvg_G5 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G5)

load(file.path(inPath, "Ca_AveragedExpression.RData"))
expressionAvg_G6 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G6)

load(file.path(inPath, "Ab_AveragedExpression.RData"))
expressionAvg_G7 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G7)

load(file.path(inPath, "Gr_AveragedExpression.RData"))
expressionAvg_G8 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G8)

load(file.path(inPath, "Mu_AveragedExpression.RData"))
expressionAvg_G9 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G9)

load(file.path(inPath, "Ne_AveragedExpression.RData"))
expressionAvg_G10 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G10)

load(file.path(inPath, "Yu_AveragedExpression.RData"))
expressionAvg_G11 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G11)

load(file.path(inPath, "DH20_AveragedExpression.RData"))
expressionAvg_G12 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G12)

load(file.path(inPath, "DH12_AveragedExpression.RData"))
expressionAvg_G13 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G13)

load(file.path(inPath, "Da_AveragedExpression.RData"))
expressionAvg_G14 = expressionAvg
rm(expressionAvg)
dimnames(expressionAvg_G14)
##skipping Ze and Qu, the permutation plots look weird!!

##Calculate the difference between Treatment and control for each genotype
  #First set up empty lists to iteratively fill
Diff1 = list()
Diff2 = list()
Diff3 = list()
Diff4 = list()
Diff5 = list()
Diff6 = list()
Diff7 = list()
Diff8 = list()
Diff9 = list()
Diff10 = list()
Diff11 = list()
Diff12 = list()
Diff13 = list()
Diff14 = list()

##### MAKE SURE EACH Diff OBJECT MATCHES THE CORRECT GENOTYPE #####
colnames(expressionAvg_G1)
Diff1$Av_Zt1 = (expressionAvg_G1[, 5] - expressionAvg_G1[, 1])
Diff1$Av_Zt7 = (expressionAvg_G1[, 6] - expressionAvg_G1[, 2])
Diff1$Av_Zt13 = (expressionAvg_G1[, 7] - expressionAvg_G1[, 3])
Diff1$Av_Zt19 = (expressionAvg_G1[, 8] - expressionAvg_G1[, 4])

Diff2$St_Zt1 = (expressionAvg_G2[, 5] - expressionAvg_G2[, 1])
Diff2$St_Zt7 = (expressionAvg_G2[, 6] - expressionAvg_G2[, 2])
Diff2$St_Zt13 = (expressionAvg_G2[, 7] - expressionAvg_G2[, 3])
Diff2$St_Zt19 = (expressionAvg_G2[, 8] - expressionAvg_G2[, 4])

Diff3$Se_Zt1 = (expressionAvg_G3[, 5] - expressionAvg_G3[, 1])
Diff3$Se_Zt7 = (expressionAvg_G3[, 6] - expressionAvg_G3[, 1])
Diff3$Se_Zt13 = (expressionAvg_G3[, 7] - expressionAvg_G3[, 3])
Diff3$Se_Zt19 = (expressionAvg_G3[, 8] - expressionAvg_G3[, 4])

Diff4$Al_Zt1 = (expressionAvg_G4[, 5] - expressionAvg_G4[, 1])
Diff4$Al_Zt7 = (expressionAvg_G4[, 6] - expressionAvg_G4[, 2])
Diff4$Al_Zt13 = (expressionAvg_G4[, 7] - expressionAvg_G4[, 3])
Diff4$Al_Zt19 = (expressionAvg_G4[, 8] - expressionAvg_G4[, 4])

Diff5$Br_Zt1 = (expressionAvg_G5[, 5] - expressionAvg_G5[, 1])
Diff5$Br_Zt7 = (expressionAvg_G5[, 6] - expressionAvg_G5[, 2])
Diff5$Br_Zt13 = (expressionAvg_G5[, 7] - expressionAvg_G5[, 3])
Diff5$Br_Zt19 = (expressionAvg_G5[, 8] - expressionAvg_G5[, 4])

Diff6$Ca_Zt1 = (expressionAvg_G6[, 5] - expressionAvg_G6[, 1])
Diff6$Ca_Zt7 = (expressionAvg_G6[, 6] - expressionAvg_G6[, 2])
Diff6$Ca_Zt13 = (expressionAvg_G6[, 7] - expressionAvg_G6[, 3])
Diff6$Ca_Zt19 = (expressionAvg_G6[, 8] - expressionAvg_G6[, 4])

Diff7$Ab_Zt1 = (expressionAvg_G7[, 5] - expressionAvg_G7[, 1])
Diff7$Ab_Zt7 = (expressionAvg_G7[, 6] - expressionAvg_G7[, 2])
Diff7$Ab_Zt13 = (expressionAvg_G7[, 7] - expressionAvg_G7[, 3])
Diff7$Ab_Zt19 = (expressionAvg_G7[, 8] - expressionAvg_G7[, 4])

Diff8$Gr_Zt1 = (expressionAvg_G8[, 5] - expressionAvg_G8[, 1])
Diff8$Gr_Zt7 = (expressionAvg_G8[, 6] - expressionAvg_G8[, 2])
Diff8$Gr_Zt13 = (expressionAvg_G8[, 7] - expressionAvg_G8[, 3])
Diff8$Gr_Zt19 = (expressionAvg_G8[, 8] - expressionAvg_G8[, 4])

Diff9$Mu_Zt1 = (expressionAvg_G9[, 5] - expressionAvg_G9[, 1])
Diff9$Mu_Zt7 = (expressionAvg_G9[, 6] - expressionAvg_G9[, 2])
Diff9$Mu_Zt13 = (expressionAvg_G9[, 7] - expressionAvg_G9[, 3])
Diff9$Mu_Zt19 = (expressionAvg_G9[, 8] - expressionAvg_G9[, 4])

Diff10$Ne_Zt1 = (expressionAvg_G10[, 5] - expressionAvg_G10[, 1])
Diff10$Ne_Zt7 = (expressionAvg_G10[, 6] - expressionAvg_G10[, 2])
Diff10$Ne_Zt13 = (expressionAvg_G10[, 7] - expressionAvg_G10[, 3])
Diff10$Ne_Zt19 = (expressionAvg_G10[, 8] - expressionAvg_G10[, 4])

Diff11$Yu_Zt1 = (expressionAvg_G11[, 5] - expressionAvg_G11[, 1])
Diff11$Yu_Zt7 = (expressionAvg_G11[, 6] - expressionAvg_G11[, 2])
Diff11$Yu_Zt13 = (expressionAvg_G11[, 7] - expressionAvg_G11[, 3])
Diff11$Yu_Zt19 = (expressionAvg_G11[, 8] - expressionAvg_G11[, 4])

Diff12$DH20_Zt1 = (expressionAvg_G12[, 5] - expressionAvg_G12[, 1])
Diff12$DH20_Zt7 = (expressionAvg_G12[, 6] - expressionAvg_G12[, 2])
Diff12$DH20_Zt13 = (expressionAvg_G12[, 7] - expressionAvg_G12[, 3])
Diff12$DH20_Zt19 = (expressionAvg_G12[, 8] - expressionAvg_G12[, 4])

Diff13$DH12_Zt1 = (expressionAvg_G13[, 5] - expressionAvg_G13[, 1])
Diff13$DH12_Zt7 = (expressionAvg_G13[, 6] - expressionAvg_G13[, 2])
Diff13$DH12_Zt13 = (expressionAvg_G13[, 7] - expressionAvg_G13[, 3])
Diff13$DH12_Zt19 = (expressionAvg_G13[, 8] - expressionAvg_G13[, 4])

Diff14$Da_Zt1 = (expressionAvg_G14[, 5] - expressionAvg_G14[, 1])
Diff14$Da_Zt7 = (expressionAvg_G14[, 6] - expressionAvg_G14[, 2])
Diff14$Da_Zt13 = (expressionAvg_G14[, 7] - expressionAvg_G14[, 3])
Diff14$Da_Zt19 = (expressionAvg_G14[, 8] - expressionAvg_G14[, 4])

#pull out expressed genes that are shared between the genotypes
shared_genes1 = expressionAvg_G1[which(rownames(expressionAvg_G1) %in% rownames(expressionAvg_G2)), ]
shared_genes2= shared_genes1[which(rownames(shared_genes1) %in% rownames(expressionAvg_G3)), ]
shared_genes3= shared_genes2[which(rownames(shared_genes2) %in% rownames(expressionAvg_G4)), ]
shared_genes4= shared_genes3[which(rownames(shared_genes3) %in% rownames(expressionAvg_G5)), ]
shared_genes5= shared_genes4[which(rownames(shared_genes4) %in% rownames(expressionAvg_G6)), ]
shared_genes6= shared_genes5[which(rownames(shared_genes5) %in% rownames(expressionAvg_G7)), ]
shared_genes7= shared_genes6[which(rownames(shared_genes6) %in% rownames(expressionAvg_G8)), ]
shared_genes8= shared_genes7[which(rownames(shared_genes7) %in% rownames(expressionAvg_G9)), ]
shared_genes9= shared_genes8[which(rownames(shared_genes8) %in% rownames(expressionAvg_G10)), ]
shared_genes10= shared_genes9[which(rownames(shared_genes9) %in% rownames(expressionAvg_G11)), ]
shared_genes11= shared_genes10[which(rownames(shared_genes10) %in% rownames(expressionAvg_G12)), ]
shared_genes12= shared_genes11[which(rownames(shared_genes11) %in% rownames(expressionAvg_G13)), ]

shared_genesAll= shared_genes12[which(rownames(shared_genes12) %in% rownames(expressionAvg_G14)), ]
#make sure to retain only genes that start with Bra and contain Bo in the name
  ##NOTE: this removes ENSRA genes, which may be things like small RNAs that we are not necessarily interested in 
shared_genesAll = shared_genesAll[str_detect(rownames(shared_genesAll), "^Bra|Bo"),]
tail(shared_genesAll)
#create a vector of gene names that are shared between the genotypes
shared_gene_names = rownames(shared_genesAll)
tail(shared_gene_names)

#save to use later
getwd()
write.csv(shared_gene_names, file.path(outPath, file = "shared_gene_names.csv"), row.names = FALSE)

#subset each expression matrix to only retain these shared genes
sharedGenes_G1 = list()

for(i in 1:length(Diff1)){
  sharedGenes_G1[[i]] = Diff1[[i]][shared_gene_names]
  #name the elements of each list
  names(sharedGenes_G1)[[i]] = names(Diff1)[[i]]
  #add genotype descriptor to the gene names
  names(sharedGenes_G1[[i]]) = paste0("Av_", names(sharedGenes_G1[[i]]))
}

sharedGenes_G2 = list()

for(i in 1:length(Diff2)){
  sharedGenes_G2[[i]] = Diff2[[i]][shared_gene_names]
  names(sharedGenes_G2)[[i]] = names(Diff2)[[i]]
  names(sharedGenes_G2[[i]]) = paste0("St_", names(sharedGenes_G2[[i]]))
}

sharedGenes_G3 = list()

for(i in 1:length(Diff3)){
  sharedGenes_G3[[i]] = Diff3[[i]][shared_gene_names]
  names(sharedGenes_G3)[[i]] = names(Diff3)[[i]]
  names(sharedGenes_G3[[i]]) = paste0("Se_", names(sharedGenes_G3[[i]]))
}

sharedGenes_G4 = list()

for(i in 1:length(Diff4)){
  sharedGenes_G4[[i]] = Diff4[[i]][shared_gene_names]
  names(sharedGenes_G4)[[i]] = names(Diff4)[[i]]
  names(sharedGenes_G4[[i]]) = paste0("Al_", names(sharedGenes_G4[[i]]))
}

sharedGenes_G5 = list()

for(i in 1:length(Diff5)){
  sharedGenes_G5[[i]] = Diff5[[i]][shared_gene_names]
  names(sharedGenes_G5)[[i]] = names(Diff5)[[i]]
  names(sharedGenes_G5[[i]]) = paste0("Br_", names(sharedGenes_G5[[i]]))
}

sharedGenes_G6 = list()

for(i in 1:length(Diff6)){
  sharedGenes_G6[[i]] = Diff6[[i]][shared_gene_names]
  names(sharedGenes_G6)[[i]] = names(Diff6)[[i]]
  names(sharedGenes_G6[[i]]) = paste0("Ca_", names(sharedGenes_G6[[i]]))
}

sharedGenes_G7 = list()

for(i in 1:length(Diff7)){
  sharedGenes_G7[[i]] = Diff7[[i]][shared_gene_names]
  names(sharedGenes_G7)[[i]] = names(Diff7)[[i]]
  names(sharedGenes_G7[[i]]) = paste0("Ab_", names(sharedGenes_G7[[i]]))
}

sharedGenes_G8 = list()

for(i in 1:length(Diff8)){
  sharedGenes_G8[[i]] = Diff8[[i]][shared_gene_names]
  names(sharedGenes_G8)[[i]] = names(Diff8)[[i]]
  names(sharedGenes_G8[[i]]) = paste0("Gr_", names(sharedGenes_G8[[i]]))
}

sharedGenes_G9 = list()

for(i in 1:length(Diff9)){
  sharedGenes_G9[[i]] = Diff9[[i]][shared_gene_names]
  names(sharedGenes_G9)[[i]] = names(Diff9)[[i]]
  names(sharedGenes_G9[[i]]) = paste0("Mu_", names(sharedGenes_G9[[i]]))
}

sharedGenes_G10 = list()

for(i in 1:length(Diff10)){
  sharedGenes_G10[[i]] = Diff10[[i]][shared_gene_names]
  names(sharedGenes_G10)[[i]] = names(Diff10)[[i]]
  names(sharedGenes_G10[[i]]) = paste0("Ne_", names(sharedGenes_G10[[i]]))
}

sharedGenes_G11 = list()

for(i in 1:length(Diff11)){
  sharedGenes_G11[[i]] = Diff11[[i]][shared_gene_names]
  names(sharedGenes_G11)[[i]] = names(Diff11)[[i]]
  names(sharedGenes_G11[[i]]) = paste0("Yu_", names(sharedGenes_G11[[i]]))
}

sharedGenes_G12 = list()

for(i in 1:length(Diff12)){
  sharedGenes_G12[[i]] = Diff12[[i]][shared_gene_names]
  names(sharedGenes_G12)[[i]] = names(Diff12)[[i]]
  names(sharedGenes_G12[[i]]) = paste0("DH20_", names(sharedGenes_G12[[i]]))
}

sharedGenes_G13 = list()

for(i in 1:length(Diff13)){
  sharedGenes_G13[[i]] = Diff13[[i]][shared_gene_names]
  names(sharedGenes_G13)[[i]] = names(Diff13)[[i]]
  names(sharedGenes_G13[[i]]) = paste0("DH12_", names(sharedGenes_G13[[i]]))
}

sharedGenes_G14 = list()

for(i in 1:length(Diff14)){
  sharedGenes_G14[[i]] = Diff14[[i]][shared_gene_names]
  names(sharedGenes_G14)[[i]] = names(Diff14)[[i]]
  names(sharedGenes_G14[[i]]) = paste0("Da_", names(sharedGenes_G14[[i]]))
}

#bind list elements together into a matrix so can run through WGCNA later
G1 = do.call(cbind, sharedGenes_G1)
G2 = do.call(cbind, sharedGenes_G2)
G3 = do.call(cbind, sharedGenes_G3)
G4 = do.call(cbind, sharedGenes_G4)
G5 = do.call(cbind, sharedGenes_G5)
G6 = do.call(cbind, sharedGenes_G6)
G7 = do.call(cbind, sharedGenes_G7)
G8 = do.call(cbind, sharedGenes_G8)
G9 = do.call(cbind, sharedGenes_G9)
G10 = do.call(cbind, sharedGenes_G10)
G11 = do.call(cbind, sharedGenes_G11)
G12 = do.call(cbind, sharedGenes_G12)
G13 = do.call(cbind, sharedGenes_G13)
G14 = do.call(cbind, sharedGenes_G14)

#bind the matrices together
TCs_Mat_shared = rbind(G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12, G13, G14)
head(TCs_Mat_shared)
#fix the column names
colnames(TCs_Mat_shared) = c("Zt1", "Zt7", "Zt13", "Zt19")
tail(TCs_Mat_shared) #Da was the last genotype to add so should be in the tail

#WGCNA cannot take negative values
#standardize the response scores such that the lowest value is greater than zero 
minVal = abs(min(TCs_Mat_shared)) + 1
Standard_TCs_Mat_shared = minVal + TCs_Mat_shared
#double check
head(Standard_TCs_Mat_shared)
head(TCs_Mat_shared)
range(Standard_TCs_Mat_shared)
range(TCs_Mat_shared)
min(Standard_TCs_Mat_shared)
min(TCs_Mat_shared)
which(Standard_TCs_Mat_shared<= 0)
which(TCs_Mat_shared<= 0)

#make sure there aren't any NAs
which(is.na(TCs_Mat_shared)==TRUE)


#save
setwd(outPath)
save(TCs_Mat_shared, file = "AllButZeQu_TCs_Mat_March3.RData")
save(Standard_TCs_Mat_shared, file = "AllButZeQu_Standard_TCs_Mat_March3.RData")










