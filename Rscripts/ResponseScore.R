### Response Score to cluster multiple genotypes
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

### CHANGE BASED ON GENOTYPE ###
# Calculating the difference of expression between Watered and Drought at each time point
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
##############################

#pull out expressed genes that are shared between the genotypes
  #found some close examples here: https://stackoverflow.com/questions/39084278/merge-multiple-tables-by-row-and-column-in-r

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
shared_genesAll = shared_genesAll[str_detect(rownames(shared_genesAll), "^Bra|Bo"),]
tail(shared_genesAll)

#create a vector of gene names that are shared between the genotypes
shared_gene_names = rownames(shared_genesAll)
tail(shared_gene_names)

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
  ### NOTE; THIS MAKES COLUMN NAMES ARBITRARY....SHOULD CHANGE TO JUST ZT
TCs_Mat_shared = rbind(G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12, G13, G14)
head(TCs_Mat_shared)
tail(TCs_Mat_shared) #Da was the last genotype to add so should be in the tail

#standardize the response scores such that the lowest value is greater than zero 
minVal = abs(min(TCs_Mat_shared)) + 1
Standard_TCs_Mat_shared = minVal + TCs_Mat_shared
head(Standard_TCs_Mat_shared)
head(TCs_Mat_shared)
range(Standard_TCs_Mat_shared)
range(TCs_Mat_shared)
min(Standard_TCs_Mat_shared)
min(TCs_Mat_shared)
which(Standard_TCs_Mat_shared<= 0)
which(TCs_Mat_shared<= 0)

which(is.na(TCs_Mat_shared)==TRUE)
tail(TCs_Mat_shared)

setwd(outPath)
save(TCs_Mat_shared, file = "AllButZeQu_TCs_Mat_Jan7.RData")
save(Standard_TCs_Mat_shared, file = "AllButZeQu_Standard_TCs_Mat_Jan7.RData")

#run response scores through network analyses
#BlockModsShared <- blockwiseModules(datExpr = t(TCs_Mat_shared), power = 18, networkType = "signed", corType="bicor", TOMType="signed", minModuleSize=100, mergeCutHeight=0.30, deepSplit=1, pamRespectsDendro = F, nThreads = 4, verbose=3)
#save(BlockModsShared, file = "AvSt_BlockModsShared.RData")


### OR JUST LOAD IN IF PREVIOUSLY RAN ### 
load(file.path(outPath, "AllButZeQu_Standard_BlockModsShared18_Dec29.RData"))

mergedColors = labels2colors(BlockModsShared$colors)
plotDendroAndColors(BlockModsShared$dendrograms[[1]], mergedColors[BlockModsShared$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                    main = "Module Relationships")

#identify how many modules were detected
dim(BlockModsShared[[3]])
MEs = (BlockModsShared[[3]])


#plot the module eigengenes to inspect patterns
#Add Zts to plot
MEs_Zts = rownames_to_column(MEs, var = "Zts")
colnames(MEs_Zts)
MEs_Zts$Zts

#need to add descriptor columns to plot with ggplot; will fix Geno later
MEs_Zts = MEs_Zts %>% separate(Zts, c("Geno", "Zts"))
#double check 
MEs_Zts$Zts

#lengthen
plottingMEs = pivot_longer(MEs_Zts, cols = starts_with("ME"), names_to = "Modules", values_to = "Expression")
str(plottingMEs)

plottingMEs$Zts = gsub("Zt*", "", plottingMEs$Zts)

plottingMEs = plottingMEs %>%
  mutate(Zts = as.numeric(Zts)) %>%
  arrange(plottingMEs)

#plot
plottingMEs %>%
  group_by(Modules) %>%
  ggplot(aes(x = as.numeric(Zts), 
             y = Expression)) +
  geom_line() +
  facet_wrap(~Modules) +
  xlab("Zt") +
  ylab("Eigengene....")


#pull gene per module information
geneByMod = as.data.frame(BlockModsShared$colors)
levels(as.factor(geneByMod$`BlockModsShared$colors`))
#change colors to a numeric ID of each module

### DOING THIS MANUALLY FOR NOW ###
geneByMod$`BlockModsShared$colors` = geneByMod$`BlockModsShared$colors` %>%
  str_replace_all(c("black" = "1", "blue" = "2", "brown" = "3", 
                    "green" = "4", "greenyelow" = "5", "grey" = "6",
                    "magenta" = "7", "pink" = "8", "purple" = "9", 
                    "red" = "10", "tan" = "11","turquoise" = "12", 
                    "yellow" = "13"))

unique(geneByMod$`BlockModsShared$colors`)

#change column name
colnames(geneByMod) = "Module"
#make numeric
geneByMod$Module = as.numeric(geneByMod$Module)
str(geneByMod)

#sort by module
geneByMod = geneByMod %>% arrange(geneByMod)
levels(as.factor(geneByMod$Module))
#check
which(geneByMod$Module == 1)
which(geneByMod$Module == 3)
which(geneByMod$Module == 7)

#make a gene id and genotype column
geneByMod = rownames_to_column(geneByMod, var = "Gene")
#remove the "gene:" prefix to oleracea gene names
tail(geneByMod)
geneByMod$Gene = gsub("gene:", "", geneByMod$Gene)
tail(geneByMod)

geneByMod = geneByMod %>% separate(Gene, c("Genotype", "Gene"))
#check
levels(as.factor(geneByMod$Genotype))
levels(as.factor(geneByMod$Module))
str(geneByMod)


#what proportion of each module belongs to each genotype? 
geneByMod = geneByMod %>% 
  group_by(Module) %>%
  mutate(GenesPerModule = length(Module)) %>%
  ungroup() %>%
  group_by(Genotype, Module) %>%
  mutate(GenesPerGeno = length(Module), 
         PropPerModule = GenesPerGeno/GenesPerModule) %>%
  ungroup() 

getwd()
write.csv(geneByMod, file = "AllButZeQu_Standardized_geneByMod.csv", row.names = FALSE)

PropPerMod = unique(geneByMod[, c(1, 3, 6)])
range(PropPerMod$PropPerModule)
write.csv(PropPerMod, file = "AllButZeQu_Standardized_ProportionPerMod.csv", row.names =  FALSE)


#how many of the sig. kME genes are now in each module? 
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Ab_Da/R_output/Ab", "Ab_SigkMEs.RData"))
Ab_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Ab_Da/R_output/Da", "Da_SigkMEs.RData"))
Al_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Al_Br_Ca/R_output/Br", "Br_SigkMEs.RData"))
Br_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Al_Br_Ca/R_output/Ca", "Ca_SigkMEs.RData"))
Ca_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Av_Se_St_Ze/R_output/Av", "Av_SigkMEs.RData"))
Av_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Av_Se_St_Ze/R_output/Se", "Se_SigkMEs.RData"))
Se_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Av_Se_St_Ze/R_output/St", "St_SigkMEs.RData"))
St_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Gr_Mu_Ne/R_output/Gr", "Gr_SigkMEs.RData"))
Gr_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Gr_Mu_Ne/R_output/Mu", "Mu_SigkMEs.RData"))
Mu_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Gr_Mu_Ne/R_output/Mu", "Mu_SigkMEs.RData"))
Mu_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Gr_Mu_Ne/R_output/Ne", "Ne_SigkMEs.RData"))
Ne_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Yu_Qu_DH20_DH12/R_output/Yu", "Yu_SigkMEs.RData"))
Yu_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Yu_Qu_DH20_DH12/R_output/Qu", "Qu_SigkMEs.RData"))
DH20_SigkMEs = names(SigkMEs)
load(file.path("~/Desktop/Git/RNAseq/Bnapus/Yu_Qu_DH20_DH12/R_output/DH12", "DH12_SigkMEs.RData"))
DH12_SigkMEs = names(SigkMEs)

rm(SigkMEs)

#bind together into a list
SigkMEs_List = list(Ab_SigkMEs, Al_SigkMEs, Av_SigkMEs, Br_SigkMEs, Ca_SigkMEs, DH12_SigkMEs, 
                    DH20_SigkMEs, Gr_SigkMEs, Mu_SigkMEs, Ne_SigkMEs, Se_SigkMEs, 
                    St_SigkMEs, Yu_SigkMEs, Da_SigkMEs)
names(SigkMEs_List) = c("Ab", "Al", "Av", "Br", "Ca", "DH12", "DH20", "Gr", "Mu", "Ne", "Se", "St", "Yu", "Da")


#loop through the list and pull out sig genes that match genes from the response network
  #first remove the "gene:" from the oleracea names
SigkMEs_List = lapply(SigkMEs_List, function(x) gsub("gene:", "", x))
SigkMEs_List[[1]][1892]

SigkME_responseGenes = lapply(SigkMEs_List, function(x) intersect(geneByMod$Gene, x))

#make sure the names match up between lists
names(SigkMEs_List)
names(SigkME_responseGenes)

#How many genes are not drought responsive
lengths(SigkMEs_List) - lengths(SigkME_responseGenes)

#what proportion are drought responsive
Prop_DrResp = lengths(SigkME_responseGenes) / lengths(SigkMEs_List)


#create and save some tables
SigkME_lengths = lengths(SigkMEs_List)
DrResp_lengths = lengths(SigkME_responseGenes)
shared_props = data.frame(SigkME_lengths, DrResp_lengths, Prop_DrResp)

#of the shared genes that are drought responsive, how many are rapa vs oleracea? 
Bra_DrResp = list()
Bo_DrResp = list()

for(i in 1:length(SigkME_responseGenes)) {
  Bra_DrResp[[i]] = SigkME_responseGenes[[i]][which(str_count(SigkME_responseGenes[[i]], "Bra") == 1)]
  Bo_DrResp[[i]] = SigkME_responseGenes[[i]][which(str_count(SigkME_responseGenes[[i]], "Bo") == 1)]
}

Bra_lengths = lengths(Bra_DrResp) / lengths(SigkME_responseGenes)
Bo_lengths = lengths(Bo_DrResp) / lengths(SigkME_responseGenes)

shared_props$BraProp = Bra_lengths
shared_props$BoProp = Bo_lengths

write.csv(shared_props, file = "shared_DrResp_props.csv")


#this will save all genes into a single list
lapply(SigkME_responseGenes, function(x) write.table(data.frame(x), 'SigkME_ResponseNetwork.csv', sep=',', row.names = FALSE))

#or can save each element into their own csv 
files <- c("Ab_SharedStandardGenes.csv", "Al_SharedStandardGenes.csv", "Av_SharedStandardGenes.csv", 
           "Br_SharedStandardGenes.csv", "Ca_SharedStandardGenes.csv", "DH12_SharedStandardGenes.csv",
           "DH20_SharedStandardGenes.csv", "Gr_SharedStandardGenes.csv", "Mu_SharedStandardGenes.csv", 
           "Gr_SharedStandardGenes.csv", "Ne_SharedStandardGenes.csv", "Se_SharedStandardGenes.csv", 
           "St_SharedStandardGenes.csv", "Yu_SharedStandardGenes.csv")
for(i in 1:length(SigkME_responseGenes)){
  write.csv(SigkME_responseGenes[[i]], file = files[i], row.names = FALSE)
}



#make another table with gene id, what geno that belongs to, which module it belongs to, and whether or not it is drought responsive 
test = data.frame(unlist(SigkMEs_List))
colnames(test) = "Gene"
test = rownames_to_column(test, var = "Genotype")
#get rid of the row number after each genotype ID
test$Genotype = gsub("[0-9]+$", "", test$Genotype)

levels(as.factor(test$Genotype))
str(test)
responseGeneNames = unlist(SigkME_responseGenes)
head(responseGeneNames)
length(responseGeneNames)
  ### DOUBLE CHECK THAT THIS NUMBER MAKES SENSE... ###

#add a binary column with a 1 if the gene is drought responsive or a 0 if not
test2$DroughtResponse <- ifelse(test2$Gene %in% responseGeneNames, 1 , 0)

  ### Obviously there are redundant genes in here...will have to pivot wider? Or just subset/grab unique? ###
length(which(test$DroughtResponse == 1))
length(which(test$DroughtResponse == 0))


#Finally add the module number for each gene
test3 = merge(test2, geneByMod, by = c("Gene", "Genotype"))


### taking a look at the ENSRNA genes..... ###
ens_genes = list()
for(i in 1:length(SigkMEs_List)) {
  ens_genes[[i]] = SigkMEs_List[[i]][str_detect(SigkMEs_List[[i]], "ENS")]
}

names(ens_genes) = c("Ab", "Da", "Al", "Av", "Br", "Ca", "DH12", "DH20", "Gr", "Mu", "Ne", "Se", "St", "Yu")
#save all genes 
files <- c("Ab_ENSRNAgenes.csv", "Da_ENSRNAgenes.csv",
           "Al_ENSRNAgenes.csv", "Av_ENSRNAgenes.csv", 
           "Br_ENSRNAgenes.csv", "Ca_ENSRNAgenes.csv", 
           "DH12_ENSRNAgenes.csv","DH20_ENSRNAgenes.csv", 
           "Gr_ENSRNAgenes.csv", "Mu_ENSRNAgenes.csv", 
           "Gr_ENSRNAgenes.csv", "Ne_ENSRNAgenes.csv", 
           "Se_ENSRNAgenes.csv", 
           "St_ENSRNAgenes.csv", "Yu_ENSRNAgenes.csv")
for(i in 1:length(ens_genes)){
  write.csv(ens_genes[[i]], file = files[i], row.names = FALSE)
}


### making sure there isn't any other weird annotation in here ###
weird_genes = list()
clean_SigkMEs_List = list()
for(i in 1: length(SigkMEs_List)){
  #pull out anything that doesn't start with Bra or Bo
  weird_genes[[i]] = SigkMEs_List[[i]][!str_detect(SigkMEs_List[[i]], "^Bra|^Bo*")]
  #remove any "weird genes"; probs all ENSRNA names
  clean_SigkMEs_List[[i]] = SigkMEs_List[[i]][which(SigkMEs_List[[i]] %in% weird_genes[[i]] == FALSE)]
}










