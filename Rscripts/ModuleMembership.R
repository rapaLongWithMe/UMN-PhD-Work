### MODULE MEMBERSHIP

#Set Up
require(WGCNA)
require(DiPALM)
require(tidyverse)

##### CHANGE BASED ON GENOTYPE #####
setwd("~/Desktop/Git/RNAseq/Brapa")
outPath = "~/Desktop/Git/RNAseq/Brapa/ACC28/R_output"
FigPath = "~/Desktop/Git/RNAseq/Brapa/ACC28/Figures"

#load the pre-processed data
load(file.path(outPath, "ACC28_BlockModsAll.RData"))
load(file.path(outPath, "ACC28_TCs_Mat.RData"))
load(file.path(outPath, "ACC28_TCs_List.RData"))
##### ##### ##### ##### ##### #####

#extract eigengenes
MEs<-BlockModsAll[[3]]
#check that the number of expressed matches the MetaData table
  #divide by the #replicates*#treatments
length(names(BlockModsAll$colors))/8
#check that the number of co-expression modules matches the MetaData table
ncol(MEs)

#calculate membership
## where kME is the pattern change difference, and Med is the difference in median expression
kMEsList<-BuildModMembership(MeMat = MEs, TCsLst = TCs_List)
Med<-sapply(TCs_List,function(x) apply(x,1,function(y) median(y,na.rm = T)))
#permute non-merged dataset to estimate the null distribution
Perm<-lapply(TCs_List,function(x) x[sample(1:nrow(x),nrow(x),replace = T),])
#calculate kME and kMed scores from permuted data
kMEsPerm<-BuildModMembership(MeMat = MEs, TCsLst = Perm)
MedPerm<-sapply(Perm,function(x) apply(x,1,function(y) median(y,na.rm = T)))

#construct linear contrasts
#make sure the Treat vector matches the samples (replicates) 
Treat<-as.factor(c("Cold","Cold","Cold", "Cold",
                   "Control","Control","Control", "Control"))
design<-model.matrix(~0+Treat)
colnames(design)<-levels(Treat)
contr<-"Cold-Control"

LimmaModskMEs<-lapply(kMEsList, function(x) BuildLimmaLM(dataMat = x, designMat = design, contrastStr = contr))
LimmaModsMed<-BuildLimmaLM(dataMat = Med, designMat = design, contrastStr = contr)

#pull out t-scores; clear memory
LimmaModskMEs<-do.call(cbind,lapply(LimmaModskMEs,function(x) x$t))
LimmaModsMed<-LimmaModsMed$t
gc()

#repeat limma tests on permuted data 
LimmaModskMEsPerm<-lapply(kMEsPerm, function(x) BuildLimmaLM(dataMat = x, designMat = design, contrastStr = contr))
LimmaModsMedPerm<-BuildLimmaLM(dataMat = MedPerm, designMat = design, contrastStr = contr)
LimmaModskMEsPerm<-do.call(cbind,lapply(LimmaModskMEsPerm,function(x) x$t))
LimmaModsMedPerm<-LimmaModsMedPerm$t
gc()

#get the absolute value of the test stats (used to plot distributions, below)
TestSumskMEs<-apply(LimmaModskMEs,1, function(x) sum(abs(x),na.rm = T))
TestSumsMed<-abs(LimmaModsMed[,1])
#repeat on permuted data
PermSumskMEs<-apply(LimmaModskMEsPerm,1, function(x) sum(abs(x),na.rm = T))
PermSumsMed<-abs(LimmaModsMedPerm[,1])


##### CHANGED BASED ON GENOTYPE #####
setwd(outPath)

save(TestSumskMEs, file = "ACC28_TestSumskMEs.RData")
save(TestSumsMed, file = "ACC28_TestSumsMed.RData")
save(LimmaModskMEs, file = "ACC28_LimmaModskMEs.RData")
save(LimmaModsMed, file = "ACC28_LimmaModsMed.RData")
##### ##### ##### ##### ##### #####

#plot both patterns datasets to determine a cut off
ggPlotMultiDensities(denslist = list(Test=TestSumskMEs,Permuted=PermSumskMEs), main = "Pattern Change Scores", xlab = "Differential Pattern Score",lwidth = 1, cols = c("chartreuse4", "darkgoldenrod1"))
ggPlotMultiDensities(denslist = list(Test=TestSumsMed,Permuted=PermSumsMed), main = "Expression Change Scores", xlab = "Differential Expression Score",lwidth = 1, cols = c("chartreuse4", "darkgoldenrod1"))

#determine a significance cutoff: calculate a pValue using the individual values from each test sum and permuted test sum with an FDR correction.
  #this will be saved later and used to pull out dr genes from both kME and kMed lists
AdjkMEs<-sapply(TestSumskMEs,function(x) AdjustPvalue(tVal = x, tVec = TestSumskMEs, pVec = PermSumskMEs))
AdjMed<-sapply(TestSumsMed,function(x) AdjustPvalue(tVal = x, tVec = TestSumsMed, pVec = PermSumsMed))

  #pull out genes with p > 0.01
SigkMEs<-AdjkMEs[which(AdjkMEs<0.01)]
SigMed<-AdjMed[which(AdjMed<0.01)]
  
  #number of differently patterned genes
length(SigkMEs)
length(SigMed)

  
#pull out total drought responsive genes
  #if a gene is found in both SigkMEs and SigMed, only count once (unique genes)
dr_names = c(names(SigkMEs), names(SigMed))
  #pull out rows with duplicated genes
dr_dup_genes = dr_names[which(duplicated(dr_names) == TRUE)]
  #now retain unique gene names that are drought responsive in some way
length(dr_names) - length(dr_dup_genes)
length(unique(dr_names))
  #the length of dr_genes should match both of these
dr_genes = union(dr_names, dr_dup_genes)

#check again
head(dr_names)
head(dr_dup_genes)
  #this should look the same as head(dr_names) but NOT head(dr_dup_genes)
head(dr_genes)


#Plot a single gene to double check 
  #rapa 
r_cca1.1 = "BraA05g01930R"
r_toc1.1 = "BraA03g42600R"
  #oleracea
b_cca1.1 = "gene:Bo4g006930"
b_toc1.1 = "gene:Bo7g098570"

##### CHANGE THE FILE NAME BASED ON GENOTYPE #####
pdf(file.path(FigPath,"DH12_CCA1_rapa.pdf"),width = 4,height = 3)
PlotTCs(TClst = TCs_List,tgene = r_cca1.1, 
        scale = T, xlab="ZT Time (hours)", xAxsLabs = c(seq(1,23,4),seq(1,23,4)), 
        ledgeX="topright",
        tcols = c("darkgoldenrod1","darkgoldenrod1", "darkgoldenrod1",
                        "chartreuse4", "chartreuse4", "chartreuse4"), 
        tltys = c(1,2, 3,1,2, 3), 
        main = "B. rapa: cca1")
dev.off()

pdf(file.path(FigPath,"DH12_CCA1_oleracea.pdf"),width = 4,height = 3)
PlotTCs(TClst = TCs_List,tgene = b_cca1.1, 
        scale = T, xlab="ZT Time (hours)", xAxsLabs = c(seq(1,23,4),seq(1,23,4)), 
        ledgeX="topright",
        tcols = c("darkgoldenrod1","darkgoldenrod1", "darkgoldenrod1",
                  "chartreuse4", "chartreuse4", "chartreuse4"), 
        tltys = c(1,2, 3,1,2, 3), 
        main = "B. oleracea: cca1")
dev.off()

pdf(file.path(FigPath,"DH12_TOC1_rapa.pdf"),width = 4,height = 3)
PlotTCs(TClst = TCs_List,tgene = r_toc1.1, 
        scale = T, xlab="ZT Time (hours)", xAxsLabs = c(seq(1,23,4),seq(1,23,4)), 
        ledgeX="topright",
        tcols = c("darkgoldenrod1","darkgoldenrod1", "darkgoldenrod1",
                  "chartreuse4", "chartreuse4", "chartreuse4"), 
        tltys = c(1,2, 3,1,2, 3), 
        main = "B. rapa: toc1")
dev.off()

pdf(file.path(FigPath,"DH12_TOC1_oleracea.pdf"),width = 4,height = 3)
PlotTCs(TClst = TCs_List,tgene = b_toc1.1, 
                    scale = T, xlab="ZT Time (hours)", xAxsLabs = c(seq(1,23,4),seq(1,23,4)), 
                    ledgeX="topright",
                    tcols = c("darkgoldenrod1","darkgoldenrod1", "darkgoldenrod1",
                              "chartreuse4", "chartreuse4", "chartreuse4"), 
                    tltys = c(1,2, 3,1,2, 3), 
                    main = "B. oleracea: toc1")
dev.off()
##### ##### ##### ##### ##### ##### ##### ##### #####


##### CHANGE BASED ON GENOTYPE ##### ##### ##### #####

#save the full lists in case you want to adjust the cutoff later and the significant ones
save(AdjkMEs, file = file.path(outPath, "ACC28_AdjkMEs.RData"))
save(SigkMEs, file = file.path(outPath, "ACC28_SigkMEs.RData"))

save(AdjMed, file = file.path(outPath, "ACC28_AdjMed.RData"))
save(SigMed, file = file.path(outPath, "ACC28_SigMed.RData"))

#save which genes are both sig. patterned and sig. different in median expression
setwd("~/Desktop/Git/RNAseq/Bnapus/DroughtResponse_Genes")
save(dr_dup_genes, file = "DH12_SigInBoth_kMEandMedGenes.RData")
#and just the unique gene names that are either, or, or both
save(dr_genes, file = "DH12_DR_Genes.RData")

##### ##### ##### ##### ##### #####


