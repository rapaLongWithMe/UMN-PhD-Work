#Set Up 
require(tidyverse)

setwd("~/Desktop/Git/RNAseq/Brapa/For_KT/GO")


#GO_RData_fromStevan has Rfunctions, GO_TermDescriptions that was updated by Stevan, and the Nap_GO_BP_Annotations that Stevan also put together (contains both rapa and oleracea)
inPath = "./20221003"
#outPath = "./For_KT/20221003/"

#load some pre-defined Rfunctions
source(file = "Rfunctions.R")


##### SKIP TO LINE 41 #####
#load the drought responsive genes with module information
DR_Dat = read.csv(file.path(inPath, file = "ACC28_4hDEL_ZT0ctl_Bra.txt"), header = FALSE)
#make sure the "gene: *" prefix is not there
tail(DR_Dat)

#Load the expressed genes as the reference set
RefGenes = load("ACC28_TCs_List.RData")
RefGenes = rownames_to_column(data.frame(TCs_List$Cold_R1), var = "Gene")
RefGenes = RefGenes %>% separate(Gene, c("Geno", "Treat", "Gene")) %>% pull(Gene)
tail(RefGenes)


#Load the GO BP annotations 
EnrichGo_BP = read.csv("Nap_GO_BP_Annotations.txt",stringsAsFactors = F)

# This function is used to put the annotations in a useful format to calculate categorical enrichment
EnrichGo_BP = BuildEnrichMaps(inDF = EnrichGo_BP)

# Load the GO term descriptions
GODescriptions<-read.csv("GO_Term_Descriptions.txt",stringsAsFactors = F)
GODescriptions<-setNames(object = GODescriptions$Description, nm = GODescriptions$Goterm)


#run the enrichment
EnrichLst<-lapply(DR_Dat,function(x) CalcEnrich(MapBuild = EnrichGo_BP, testvec = x, codedesc = GODescriptions, refvec = RefGenes))

# Filter out all non-significant categories
EnrichLst_Sig.01<-lapply(EnrichLst,function(x) x@CatList[which(as.numeric(x@CatList[,"Adj_P-Value"]) <= 0.01),,drop=F])
EnrichLst_Sig.05<-lapply(EnrichLst,function(x) x@CatList[which(as.numeric(x@CatList[,"Adj_P-Value"]) <= 0.05),,drop=F])


#write the output to the output folder
write.csv(EnrichLst_Sig.01, file.path(outPath, file = "ACC28_EnrichLst_Sig.01.csv"))
write.csv(EnrichLst_Sig.05, file.path(outPath, file = "ACC28_EnrichLst_Sig.05.csv"))



