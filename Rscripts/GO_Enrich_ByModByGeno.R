#Set Up 
require(tidyverse)

setwd("~/Desktop/Git/RNAseq/Bnapus/SharedResponse")


#GO_RData_fromStevan has Rfunctions, GO_TermDescriptions that was updated by Stevan, and the Nap_GO_BP_Annotations that Stevan also put together (contains both rapa and oleracea)
inPath = "~/Desktop/Git/RNAseq/Bnapus/SharedResponse/GO_Enrichment"
outPath = "~/Desktop/Git/RNAseq/Bnapus/SharedResponse/GO_Enrichment/20220307"

#load some pre-defined Rfunctions
source(file.path(inPath, file = "Rfunctions.R"))


##### SKIP TO LINE 41 #####
#load the drought responsive genes with module information
DR_Dat = read.csv("DRgenes_ByGeno_WithResponses_AndMods.csv")
  #remove the first column of erroneous rownames; will have to fix that
DR_Dat = DR_Dat[, -1]

#split the ID column to pull out the Genotype of interest
DR_Dat = DR_Dat %>% separate(IDs, c("Geno", "Gene"))

#make sure the "gene: *" prefix is not there
tail(DR_Dat)
#double check a few other things
levels(as.factor(DR_Dat$Geno))
  #should NOT have a grey module in here...
levels(as.factor(DR_Dat$Module))

#remove the response scores
DR_Dat = DR_Dat[, c(1:2, 7)]

#write out the DR_Dat into the GO_Enrichment folder so don't have to repeat all that again
getwd()
write.csv(DR_Dat, file.path(outPath, file = "DR_Dat_GO_Enrich_March7.csv"), row.names = FALSE)
######################################
##### AFTER FIRST RUN START HERE #####
######################################

DR_Dat = read.csv(file.path(outPath, file = "DR_Dat_GO_Enrich_March7.csv"))


##### CHANGE BASED ON GENOTYPE #####
DR_byGeno = DR_Dat %>% filter(Geno == "Yu")
  ## make sure the number of obs matches the MetaData sheet! 

#and break this into DR genes by module
ModList = split(DR_byGeno, DR_byGeno$Module)
  #check the number of genes by module in the PropTable csv (under MS_Figures)


#Now that they're split up, just pull out the gene names
DR_List = lapply(ModList, function(x) x$Gene)

#Load the GO BP annotations 
EnrichGo_BP = read.csv(file.path(inPath, "Nap_GO_BP_Annotations.txt"),stringsAsFactors = F)

# This function is used to put the annotations in a useful format to calculate categorical enrichment
EnrichGo_BP = BuildEnrichMaps(inDF = EnrichGo_BP)

# Load the GO term descriptions
GODescriptions<-read.csv(file.path(inPath, "GO_Term_Descriptions.txt"),stringsAsFactors = F)
GODescriptions<-setNames(object = GODescriptions$Description, nm = GODescriptions$Goterm)


#pull all expressed genes as a reference set to test against
DR_GeneNames <- DR_Dat$Gene
  #again make sure that "gene: " prefix is gone otherwise won't work properly 
tail(DR_GeneNames)


#run the enrichment
EnrichLst<-lapply(DR_List,function(x) CalcEnrich(MapBuild = EnrichGo_BP, testvec = x, codedesc = GODescriptions, refvec = DR_GeneNames))

# Filter out all non-significant categories
EnrichLst_Sig.01<-lapply(EnrichLst,function(x) x@CatList[which(as.numeric(x@CatList[,"Adj_P-Value"]) <= 0.01),,drop=F])


#write the output to the output folder
  ##### CHANGE BASED ON GENOTYPE #####
dir.create(path = file.path(outPath,"Yu_DR_GeneEnrichment_March7"))
sapply(names(EnrichLst_Sig.01), function(x) 
  write.csv(EnrichLst_Sig.01[[x]],file=file.path(outPath,"Yu_DR_GeneEnrichment_March7",paste(x,".csv",sep=""))))
#####  #####  #####  #####  #####  #####


#read the file back in but just keep the first five rows. this will be used to compile a master sheet of all genos by mod and the top five GO terms
 
##### CHANGE BASED ON GENOTYPE #####

setwd("~/Desktop/Git/RNAseq/Bnapus/SharedResponse/GO_Enrichment/20220307/Yu_DR_GeneEnrichment_March7")
##### ##### ##### ##### ##### ##### #####


enrich_files <- list.files(".", full.names = T, recursive = T)
enrich_files
  #grab just the top five (in this case first five) GO terms
enrich_list <- lapply(enrich_files, function(x) read.csv(paste0(x), nrow = 5))
  names(enrich_list) = names(DR_List)
  
#Add a module id column to each element of the list
enrich_master = list()
for(i in 1:length(enrich_list)) {
enrich_master[[i]] = enrich_list[[i]] %>%
  mutate(Module = names(enrich_list[i]))
}  
names(enrich_master) = names(enrich_list)

#bind all the elements together
enrich_df <- do.call("rbind", enrich_master)
glimpse(enrich_df)
  #rename X to be GO Term
colnames(enrich_df)[1] = "GO_Term"
colnames(enrich_df)
rownames(enrich_df)
  #probs want to not write out rownames since they'll end up duplicated once you bind all the genos together

##### CHANGE BASED ON GENOTYPE #####
#add a genotype id column
enrich_df$Geno = "Yu"

#save this with a genotype prefix so can bring all in later and bind into one master sheet
setwd(outPath)
write.csv(enrich_df, file = "Yu_enrich_df_March7.csv", row.names = FALSE)
##### ##### ##### ##### ##### #####



##### NOW WANT TO COMBINE ALL THE INDIVIDUAL CSV's INTO A SINGLE MASTER SHEET 

#compile 
  ##### MAKE SURE the files you want to combine are the only csv files in this wd!! #####
data_all <- list.files(path = "~/Desktop/Git/RNAseq/Bnapus/EnrichedGenes_GO_Terms/20220307", pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                                            
  bind_rows   

levels(as.factor(data_all$Geno))
levels(as.factor(data_all$Module))
levels(as.factor(data_all$GO_Term))


#write it back out as a final master sheet
getwd()
write.csv(data_all, file = "20220307_GO_MasterSheet_TopFive.csv", row.names = FALSE)

#load back in if starting here
topFive = read.csv(file = "~/Desktop/Git/RNAseq/Bnapus/EnrichedGenes_GO_Terms/AnnaC_20220323_GO_MasterSheet_TopFive.csv")
#make sure these are all over enriched
levels(as.factor(topFive$Enrichment.Direction))

#pull out just the columns you're interested in 
topFive = topFive %>% select(GO_Term, Module, Geno) 
#split off the GO term and the description into different columns to make it easier to sort
topFive = topFive %>% separate(GO_Term, into = c("GO_Term", "Description"), sep = " - ")


#reorganize
topFive = topFive %>% pivot_wider(names_from = Geno, values_from = Description)
#split such that each module is its own sheet
topFive_mods = split(topFive, topFive$Module)

#write out so each mod is it's own sheet 
names(topFive_mods)
getwd()

for(i in names(topFive_mods)){
  write.csv(topFive_mods[[i]], paste0(i,".csv"), row.names = FALSE)
}
