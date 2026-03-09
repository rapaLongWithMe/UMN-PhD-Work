## SetUp
require(tidyverse)

##### CHANGE BASED ON GENOTYPE #####
setwd("~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses/PreOct2022/Yu_Qu_DH20_DH12")

#Define some path variables
inPath = "./R_input"
outPath = "./R_output"

## Load and clean data 
### CHANGE THE FILENAME TO MATCH GENOTYPE ### 
load(file.path(inPath, file = "Yu_NormCounts.RData"))
##### ##### ##### ##### ##### ##### ##### #####

which(is.na(NormCounts))

#fix the column names
DF = data.frame(NormCounts)
colnms<-colnames(data.frame(NormCounts))
colnms

#remove the bam. and .bam 
colnms<-gsub("bam.","",colnms)
colnms<-gsub(".bam","",colnms)
#C is cold, W is warm, and TP is time point
colnms<-gsub("100","Watered_",colnms)
colnms<-gsub("20","Drought_",colnms)

#split into pieces, each piece will become a list used to index
rawnames = strsplit(x = colnms, split = "_")
#making a new vector with names ordered as Geno_Replicate#_Treatment_Timepoint
short_names = sapply(rawnames,function(x) paste(x[c(1,4,2,3)],collapse = "_"))

#set cleaned column names to previous df
colnames(DF)<-short_names
colnames(DF)

#Remove time point 7 and 8 since those are actually recovery
clean = DF %>%
  select(!contains(c("_8", "_7")))
#make GeneName the rowname
rownames(clean)

clean = clean[which(str_detect(rownames(clean), "Bra") == TRUE), ]

#make a matrix; the number of elements should match the number of elements present in NormCounts
CleanGenes = as.matrix(clean)


##Remove non-expressed genes
  #plot to determine a cutoff
geneMeans<-apply(CleanGenes,1,function(x) mean(x,na.rm = T))
hist(geneMeans,col="skyblue", breaks = seq(-12,16,0.25), main = "GeneMeans")
abline(v=0,col="red",lwd=3,lty=2)

  #remove genes that do not have at least one sample with mean log2 FPKM >0
FilteredGenes<-CleanGenes[which(geneMeans>0),]
  #plot again to make sure retaining genes with mean >0
geneMeans<-apply(FilteredGenes,1,function(x) mean(x,))
hist(geneMeans,col="skyblue", breaks = seq(-12,16,0.25), main = "FPKM >0")
abline(v=0,col="red",lwd=3,lty=2)
  #check that the number/proportion of genes removed makes sense
(nrow(CleanGenes) - nrow(FilteredGenes)) / nrow(CleanGenes)

  #replace NaNs resulting from taking the log of 0 (when normalizing counts)
ProcessedGenes = FilteredGenes
minVal<-min(FilteredGenes[!is.na(FilteredGenes)])-1
ProcessedGenes[is.na(ProcessedGenes)]<-minVal
  #double check
which(is.na(FilteredGenes))
which(is.na(ProcessedGenes))
  #Final plot to double check; depending on the number of NaNs you start with you may see a small tail below zero
geneMeans<-apply(ProcessedGenes,1,function(x) mean(x,na.rm = T))
hist(geneMeans,col="skyblue", breaks = seq(-12,16,0.25), main = "Processed Genes")
abline(v=0,col="red",lwd=3,lty=2)

#Save the processed data
  #want genes as rows and samples as columns
rownames(ProcessedGenes)
colnames(ProcessedGenes)

##### CHANGE BASED ON GENOTYPE #####
write.csv(ProcessedGenes, file = file.path(outPath, "ACC28_ProcessedGenes.csv"))
##### ##### ##### ##### ##### ##### #####


#can just load the pre-processed genes if starting here
#load(file.path(outPath, "Av_ProcessedGenes.RData"))

## Impute missing samples if necessary
#### CHANGE THIS BASED ON GENOTYPE!! ####
ImputedGenes = data.frame(ProcessedGenes)

colnames(ImputedGenes)
ImputedGenes$Ne_7_Drought_R1 = ImputedGenes %>%
  select(Ne_7_Drought_R2, Ne_7_Drought_R3) %>%
  apply(1, median)
colnames(ImputedGenes)

ImputedGenes = as.matrix(ImputedGenes)

#remove genes that do not have at least 3 samples with log2(fpkm)>1
  #make a logical table (T/F) to ID when FPKM is not >1
#ImputedGenes = as.matrix(ImputedGenes)
  ### CHANGE TO MATCH YOUR REQUIRIMENT ###
#fpkm1 = ImputedGenes > 1
  #for each gene, sum the number of TRUE occurances
#fpkm_sum = apply(fpkm1, 1, function(x) sum(x, na.rm = TRUE))
  #now subset to only include 3+ TRUE occurances
#fpkm2 = ImputedGenes[fpkm_sum >= 3, ]
#colnames(fpkm2)


#remove genes with no variation in expression
  ##first separate data into TimeCourses by treatment DroughtR1, DroughtR2, DroughtR3; same for watered
tmp<-ProcessedGenes

#split the names by a separator
spNms<-strsplit(x = colnames(tmp), split = "_")
spNms
#take the full names (Geno_Rep_Treat_Zt) and collapse them to Treatment.R1
tnms<-sapply(spNms,function(x) paste(x[c(3,2)],collapse = "_"))
#separate by sample type
  ### MAKE SURE THIS TABLE HAS EXPRESSIOIN DATA AND NOT TRUE/FALSE ###
FilteredTCs<-tapply(1:length(tnms),INDEX = tnms, function(x) tmp[,x])

#find what genes = TRUE for var > 0
varFiltered<-lapply(FilteredTCs,function(x) apply(x,1,function(y) var(y)>0))
#bind into a matrix of TRUE instances (for above) by sample type
varFiltered<-do.call(cbind,varFiltered)
#takes logical data and finds which vectors have all TRUEs
varFiltered<-apply(varFiltered,1,function(x) all(x))
#adding gene names back
varFiltered<-names(varFiltered)[which(varFiltered)]
#applying the above list to the previous data (seperated by TC; all tps for each replicate by treatment)
TCs<-lapply(FilteredTCs,function(x) x[varFiltered,])

## Reorder by Zt (here y = 4; or the piece of the colname that indicates Zt)
  #downstream will need both a list with elements as TCs (colnames as timepoint)
TCs_List<-lapply(TCs,function(x) x[,order(as.numeric(sapply(strsplit(colnames(x),split = "_"),function(y) y[4])))])
#want one TC as an element and all timepoints as columns for each element
sapply(TCs_List,colnames)

#compile into a large matrix
TCs_Mat<-do.call(rbind,TCs_List)

#check the gene names
head(rownames(TCs_Mat))
tail(rownames(TCs_Mat))
#remove any annotation (ie ENSRNA genes) that are not Bra 
which(!str_detect(rownames(TCs_Mat), "Bra")== TRUE)

#Add a treatment prefix to every gene
Cold = TCs_List[1:4]
names(Cold)
Control = TCs_List[5:8]
names(Control)

#Create rownames for the two lists
Cold_names = paste("Cold", rownames(Cold[[1]]), sep = "_")
Control_names = paste("Control", rownames(Control[[1]]), sep = "_")

#update the rownames
Cold = lapply(Cold, function(x) {rownames(x) = Cold_names; x})
Control = lapply(Control, function(x) {rownames(x) = Control_names; x})


#bind them back together
test = list(Cold, Control)
TCs_List = unlist(test, recursive = FALSE)

#want one TC as an element and all timepoints as columns for each element
sapply(TCs_List,colnames)

#compile into a large matrix
TCs_Mat<-do.call(rbind,TCs_List)

#check the gene names
head(rownames(TCs_Mat))
tail(rownames(TCs_Mat))



##### CHANGE TO MATCH GENOTYPE ##### 
filepath = file.path("~/Desktop/Git/RNAseq/Brapa/ExpressedGenes/TCs_Lists")
save(TCs_List, file = file.path(filepath, "ACC28_TCs_List.RData"))

filepath = file.path("~/Desktop/Git/RNAseq/Brapa/ExpressedGenes/TCs_Mat")
save(TCs_Mat, file = "ACC28_TCs_Mat.RData")
##### ##### ##### ##### ##### #####



