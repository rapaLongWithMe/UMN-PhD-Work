## SetUp
require(tidyverse)

##### CHANGE BASED ON GENOTYPE #####
setwd("~/Desktop/JGI_Cold/TestingPermPipeline/Pull_BestAligned_TPMs")

# Set the accession variable and the reference
Geno = "PC185"
Reference = "L58"
##### 

## Load and clean data 
tpms = read.delim(paste0("out/", Geno, "/", Reference, "_BestAlignment_tpm_counts", Geno, ".txt"))

# Fix the column names
colnms<-colnames(data.frame(tpms))
colnms

# C is freeze, W is warm, and TP is time point  
  # Use this one when running any of the 'ACC's through or the CC168 line
colnms<-gsub("CTP","Freeze_TP",colnms)
#colnms<-gsub("C","Freeze_",colnms)
colnms<-gsub("WTP","Control_TP",colnms)

# Split into pieces, each piece will become a list used to index
rawnames = strsplit(x = colnms, split = "_")
# Making a new vector with names ordered as Geno_Replicate#_Treatment_Timepoint
short_names = sapply(rawnames,function(x) paste(x[c(1,4,2,3)],collapse = "_"))

#set cleaned column names to previous df
colnames(tpms)<-short_names
colnames(tpms)

# Remove time point 7 and 8 since those are actually recovery
clean = tpms %>%
  select(!contains(c("TP8", "TP7")))
colnames(clean)

# Make GeneName the rowname
rownames(clean)
  # Make sure there aren't any genes from contamination
clean = clean[which(str_detect((clean$Geneid_NA_NA_NA), "Br") == TRUE), ]
clean = column_to_rownames(clean, var = "Geneid_NA_NA_NA")

# Write out
norm.counts = clean 
colnames(norm.counts) = gsub("TP", "", colnames(norm.counts))
rownames(norm.counts)
colnames(norm.counts)
write.csv(norm.counts, file = paste0("ReadyForPermPipe/", Geno, "_norm.counts.csv"))



### SPLIT THE dipalm OBJECT AFTER PERMUTING ### 
load("/Users/ricon001/Desktop/Git/BrapaPangenome_ColdAcclimation/DiPALM_PermResults/ZCT/dipalm.RData")
setwd("~/Desktop/Git/BrapaPangenome_ColdAcclimation/DiPALM_PermResults/ZCT/")
write.csv(sig.genes, file = "ZCT_SigDiPALMgenes.csv", row.names = FALSE)






