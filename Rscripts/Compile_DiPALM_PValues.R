# Compiling all of the individual kME and kMed lists (pvalues for every expressed gene from DiPALM; per genotype)

# Set up
setwd("~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes")

kMEPath = "~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes/kME_genes/AdjkMEs"
outPath = "~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes/kME_genes"

#bring in all the files and add the file name to separate
all.files = list.files(file.path(kMEPath))
all.files

kMEs <- lapply(all.files, function(x) get(load(file.path(kMEPath, file = x))))
#really just interested in the Genotype ID
all.files = gsub("_AdjkMEs.RData", "", all.files)
all.files
names(kMEs) <- all.files 

#for now, not interested in Qu or Ze as they had odd permutation plots
kMEs = kMEs[-c(12, 16)]

#pull out gene name and geno information into dataframes
pValues_kME = list()
for(i in 1:length(kMEs)) {
  #make each element a dataframe
  pValues_kME[[i]] = data.frame(AdjkME_p = kMEs[[i]])
  #pull the rownames in to get a column of gene names
  pValues_kME[[i]] = rownames_to_column(pValues_kME[[i]], var = "Gene")
  #paste the genotype ID as a prefix to each gene
  pValues_kME[[i]]$Geno = paste0(names(kMEs)[i])
}
lapply(pValues_kME, function(x) levels(as.factor(x$Geno)))
names(pValues_kME) = names(kMEs)

#unite the Geno and Gene columns into a unique identifier
pValues_kME = lapply(pValues_kME, function(x) x %>% unite(IDs, c("Geno", "Gene"), sep = "_"))

pValues_kMETable = do.call(rbind, pValues_kME)

#save the table and the RObject
setwd(outPath)
save(pValues_kMETable, file = "All_AdjkME_pValues_April25.RData")
write.csv(pValues_kMETable, file = "All_AdjkME_pValues_Apri25.csv", row.names = FALSE)


## Now bring in the kMed data and repeat
MedPath = "~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes/Med_genes/AdjMed"
outPath = "~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes/Med_genes"

#bring in all the files and add the file name to separate
all.files = list.files(file.path(MedPath))
all.files

Meds <- lapply(all.files, function(x) get(load(file.path(MedPath, file = x))))
#really just interested in the Genotype ID
all.files = gsub("_AdjMed.RData", "", all.files)
all.files
names(Meds) <- all.files 

#remove Qu or Ze 
Meds = Meds[-c(12, 16)]

pValues_Med = list()
for(i in 1:length(Meds)) {
  #make each element a dataframe
  pValues_Med[[i]] = data.frame(AdjMed_p = Meds[[i]])
  #pull the rownames in to get a column of gene names
  pValues_Med[[i]] = rownames_to_column(pValues_Med[[i]], var = "Gene")
  #paste the genotype ID as a prefix to each gene
  pValues_Med[[i]]$Geno = paste0(names(Meds)[i])
}
lapply(pValues_Med, function(x) levels(as.factor(x$Geno)))
names(pValues_Med) = names(Meds)
#unite the Geno and Gene columns into a unique identifier
pValues_Med = lapply(pValues_Med, function(x) x %>% unite(IDs, c("Geno", "Gene"), sep = "_"))

#make a table
pValues_MedTable = do.call(rbind, pValues_Med)


#save the table and the RObject
setwd(outPath)
save(pValues_MedTable, file = "All_AdjMed_pValues_April25.RData")
write.csv(pValues_MedTable, file = "All_AdjMed_pValues_April25.csv", row.names = FALSE)

#combine the two tables into a master dataframe
  ##### MAKE SURE THIS MATCHES THE METADATA; SHOULD BE 592,429 SINCE ADJKME AND ADJMED HAVE ALL GENES
pValues = full_join(pValues_kMETable, pValues_MedTable)

  #MetaData has 111,704 total
Sig_pValues.01 = pValues %>% filter(AdjMed_p <= 0.01 | AdjkME_p <= 0.01)
tail(Sig_pValues.01)
  #remove the "gene:" prefix
Sig_pValues.01$IDs = gsub("gene:", "", Sig_pValues.01$IDs)
tail(Sig_pValues.01)
  
  #save
setwd("~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes")
save(pValues, file = "All_kME_kMed_pValues_April4.RData")
write.csv(pValues, file = "All_kME_kMed_pValues_April4.csv")
save(Sig_pValues.01, file = "Sig_pValues.01_April4.RData")
write.csv(Sig_pValues.01, file = "Sig_pValues.01_April4.csv")


