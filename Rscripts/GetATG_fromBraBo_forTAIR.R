## Take a list of Bra/Bo genes and find the At ortholog to copy into TAIR

## Set Up 
require(tidyverse)
require(readxl)
require(openxlsx)
require(readr)
BiocManager::install("xlsx", force = TRUE)

setwd("~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses/GO")

## Bring in the list of Br genes
  # Here I am using napus genes enriched in prolonged drought, by broad descript/Module/Geno
Enriched = read_xlsx("20240415_GSLgenes_ToPlot.xlsx",
                     #, col_names = FALSE)
                     ##### CHANGE TO THE SHEET OF INTEREST #####
                     sheet = "GSL_Anns")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### This chunk is useful if the Bra's are already split up and not in the original format (one string with mult. ann's separated by " | " )
  # Here I am pulling in the entire 50k GeneAnns_ByGenoMod csv that has all enriched annotations (one per row) for each broad descriptor, term, geno and module
Enriched = read.csv("GeneAnns_ByGenoMod.csv")
  # Interested in sulfur so start there
Enriched = Enriched %>% filter(BroadDescriptor == "sulfur/glucosinolate metabolism")

# Start by just pulling the unique Br names and running through
Enriched = Enriched[,3]

# Remove duplicates and sort; skip if you want to paste the ATG back to the previous DF at the end (ie keep the same rownumbers)
Enriched = unique(Enriched)
# Make a df for now since it's easier to read and check the At orthos
Enriched = data.frame(Br_Name = sort(Enriched))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



### Use this if you have the original output that has the annoying " | " between annotations (which are compressed into one string per row)
Enriched = Enriched %>% separate_rows(Gene, sep = " \\| ")

# Remove duplicates and sort; skip if you want to paste the ATG back to the previous DF at the end (ie keep the same rownumbers)
EnrichedGenes = unique(Enriched$Gene)
# Make a df for now since it's easier to read and check the At orthos
Enriched = data.frame(Br_Name = sort(Enriched))


# Bring in the At ortho's for Bra and match those in
At_orthos = read.csv("~/Desktop/RNAseq/KT_R500_Bro_Ath_Orthologues_NonsyntRemoved.csv")
At_orthos = At_orthos[, -c(1:4,7:8,10)]

# Merge to have a link of which paralog goes with which At ortholo

# Also useful for looking at paralogs from the other subgenome that may not be enriched
bra_key = At_orthos[which(At_orthos$BRA %in% Enriched$Br_Name),]
bo_key = At_orthos[which(At_orthos$B.oleracea_ortho %in% Enriched$Br_Name),]
  # keep the duplicates since the paralogs might differ
gsl_key = rbind(bra_key, bo_key)


# Pull the At ortho and write out
ATG = list()
for (i in 1:length(At_orthos)) {
  # subset the dataframe by matches for each column
  ATG[[i]] = At_orthos[which(At_orthos[[i]] %in% Enriched$Br_Name),] %>% pull(Ath_ortho)
  
}



# To use the write.xlsx function to write mult. sheets to the same file needs to be a df
# can't set colnames to NULL cause that creates a chr vector, so will have to manually remove after writing
ATG = data.frame(ATG = unique(unlist(ATG)))



##### For now save as a separate xlsx
write.xlsx(ATG, "AlL_enrichedGSLs_ATGs.xlsx",
           sheetName = "All_Smetab", rownames = FALSE, append = TRUE)




##### CLOSE TO WRITING EACH DF AS AN INDIVIUAL XLSX TAB #####

# Bring in prev file
Master_ATG = read.xlsx("Enriched_ATG_GeneDescriptions/AtOrthos_AnnByModGeno.xlsx")

## For each new tab add a new worksheet/name
wb <- createWorkbook()
addWorksheet(wb, "Smetab_blue_Ca")
addWorksheet(wb, "Smetab_blue_Mu")
addWorksheet(wb, "Smetab_blue_St")

## And write the data in appropriately
writeData(wb, "Smetab_blue_Ca", Master_ATG, startRow = 1, startCol = 1)
writeData(wb, "Smetab_blue_Mu", ATG, startRow = 1, startCol = 1)
writeData(wb, "Smetab_blue_St", ATG, startRow = 1, startCol = 1)


## Save
saveWorkbook(wb, file = "Enriched_ATG_GeneDescriptions/Master_ATG.xlsx", overwrite = TRUE)
##### ##### ##### ##### ##### ##### ##### 







