### Compiling a table of syntenic TFs for six B.rapa morphotypes

## Set up the environment
require(tidyverse)
require(readr)
require(readxl)

setwd("~/Desktop/JGI_Cold")


## Bring in the data
  # Bring in the list of clock genes to add as potential regulators
clock = read_xlsx("KT_Brapa_ClockGene_List.xlsx")
  # Remove the excess column with notes
clock = clock[, -8]
  # Pull out just the columns needed for this analyses
clock = clock[, c(4, 6)]
colnames(clock) = c("Gene", "CommonName")
  # Pull unique since there are several duplicates resulting from multiple Bra paralogs
clock = unique(clock)
#write.csv(clock, file = "SyntenicGenes/UniqueClockATGs.csv", row.names = FALSE)

## Bring in Chiemeka's 'GRN_ready' tables
A03 = read.csv("SyntenicGenes/20250314_A03Only_synTableTandemRm.csv")

# Clean up 
A03_At = A03 %>% separate_rows(ATG, sep = ", ") %>%
  # append an underscore and a sequential number after the A03 annotation to keep track of multi-mappings
  group_by(ATG) %>%
  mutate(Gene = case_when(
    n() > 1 ~ paste0(Gene, "_", row_number()),  # Add iterative numbers when there are multiple Athaliana genes
    TRUE ~ Gene)) %>% ungroup()

# Split off the copy number ID that was just created into its own column
A03_At = A03_At %>% separate(Gene, c("Gene", "A03_CopyNumID"), sep = "_")  


## How many A03 TFs do we find?
  # Bring in the newer larger list of TFs (2296 Arabidopsis TFs vs 760)
atgs = read.csv("GeneRegulatoryNetworks/Ath_TF_list.txt")
colnames(atgs) = "Gene"
# Separate
atgs = atgs %>% separate(Gene, c("isoform", "Gene", "CommonName"), sep = "\t")
# Get rid of the isoform information for now (just take the first one)
atgs = atgs[-1]
atgs = unique(atgs)
  # Which clock genes are already in the list? 
clock[which(clock$Gene %in% atgs$Gene),]

# Add any clock genes that are missing from the list 
AthTFs = rbind(atgs, clock)

# Sanity Check
unique(AthTFs$Gene) # This tells us there were already 22 clock genes in the original list

# Pull unique once more to make sure there weren't some clock genes that were already in the list
AthTFs = AthTFs[which(!(duplicated(AthTFs$Gene))), ]
#write.csv("SyntenicGenes/AthTFs_withClock.csv", row.names = FALSE)

# This is for Dennis Such's masters thesis but good to have regardless 
A03_AddedClockGenes = clock[which(clock$Gene %in% A03_At$ATG), ]
write.csv(A03_AddedClockGenes, file = "GeneRegulatoryNetworks/A03_AddedClockGenes.csv", row.names = FALSE)

## Join the two dfs together by subsetting for Arabidopsis TFs
BrTFs = left_join(AthTFs, A03_At, by = join_by(Gene == ATG), relationship = "many-to-many")
  # Remove any rows that do not have a Bra match
BrTFs = BrTFs[!is.na(BrTFs$Gene.y),] %>% pull(Gene.y)
BrTFs = unique(BrTFs)


## Write out the A03 TFs that will be used to build the GRNs for each reference line (since all annotations will be anchored to A03)
write.csv(BrTFs, file = "SyntenicGenes/20250218_BrTFs_A03v2_PlusClock.csv", row.names = FALSE)
##### ##### ##### ##### ##### ##### ##### ##### 


##### ##### ##### ##### ##### ##### ##### ##### 
## Bring in Chiemeka's 'GRN_ready' tables
R500 = read.csv("SyntenicGenes/20250314_R500Only_synTableTandemRm.csv")

# Clean up 
R500_At = R500 %>% separate_rows(ATG, sep = ", ") %>%
  # append an underscore and a sequential number after the R500 annotation to keep track of multi-mappings
  group_by(ATG) %>%
  mutate(Gene = case_when(
    n() > 1 ~ paste0(Gene, "_", row_number()),  # Add iterative numbers when there are multiple Athaliana genes
    TRUE ~ Gene)) %>% ungroup()

# Split off the copy number ID that was just created into its own column
R500_At = R500_At %>% separate(Gene, c("Gene", "R500_CopyNumID"), sep = "_")  


## How many R500 TFs do we find?
# Bring in the newer larger list of TFs (2296 Arabidopsis TFs vs 760)
atgs = read.csv("GeneRegulatoryNetworks/Ath_TF_list.txt")
colnames(atgs) = "Gene"
# Separate
atgs = atgs %>% separate(Gene, c("isoform", "Gene", "CommonName"), sep = "\t")
# Get rid of the isoform information for now (just take the first one)
atgs = atgs[-1]
atgs = unique(atgs)
# Which clock genes are already in the list? 
clock[which(clock$Gene %in% atgs$Gene),]

# Add any clock genes that are missing from the list 
AthTFs = rbind(atgs, clock)

# Sanity Check
unique(AthTFs$Gene) # This tells us there were already 22 clock genes in the original list

# Pull unique once more to make sure there weren't some clock genes that were already in the list
AthTFs = AthTFs[which(!(duplicated(AthTFs$Gene))), ]
#write.csv("SyntenicGenes/AthTFs_withClock.csv", row.names = FALSE)

# This is for Dennis Such's masters thesis but good to have regardless 
R500_AddedClockGenes = clock[which(clock$Gene %in% R500_At$ATG), ]
write.csv(R500_AddedClockGenes, file = "GeneRegulatoryNetworks/R500_AddedClockGenes.csv", row.names = FALSE)


##### ##### ##### ##### ##### ##### ##### ##### 

##### ##### ##### ##### ##### ##### ##### ##### 
## L58
clock = read.csv("SyntenicGenes/UniqueClockATGs.csv")
AthTFs = read.csv("SyntenicGenes/AthTFs_withClock.csv")


## Bring in Chiemeka's 'GRN_ready' tables
L58 = read.csv("SyntenicGenes/GRN_ReadyTables_FromChiemeka/L58_Syntenic_Table_GRN_Ready.csv")

# Clean up 
L58_At = L58 %>% separate_rows(ATG, sep = ", ") %>%
  # append an underscore and a sequential number after the L58 annotation to keep track of multi-mappings
  group_by(ATG) %>%
  mutate(Gene = case_when(
    n() > 1 ~ paste0(Gene, "_", row_number()),  # Add iterative numbers when there are multiple Athaliana genes
    TRUE ~ Gene)) %>% ungroup()

# Split off the copy number ID that was just created into its own column
L58_At = L58_At %>% separate(Gene, c("Gene", "L58_CopyNumID"), sep = "_")  

## Join the two dfs together by subsetting for Arabidopsis TFs
BrTFs = left_join(AthTFs, L58_At, by = join_by(Gene == ATG), relationship = "many-to-many")
# Remove any rows that do not have a Bra match
BrTFs = BrTFs[!is.na(BrTFs$Gene.y),] %>% pull(Gene.y)
BrTFs = unique(BrTFs)


## Write out the L58 TFs that will be used to build the GRNs for each reference line (since all annotations will be anchored to L58)
write.csv(BrTFs, file = "SyntenicGenes/20250221_BrTFs_L58v1_PlusClock.csv", row.names = FALSE)
##### ##### ##### ##### ##### ##### ##### ##### 


##### ##### ##### ##### ##### ##### ##### ##### 
## Repeat with additional lines 
# Bring in these files if starting fresh
#clock = read.csv("SyntenicGenes/UniqueClockATGs.csv")
#AthTFs = read.csv("SyntenicGenes/AthTFs_withClock.csv")


## Bring in Chiemeka's 'GRN_ready' tables
PCGlu = read.csv("SyntenicGenes/GRN_ReadyTables_FromChiemeka/PCGlu_Syntenic_Table_GRN_Ready.csv")

# Clean up 
PCGlu_At = PCGlu %>% separate_rows(ATG, sep = ", ") %>%
  # append an underscore and a sequential number after the PCGlu annotation to keep track of multi-mappings
  group_by(ATG) %>%
  mutate(Gene = case_when(
    n() > 1 ~ paste0(Gene, "_", row_number()),  # Add iterative numbers when there are multiple Athaliana genes
    TRUE ~ Gene)) %>% ungroup()

# Split off the copy number ID that was just created into its own column
PCGlu_At = PCGlu_At %>% separate(Gene, c("Gene", "PCGlu_CopyNumID"), sep = "_")  

## Join the two dfs together by subsetting for Arabidopsis TFs
BrTFs = left_join(AthTFs, PCGlu_At, by = join_by(Gene == ATG), relationship = "many-to-many")
# Remove any rows that do not have a Bra match
BrTFs = BrTFs[!is.na(BrTFs$Gene.y),] %>% pull(Gene.y)
BrTFs = unique(BrTFs)


## Write out the PCGlu TFs that will be used to build the GRNs for each reference line (since all annotations will be anchored to PCGlu)
write.csv(BrTFs, file = "SyntenicGenes/20250221_BrTFs_PCGluv1_PlusClock.csv", row.names = FALSE)
##### ##### ##### ##### ##### ##### ##### ##### 


##### ##### ##### ##### ##### ##### ##### ##### 
## Repeat with additional lines 
# Bring in these files if starting fresh
#clock = read.csv("SyntenicGenes/UniqueClockATGs.csv")
#AthTFs = read.csv("SyntenicGenes/AthTFs_withClock.csv")


## Bring in Chiemeka's 'GRN_ready' tables
VT123 = read.csv("SyntenicGenes/GRN_ReadyTables_FromChiemeka/VT123_Syntenic_Table_GRN_Ready.csv")

# Clean up 
VT123_At = VT123 %>% separate_rows(ATG, sep = ", ") %>%
  # append an underscore and a sequential number after the VT123 annotation to keep track of multi-mappings
  group_by(ATG) %>%
  mutate(Gene = case_when(
    n() > 1 ~ paste0(Gene, "_", row_number()),  # Add iterative numbers when there are multiple Athaliana genes
    TRUE ~ Gene)) %>% ungroup()

# Split off the copy number ID that was just created into its own column
VT123_At = VT123_At %>% separate(Gene, c("Gene", "VT123_CopyNumID"), sep = "_")  

## Join the two dfs together by subsetting for Arabidopsis TFs
BrTFs = left_join(AthTFs, VT123_At, by = join_by(Gene == ATG), relationship = "many-to-many")
# Remove any rows that do not have a Bra match
BrTFs = BrTFs[!is.na(BrTFs$Gene.y),] %>% pull(Gene.y)
BrTFs = unique(BrTFs)


## Write out the VT123 TFs that will be used to build the GRNs for each reference line (since all annotations will be anchored to VT123)
write.csv(BrTFs, file = "SyntenicGenes/20250221_BrTFs_VT123v1_PlusClock.csv", row.names = FALSE)
##### ##### ##### ##### ##### ##### ##### ##### 


##### ##### ##### ##### ##### ##### ##### ##### 
## Repeat with additional lines 
# Bring in these files if starting fresh
clock = read.csv("SyntenicGenes/UniqueClockATGs.csv")
AthTFs = read.csv("SyntenicGenes/AthTFs_withClock.csv")


## Bring in Chiemeka's 'GRN_ready' tables
WO83 = read.csv("SyntenicGenes/GRN_ReadyTables_FromChiemeka/WO83_Syntenic_Table_GRN_Ready.csv")

# Clean up 
## WO83 specific!! 
WO83$Gene = gsub("WO_83", "WO83", WO83$Gene)

WO83_At = WO83 %>% separate_rows(ATG, sep = ", ") %>%
  # append an underscore and a sequential number after the WO83 annotation to keep track of multi-mappings
  group_by(ATG) %>%
  mutate(Gene = case_when(
    n() > 1 ~ paste0(Gene, "_", row_number()),  # Add iterative numbers when there are multiple Athaliana genes
    TRUE ~ Gene)) %>% ungroup()

# Split off the copy number ID that was just created into its own column
WO83_At = WO83_At %>% separate(Gene, c("Gene", "WO83_CopyNumID"), sep = "_")  

## Join the two dfs together by subsetting for Arabidopsis TFs
BrTFs = left_join(AthTFs, WO83_At, by = join_by(Gene == ATG), relationship = "many-to-many")
# Remove any rows that do not have a Bra match
BrTFs = BrTFs[!is.na(BrTFs$Gene.y),] %>% pull(Gene.y)
BrTFs = unique(BrTFs)


## Write out the WO83 TFs that will be used to build the GRNs for each reference line (since all annotations will be anchored to WO83)
write.csv(BrTFs, file = "SyntenicGenes/20250219_BrTFs_WO83v1_PlusClock.csv", row.names = FALSE)
##### ##### ##### ##### ##### ##### ##### ##### 

