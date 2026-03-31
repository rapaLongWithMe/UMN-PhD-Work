### Associating CNSs with a particular expression pattern 
  
## Set up the environment
setwd("~/Desktop/JGI_Cold/Paralogs")
require(tidyverse)

## Workflow
  # Step 1: Make CNS table that has a column for: 
    # Every possible CNS across all six ref lines
    # A column of Y/N that indicates whether or not the first paralog is associated with that CNS
    # A column of Y/N that indicates whether or not the next copy of that paralog is associated with that CNS
    # A column of Y/N that indicates whether or not the third paralog, if present in three copies, is associated with that CNS

  # Step 2: Filter
    # We want to retain rows (CNSs) that have a mix of Y/N
    # So filter out rows with all Y or all N 

  # Step 3: Bring in additional information for further filtering
    # Bring in cluster assignments to look for CNSs that are good candidates for a particular cold response
    # Compare (phase?) of control expression to identify CNSs with specific diel expression patterns


##### ##### ##### ##### ##### ##### ##### ##### 
##### Step 1

# Bring in the cns data - want all results for all ref lines (no filtering)
cns1 = read.csv("~/Desktop/Git/BrapaPangenome_ColdAcclimation/CNS_BLAST/Brapassp_chinensisvar_communisPCGluCNS_BLAST_results.csv")
cns1 = cns1[, c(2,12)]
cns1$Nearest_Gene_ID = gsub("\\.v[0-9].[0-9]", "", cns1$Nearest_Gene_ID)
# Second geno
cns2 = read.csv("~/Desktop/Git/BrapaPangenome_ColdAcclimation/CNS_BLAST/Brapassp_chinensisvar_parachinensisL58CNS_BLAST_results.csv")
cns2 = cns2[, c(2,12)]
cns2$Nearest_Gene_ID = gsub("\\.v[0-9].[0-9]", "", cns2$Nearest_Gene_ID)
# Third geno
cns3 = read.csv("~/Desktop/Git/BrapaPangenome_ColdAcclimation/CNS_BLAST/Brapassp_oleiferavar_oleiferaWO_83CNS_BLAST_results.csv")
cns3 = cns3[, c(2,12)]
cns3$Nearest_Gene_ID = gsub("\\.v[0-9].[0-9]", "", cns3$Nearest_Gene_ID)
cns3$Nearest_Gene_ID = gsub("WO_83", "WO83", cns3$Nearest_Gene_ID)
# Fourth geno
cns4 = read.csv("~/Desktop/Git/BrapaPangenome_ColdAcclimation/CNS_BLAST/Brapassp_pekinensisvar_pekinensisA03CNS_BLAST_results.csv")
cns4 = cns4[, c(2,12)]
cns4$Nearest_Gene_ID = gsub("\\.v[0-9].[0-9]", "", cns4$Nearest_Gene_ID)
# Fifth geno
cns5 = read.csv("~/Desktop/Git/BrapaPangenome_ColdAcclimation/CNS_BLAST/Brapassp_rapavar_rapaVT12CNS_BLAST_results.csv")
cns5 = cns5[, c(2,12)]
cns5$Nearest_Gene_ID = gsub("\\.v[0-9].[0-9]", "", cns5$Nearest_Gene_ID)
# Sixth geno
cns6 = read.csv("~/Desktop/Git/BrapaPangenome_ColdAcclimation/CNS_BLAST/Brapassp_trilocularisR500CNS_BLAST_results.csv")
cns6 = cns6[, c(2,12)]
cns6$Nearest_Gene_ID = gsub("\\.v[0-9].[0-9]", "", cns6$Nearest_Gene_ID)


## Combine the CNS results from all ref lines
cns_results = rbind(cns1, cns2,cns3,cns4,cns5,cns6)

# Check to make sure that if a gene is found in multiple rows it is associated with different CNSs (no duplicates)
which(duplicated(cns_results))
cns_results[900,] # Use this to search the full table for this gene
  # Looked at several of them, they seem like real dups, not sure how they got in there 
  
# Remove the 84 random duplicates
cns_results = unique(cns_results)
# Clean up the environment
#rm(cns1, cns2, cns3, cns4, cns5, cns6)



##### SHOULD MAKE SURE I HAVEN'T DONE THIS PREVIOUSLY!!! #####
  # Read in the two- and three-copy gene tables
twoCopies <- read.csv("FinalSynTable_twoCopies_withATG_AndSingletons_AndSig.csv")
threeCopies <- read.csv("FinalSynTable_threeCopies_withATG_AndSingletons_AndSig.csv")

# Pivot two-copy table to long format
two_long <- twoCopies %>%
  pivot_longer(cols = -ATG, names_to = "Copy_Label", values_to = "Gene") %>%
  filter(!is.na(Gene)) %>%
  mutate(
    CopyNum = ifelse(str_detect(Copy_Label, "Bra_1"), "Copy1", "Copy2"),
    CopySetSize = 2
  ) %>%
  select(ATG, CopyNum, Gene, CopySetSize)

# Pivot three-copy table to long format
three_long <- threeCopies %>%
  pivot_longer(cols = -ATG, names_to = "Copy_Label", values_to = "Gene") %>%
  filter(!is.na(Gene)) %>%
  mutate(
    CopyNum = case_when(
      str_detect(Copy_Label, "Bra_1") ~ "Copy1",
      str_detect(Copy_Label, "Bra_2") ~ "Copy2",
      str_detect(Copy_Label, "Bra_3") ~ "Copy3"
    ),
    CopySetSize = 3
  ) %>%
  select(ATG, CopyNum, Gene, CopySetSize)

# Combine both tables
combined_long <- bind_rows(two_long, three_long) %>%
  # Extract genotype (remove Br/Bra prefix, then take everything up to the first dot)
  mutate(
    Genotype = str_remove(Gene, "^Bra|^Br"),
    Genotype = str_extract(Genotype, "^[^\\.]+")
  )


# Get the duplicated gene rows
combined_long %>%
  filter(duplicated(Gene) | duplicated(Gene, fromLast = TRUE)) %>%
  arrange(Gene)

# Check 
levels(factor(combined_long$CopySetSize))
levels(factor(combined_long$Genotype))
which(duplicated(combined_long$Gene))

# Get the duplicated gene rows
View(combined_long %>%
  filter(duplicated(Gene) | duplicated(Gene, fromLast = TRUE)) %>%
  arrange(Gene))
# Not sure how these ended up in both 2 and 3 copy sets
  # Write out to check on later and then remove for now
Why_dups = combined_long %>%
       filter(duplicated(Gene) | duplicated(Gene, fromLast = TRUE)) %>%
       arrange(Gene)
write.csv(Why_dups, file = "20250416_WhyAreTheseInTwoAndThreeCopies.csv", row.names = FALSE)
  # Since I don't know for sure if these are two or three copy, remove from both sets
Allparas = combined_long %>%
  filter(!duplicated(Gene) & !duplicated(Gene, fromLast = TRUE)) # Should remove 292*2 rows

# Write out this table for downstream analyses 
write.csv(Allparas, file = "20250416_SigAndNonSig_AllParas_TwoAndThreeCopies.csv", row.names = FALSE)


## Start building the table
cns_table_wide <- cns_results %>%
  # 1. Extract genotype from gene names
  mutate(Genotype = str_remove(Nearest_Gene_ID, "^Bra|^Br"),
         Genotype = str_extract(Genotype, "^[^\\.]+")) %>%
  # 2. Add an index per CNS_ID and Genotype
  group_by(CNS_ID, Genotype) %>%
  arrange(CNS_ID, Genotype, Nearest_Gene_ID) %>%
  mutate(Gene_Index = row_number()) %>%
  ungroup() %>%
  # 3. Create composite column names like Gene_R500_1, Gene_L58_2
  mutate(Gene_Col = paste0("Gene_", Genotype, "_", Gene_Index)) %>%
  select(CNS_ID, Gene_Col, Nearest_Gene_ID) %>%
  pivot_wider(names_from = Gene_Col, values_from = Nearest_Gene_ID)


## Use Allparas to filter out paralogs with high copy number (just interested in these anyway)
valid_genes <- Allparas$Gene
which(duplicated(valid_genes))

# Get gene columns (excluding CNS_ID)
gene_cols <- names(cns_table_wide)[names(cns_table_wide) != "CNS_ID"]

# Convert gene columns to binary presence/absence for each paralog
  # Double check that the highcopy number paralogs are all zeros; should be removed in next step
binary_matrix <- cns_table_wide %>%
  mutate(across(all_of(gene_cols), ~ ifelse(. %in% valid_genes, 1, 0)))


## Remove columns that are now all 0 (no valid genes)
  # Need to do this by genotype so split into a list 
  # Get genotype names from binary column names
genotypes <- colnames(binary_matrix) %>%
  grep("^Gene_", ., value = TRUE) %>%
  str_extract("(?<=Gene_)[^_]+") %>%
  unique()

#  Extract relevant columns (+ CNS_ID) into a list
binary_by_genotype <- setNames(
  lapply(genotypes, function(gt) {
    cols <- grep(paste0("Gene_", gt, "_"), colnames(binary_matrix), value = TRUE)
    binary_matrix %>%
      select(CNS_ID, all_of(cols))
  }),
  genotypes
)

## Add a new column that sums the number of hits per CNS
binary_by_genotype <- lapply(names(binary_by_genotype), function(gt) {
  df <- binary_by_genotype[[gt]]
  df %>%
    mutate(HitCount = rowSums(select(., -CNS_ID)))
}) %>%
  setNames(names(binary_by_genotype))





##### UPDATE THIS TO WORK ON A LIST: #####
## Check that there are no more than 5 columns per line (2copy + 3 copy)
a03_cols <- grep("A03", colnames(binary_matrix), value = TRUE)

binary_matrix %>%
  rowwise() %>%
  mutate(A03_total = sum(c_across(all_of(a03_cols)))) %>%
  ungroup() %>%
  filter(A03_total > 5)


## Check that the genes in the table are only from the list of paralogs (multicopies)
# Step 1: Extract all gene values from cns_table_wide
all_gene_values <- unlist(cns_table_wide[, -1])  # remove CNS_ID

# Step 2: Check which genes are not in Allparas
setdiff(na.omit(unique(all_gene_values)), Allparas$Gene)  # There are still 7296 ...need to do the filtering by genotype 

