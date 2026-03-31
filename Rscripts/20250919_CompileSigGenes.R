### Compiling kMEs and kMeds from the DiPALM permutation pipeline 

## Set up the environment
setwd("~/Desktop/JGI_Cold/DiPALMperm_output/")

require(tidyverse)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### COMPILE DIPALM GENES (kME + kMed)

## Pull kMeds
## We had to rerun the perm pipeline after adding in code to call kMed genes
  # pull those out now (and ignore the kMEs that came with that later run since we already have them)

load("A03/kMeds/dipalm.RData") 
  # Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "A03/A03_kMedOnly_PermOutput.RData")
A03_kMed = sig.MedGenes
  # Add a column to indicate which type of CR gene it is
A03_kMed = data.frame(Gene = A03_kMed, CR_type = "kMed")
  # Remove the isoform information
A03_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", A03_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("ACC28/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "ACC28/ACC28_kMedOnly_PermOutput.RData")
ACC28_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
ACC28_kMed = data.frame(Gene = ACC28_kMed, CR_type = "kMed")
# Remove the isoform information
ACC28_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", ACC28_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("ACC50/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "ACC50/ACC50_kMedOnly_PermOutput.RData")
ACC50_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
ACC50_kMed = data.frame(Gene = ACC50_kMed, CR_type = "kMed")
# Remove the isoform information
ACC50_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", ACC50_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("CC168/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "CC168/CC168_kMedOnly_PermOutput.RData")
CC168_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
CC168_kMed = data.frame(Gene = CC168_kMed, CR_type = "kMed")
# Remove the isoform information
CC168_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", CC168_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

## FT005 was only run once so data isn't split into folders by kME or kMed (happened together)
load("FT005/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "FT005/FT005_kMedOnly_PermOutput.RData")
FT005_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
FT005_kMed = data.frame(Gene = FT005_kMed, CR_type = "kMed")
# Remove the isoform information
FT005_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", FT005_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("HN53/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "HN53/HN53_kMedOnly_PermOutput.RData")
HN53_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
HN53_kMed = data.frame(Gene = HN53_kMed, CR_type = "kMed")
# Remove the isoform information
HN53_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", HN53_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("L58/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "L58/L58_kMedOnly_PermOutput.RData")
L58_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
L58_kMed = data.frame(Gene = L58_kMed, CR_type = "kMed")
# Remove the isoform information
L58_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", L58_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)


load("PC185/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "PC185/PC185_kMedOnly_PermOutput.RData")
PC185_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
PC185_kMed = data.frame(Gene = PC185_kMed, CR_type = "kMed")
# Remove the isoform information
PC185_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", PC185_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)


load("PCGlu/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "PCGlu/PCGlu_kMedOnly_PermOutput.RData")
PCGlu_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
PCGlu_kMed = data.frame(Gene = PCGlu_kMed, CR_type = "kMed")
# Remove the isoform information
PCGlu_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", PCGlu_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("R500/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "R500/R500_kMedOnly_PermOutput.RData")
R500_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
R500_kMed = data.frame(Gene = R500_kMed, CR_type = "kMed")
# Remove the isoform information
R500_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", R500_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("VT123/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "VT123/VT123_kMedOnly_PermOutput.RData")
VT123_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
VT123_kMed = data.frame(Gene = VT123_kMed, CR_type = "kMed")
# Remove the isoform information
VT123_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", VT123_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("WO83/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "WO83/WO83_kMedOnly_PermOutput.RData")
WO83_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
WO83_kMed = data.frame(Gene = WO83_kMed, CR_type = "kMed")
# Remove the isoform information
WO83_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", WO83_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("ZCT/kMeds/dipalm.RData") 
# Save just the kMed data
#save(AdjMed, sig.MedGenes, file = "ZCT/ZCT_kMedOnly_PermOutput.RData")
ZCT_kMed = sig.MedGenes
# Add a column to indicate which type of CR gene it is
ZCT_kMed = data.frame(Gene = ZCT_kMed, CR_type = "kMed")
# Remove the isoform information
ZCT_kMed$Gene = gsub("\\.v[0-9].[0-9]", "", ZCT_kMed$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)


##### ##### ##### ##### ##### ##### ##### 
## Compile into a single list 
  # Load if needed 
load("DiPALMperm_Output/A03/A03_kMedOnly_PermOutput.RData")
A03_kMed = sig.MedGenes
load("DiPALMperm_Output/HN53/HN53_kMedOnly_PermOutput.RData")
HN53_kMed = sig.MedGenes
load("DiPALMperm_Output/CC168/CC168_kMedOnly_PermOutput.RData")
CC168_kMed = sig.MedGenes
load("DiPALMperm_Output/R500/R500_kMedOnly_PermOutput.RData")
R500_kMed = sig.MedGenes
load("DiPALMperm_Output/ACC50/ACC50_kMedOnly_PermOutput.RData")
ACC50_kMed = sig.MedGenes
load("DiPALMperm_Output/ACC28/ACC28_kMedOnly_PermOutput.RData")
ACC28_kMed = sig.MedGenes
load("DiPALMperm_Output/L58/L58_kMedOnly_PermOutput.RData")
L58_kMed = sig.MedGenes
load("DiPALMperm_Output/PC185/PC185_kMedOnly_PermOutput.RData")
PC185_kMed = sig.MedGenes
load("DiPALMperm_Output/VT123/VT123_kMedOnly_PermOutput.RData")
VT123_kMed = sig.MedGenes
load("DiPALMperm_Output/FT005/FT005_kMedOnly_PermOutput.RData")
FT005_kMed = sig.MedGenes
load("DiPALMperm_Output/PCGlu/PCGlu_kMedOnly_PermOutput.RData")
PCGlu_kMed = sig.MedGenes
load("DiPALMperm_Output/ZCT/ZCT_kMedOnly_PermOutput.RData")
ZCT_kMed = sig.MedGenes
load("DiPALMperm_Output/WO83/WO83_kMedOnly_PermOutput.RData")
WO83_kMed = sig.MedGenes
##### ##### ##### ##### ##### ##### ##### 


##### ##### ##### ##### #####
## Bring in the kMEs
load("A03/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "A03/A03_kMEOnly_PermOutput.RData")
A03_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
A03_kMEs = data.frame(Gene = A03_kMEs, CR_type = "kME")
# Remove the isoform information
A03_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", A03_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("ACC28/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "ACC28/ACC28_kMEOnly_PermOutput.RData")
ACC28_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
ACC28_kMEs = data.frame(Gene = ACC28_kMEs, CR_type = "kME")
# Remove the isoform information
ACC28_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", ACC28_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("ACC50/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "ACC50/ACC50_kMEOnly_PermOutput.RData")
ACC50_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
ACC50_kMEs = data.frame(Gene = ACC50_kMEs, CR_type = "kME")
# Remove the isoform information
ACC50_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", ACC50_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("CC168/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "CC168/CC168_kMEOnly_PermOutput.RData")
CC168_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
CC168_kMEs = data.frame(Gene = CC168_kMEs, CR_type = "kME")
# Remove the isoform information
CC168_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", CC168_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("FT005/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "FT005/FT005_kMEOnly_PermOutput.RData")
FT005_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
FT005_kMEs = data.frame(Gene = FT005_kMEs, CR_type = "kME")
# Remove the isoform information
FT005_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", FT005_kMEs$Gene)
rm(AdjkMEs, AdjMed, patternCor, sig.genes, sig.MedGenes)

load("HN53/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "HN53/HN53_kMEOnly_PermOutput.RData")
HN53_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
HN53_kMEs = data.frame(Gene = HN53_kMEs, CR_type = "kME")
# Remove the isoform information
HN53_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", HN53_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("L58/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "L58/L58_kMEOnly_PermOutput.RData")
L58_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
L58_kMEs = data.frame(Gene = L58_kMEs, CR_type = "kME")
# Remove the isoform information
L58_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", L58_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("PCGlu/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "PCGlu/PCGlu_kMEOnly_PermOutput.RData")
PCGlu_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
PCGlu_kMEs = data.frame(Gene = PCGlu_kMEs, CR_type = "kME")
# Remove the isoform information
PCGlu_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", PCGlu_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("R500/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "R500/R500_kMEOnly_PermOutput.RData")
R500_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
R500_kMEs = data.frame(Gene = R500_kMEs, CR_type = "kME")
# Remove the isoform information
R500_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", R500_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("VT123/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "VT123/VT123_kMEOnly_PermOutput.RData")
VT123_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
VT123_kMEs = data.frame(Gene = VT123_kMEs, CR_type = "kME")
# Remove the isoform information
VT123_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", VT123_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("WO83/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "WO83/WO83_kMEOnly_PermOutput.RData")
WO83_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
WO83_kMEs = data.frame(Gene = WO83_kMEs, CR_type = "kME")
# Remove the isoform information
WO83_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", WO83_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("ZCT/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "ZCT/ZCT_kMEOnly_PermOutput.RData")
ZCT_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
ZCT_kMEs = data.frame(Gene = ZCT_kMEs, CR_type = "kME")
# Remove the isoform information
ZCT_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", ZCT_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)

load("PC185/kMEs/dipalm.RData") 
# Save just the Sig kMEs
save(AdjkMEs, sig.genes, file = "PC185/PC185_kMEOnly_PermOutput.RData")
PC185_kMEs = sig.genes
# Add a column to indicate which type of CR gene it is
PC185_kMEs = data.frame(Gene = PC185_kMEs, CR_type = "kME")
# Remove the isoform information
PC185_kMEs$Gene = gsub("\\.v[0-9].[0-9]", "", PC185_kMEs$Gene)
rm(AdjkMEs, patternCor, sig.genes)


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### PULL UNIQUE KME AND KME (UNION)
## Pull the union of kME and kMed (which are unique to both lists)
A03_UnionkME_kMeds = union(A03_kMEs$Gene, A03_kMed$Gene)
#write.csv(A03_UnionkME_kMeds, file = "A03/A03_UnionkME_kMeds.csv", row.names = FALSE)
  # Pull the genes that are unique to kME
A03_kME_only  <- setdiff(A03_kMEs$Gene, A03_kMed$Gene)
  # Pull the genes that are unique to kME
A03_kMed_only <- setdiff(A03_kMed$Gene, A03_kMEs$Gene)
  # Add a CR_type column
A03_kME_only = data.frame(Gene = A03_kME_only, CR_type = "kME")
A03_kMed_only = data.frame(Gene = A03_kMed_only, CR_type = "kMed")
  # Check that the sum of these two and the intersection add up
A03_shared = intersect(A03_kMEs$Gene, A03_kMed$Gene)
  # Save the lists for later
write.csv(A03_kME_only, file = "A03/A03_kMEsOnly.csv", row.names = FALSE)
write.csv(A03_kMed_only, file = "A03/A03_kMedOnly.csv", row.names = FALSE)
  # Can keep this as a vector since it includes both kME and kMed genes
write.csv(A03_shared, file = "A03/A03_sharedkMEkMed.csv", row.names = FALSE)
rm(A03_kME_only, A03_kMed_only, A03_kMEs, A03_kMed, A03_shared)


## Pull the union of kME and kMed (which are unique to both lists)
ACC28_UnionkME_kMeds = union(ACC28_kMEs$Gene, ACC28_kMed$Gene)
#write.csv(ACC28_UnionkME_kMeds, file = "ACC28/ACC28_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
ACC28_kME_only  <- setdiff(ACC28_kMEs$Gene, ACC28_kMed$Gene)
# Pull the genes that are unique to kME
ACC28_kMed_only <- setdiff(ACC28_kMed$Gene, ACC28_kMEs$Gene)
# Add a CR_type column
ACC28_kME_only = data.frame(Gene = ACC28_kME_only, CR_type = "kME")
ACC28_kMed_only = data.frame(Gene = ACC28_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
ACC28_shared = intersect(ACC28_kMEs$Gene, ACC28_kMed$Gene)
# Save the lists for later
write.csv(ACC28_kME_only, file = "ACC28/ACC28_kMEsOnly.csv", row.names = FALSE)
write.csv(ACC28_kMed_only, file = "ACC28/ACC28_kMedOnly.csv", row.names = FALSE)
write.csv(ACC28_shared, file = "ACC28/ACC28_sharedkMEkMed.csv", row.names = FALSE)
rm(ACC28_kME_only, ACC28_kMed_only, ACC28_kMEs, ACC28_kMed, ACC28_shared)


ACC50_UnionkME_kMeds = union(ACC50_kMEs$Gene, ACC50_kMed$Gene)
write.csv(ACC50_UnionkME_kMeds, file = "ACC50/ACC50_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
ACC50_kME_only  <- setdiff(ACC50_kMEs$Gene, ACC50_kMed$Gene)
# Pull the genes that are unique to kME
ACC50_kMed_only <- setdiff(ACC50_kMed$Gene, ACC50_kMEs$Gene)
# Add a CR_type column
ACC50_kME_only = data.frame(Gene = ACC50_kME_only, CR_type = "kME")
ACC50_kMed_only = data.frame(Gene = ACC50_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
ACC50_shared = intersect(ACC50_kMEs$Gene, ACC50_kMed$Gene)
# Save the lists for later
write.csv(ACC50_kME_only, file = "ACC50/ACC50_kMEsOnly.csv", row.names = FALSE)
write.csv(ACC50_kMed_only, file = "ACC50/ACC50_kMedOnly.csv", row.names = FALSE)
write.csv(ACC50_shared, file = "ACC50/ACC50_sharedkMEkMed.csv", row.names = FALSE)
rm(ACC50_kME_only, ACC50_kMed_only, ACC50_kMEs, ACC50_kMed, ACC50_shared)


CC168_UnionkME_kMeds = union(CC168_kMEs$Gene, CC168_kMed$Gene)
#write.csv(CC168_UnionkME_kMeds, file = "CC168/CC168_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
CC168_kME_only  <- setdiff(CC168_kMEs$Gene, CC168_kMed$Gene)
# Pull the genes that are unique to kME
CC168_kMed_only <- setdiff(CC168_kMed$Gene, CC168_kMEs$Gene)
# Add a CR_type column
CC168_kME_only = data.frame(Gene = CC168_kME_only, CR_type = "kME")
CC168_kMed_only = data.frame(Gene = CC168_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
CC168_shared = intersect(CC168_kMEs$Gene, CC168_kMed$Gene)
# Save the lists for later
write.csv(CC168_kME_only, file = "CC168/CC168_kMEsOnly.csv", row.names = FALSE)
write.csv(CC168_kMed_only, file = "CC168/CC168_kMedOnly.csv", row.names = FALSE)
write.csv(CC168_shared, file = "CC168/CC168_sharedkMEkMed.csv", row.names = FALSE)
rm(CC168_kME_only, CC168_kMed_only, CC168_kMEs, CC168_kMed, CC168_shared)


## Had to update this since FT005 is pulling from one dipalm object not two like the others (run together)
FT005_UnionkME_kMeds = union(FT005_kMEs$Gene, FT005_kMed$Gene)
#write.csv(FT005_UnionkME_kMeds, file = "DiPALMperm_output/FT005/FT005_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
FT005_kME_only  <- setdiff(FT005_kMEs$Gene, FT005_kMed$Gene)
# Pull the genes that are unique to kMED
FT005_kMed_only <- setdiff(FT005_kMed$Gene, FT005_kMEs$Gene)
# Add a CR_type column
FT005_kME_only = data.frame(Gene = FT005_kME_only, CR_type = "kME")
FT005_kMed_only = data.frame(Gene = FT005_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
FT005_shared = intersect(FT005_kMEs$Gene, FT005_kMed$Gene)
# Save the lists for later
write.csv(FT005_kME_only, file = "FT005/FT005_kMEsOnly.csv", row.names = FALSE)
write.csv(FT005_kMed_only, file = "FT005/FT005_kMedOnly.csv", row.names = FALSE)
write.csv(FT005_shared, file = "FT005/FT005_sharedkMEkMed.csv", row.names = FALSE)
rm(FT005_kME_only, FT005_kMed_only, FT005_kMEs, FT005_kMed, FT005_shared)


HN53_UnionkME_kMeds = union(HN53_kMEs$Gene, HN53_kMed$Gene)
#write.csv(HN53_UnionkME_kMeds, file = "HN53/HN53_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
HN53_kME_only  <- setdiff(HN53_kMEs$Gene, HN53_kMed$Gene)
# Pull the genes that are unique to kMed
HN53_kMed_only <- setdiff(HN53_kMed$Gene, HN53_kMEs$Gene)
# Add a CR_type column
HN53_kME_only = data.frame(Gene = HN53_kME_only, CR_type = "kME")
HN53_kMed_only = data.frame(Gene = HN53_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
HN53_shared = intersect(HN53_kMEs$Gene, HN53_kMed$Gene)
# Save the lists for later
write.csv(HN53_kME_only, file = "HN53/HN53_kMEsOnly.csv", row.names = FALSE)
write.csv(HN53_kMed_only, file = "HN53/HN53_kMedOnly.csv", row.names = FALSE)
write.csv(HN53_shared, file = "HN53/HN53_sharedkMEkMed.csv", row.names = FALSE)
rm(HN53_kME_only, HN53_kMed_only, HN53_kMEs, HN53_kMed, HN53_shared)


L58_UnionkME_kMeds = union(L58_kMEs$Gene, L58_kMed$Gene)
#write.csv(L58_UnionkME_kMeds, file = "L58/L58_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
L58_kME_only  <- setdiff(L58_kMEs$Gene, L58_kMed$Gene)
# Pull the genes that are unique to kMed
L58_kMed_only <- setdiff(L58_kMed$Gene, L58_kMEs$Gene)
# Add a CR_type column
L58_kME_only = data.frame(Gene = L58_kME_only, CR_type = "kME")
L58_kMed_only = data.frame(Gene = L58_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
L58_shared = intersect(L58_kMEs$Gene, L58_kMed$Gene)
# Save the lists for later
write.csv(L58_kME_only, file = "L58/L58_kMEsOnly.csv", row.names = FALSE)
write.csv(L58_kMed_only, file = "L58/L58_kMedOnly.csv", row.names = FALSE)
write.csv(L58_shared, file = "L58/L58_sharedkMEkMed.csv", row.names = FALSE)
rm(L58_kME_only, L58_kMed_only, L58_kMEs, L58_kMed, L58_shared)


PCGlu_UnionkME_kMeds = union(PCGlu_kMEs$Gene, PCGlu_kMed$Gene)
#write.csv(PCGlu_UnionkME_kMeds, file = "PCGlu/PCGlu_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
PCGlu_kME_only  <- setdiff(PCGlu_kMEs$Gene, PCGlu_kMed$Gene)
# Pull the genes that are unique to kMed
PCGlu_kMed_only <- setdiff(PCGlu_kMed$Gene, PCGlu_kMEs$Gene)
# Add a CR_type column
PCGlu_kME_only = data.frame(Gene = PCGlu_kME_only, CR_type = "kME")
PCGlu_kMed_only = data.frame(Gene = PCGlu_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
PCGlu_shared = intersect(PCGlu_kMEs$Gene, PCGlu_kMed$Gene)
# Save the lists for later
write.csv(PCGlu_kME_only, file = "PCGlu/PCGlu_kMEsOnly.csv", row.names = FALSE)
write.csv(PCGlu_kMed_only, file = "PCGlu/PCGlu_kMedOnly.csv", row.names = FALSE)
write.csv(PCGlu_shared, file = "PCGlu/PCGlu_sharedkMEkMed.csv", row.names = FALSE)
rm(PCGlu_kME_only, PCGlu_kMed_only, PCGlu_kMEs, PCGlu_kMed, PCGlu_shared)


R500_UnionkME_kMeds = union(R500_kMEs$Gene, R500_kMed$Gene)
#write.csv(R500_UnionkME_kMeds, file = "R500/R500_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
R500_kME_only  <- setdiff(R500_kMEs$Gene, R500_kMed$Gene)
# Pull the genes that are unique to kMed
R500_kMed_only <- setdiff(R500_kMed$Gene, R500_kMEs$Gene)
# Add a CR_type column
R500_kME_only = data.frame(Gene = R500_kME_only, CR_type = "kME")
R500_kMed_only = data.frame(Gene = R500_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
R500_shared = intersect(R500_kMEs$Gene, R500_kMed$Gene)
# Save the lists for later
write.csv(R500_kME_only, file = "R500/R500_kMEsOnly.csv", row.names = FALSE)
write.csv(R500_kMed_only, file = "R500/R500_kMedOnly.csv", row.names = FALSE)
write.csv(R500_shared, file = "R500/R500_sharedkMEkMed.csv", row.names = FALSE)
rm(R500_kME_only, R500_kMed_only, R500_kMEs, R500_kMed, R500_shared)


VT123_UnionkME_kMeds = union(VT123_kMEs$Gene, VT123_kMed$Gene)
#write.csv(VT123_UnionkME_kMeds, file = "VT123/VT123_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
VT123_kME_only  <- setdiff(VT123_kMEs$Gene, VT123_kMed$Gene)
# Pull the genes that are unique to kMed
VT123_kMed_only <- setdiff(VT123_kMed$Gene, VT123_kMEs$Gene)
# Add a CR_type column
VT123_kME_only = data.frame(Gene = VT123_kME_only, CR_type = "kME")
VT123_kMed_only = data.frame(Gene = VT123_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
VT123_shared = intersect(VT123_kMEs$Gene, VT123_kMed$Gene)
# Save the lists for later
write.csv(VT123_kME_only, file = "VT123/VT123_kMEsOnly.csv", row.names = FALSE)
write.csv(VT123_kMed_only, file = "VT123/VT123_kMedOnly.csv", row.names = FALSE)
write.csv(VT123_shared, file = "VT123/VT123_sharedkMEkMed.csv", row.names = FALSE)
rm(VT123_kME_only, VT123_kMed_only, VT123_kMEs, VT123_kMed, VT123_shared)


WO83_UnionkME_kMeds = union(WO83_kMEs$Gene, WO83_kMed$Gene)
#write.csv(WO83_UnionkME_kMeds, file = "WO83/WO83_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
WO83_kME_only  <- setdiff(WO83_kMEs$Gene, WO83_kMed$Gene)
# Pull the genes that are unique to kMed
WO83_kMed_only <- setdiff(WO83_kMed$Gene, WO83_kMEs$Gene)
# Add a CR_type column
WO83_kME_only = data.frame(Gene = WO83_kME_only, CR_type = "kME")
WO83_kMed_only = data.frame(Gene = WO83_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
WO83_shared = intersect(WO83_kMEs$Gene, WO83_kMed$Gene)
# Save the lists for later
write.csv(WO83_kME_only, file = "WO83/WO83_kMEsOnly.csv", row.names = FALSE)
write.csv(WO83_kMed_only, file = "WO83/WO83_kMedOnly.csv", row.names = FALSE)
write.csv(WO83_shared, file = "WO83/WO83_sharedkMEkMed.csv", row.names = FALSE)
rm(WO83_kME_only, WO83_kMed_only, WO83_kMEs, WO83_kMed, WO83_shared)


ZCT_UnionkME_kMeds = union(ZCT_kMEs$Gene, ZCT_kMed$Gene)
#write.csv(ZCT_UnionkME_kMeds, file = "ZCT/ZCT_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
ZCT_kME_only  <- setdiff(ZCT_kMEs$Gene, ZCT_kMed$Gene)
# Pull the genes that are unique to kMed
ZCT_kMed_only <- setdiff(ZCT_kMed$Gene, ZCT_kMEs$Gene)
# Add a CR_type column
ZCT_kME_only = data.frame(Gene = ZCT_kME_only, CR_type = "kME")
ZCT_kMed_only = data.frame(Gene = ZCT_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
ZCT_shared = intersect(ZCT_kMEs$Gene, ZCT_kMed$Gene)
# Save the lists for later
write.csv(ZCT_kME_only, file = "ZCT/ZCT_kMEsOnly.csv", row.names = FALSE)
write.csv(ZCT_kMed_only, file = "ZCT/ZCT_kMedOnly.csv", row.names = FALSE)
write.csv(ZCT_shared, file = "ZCT/ZCT_sharedkMEkMed.csv", row.names = FALSE)
rm(ZCT_kME_only, ZCT_kMed_only, ZCT_kMEs, ZCT_kMed, ZCT_shared)


PC185_UnionkME_kMeds = union(PC185_kMEs$Gene, PC185_kMed$Gene)
#write.csv(PC185_UnionkME_kMeds, file = "PC185/PC185_UnionkME_kMeds.csv", row.names = FALSE)
# Pull the genes that are unique to kME
PC185_kME_only  <- setdiff(PC185_kMEs$Gene, PC185_kMed$Gene)
# Pull the genes that are unique to kMed
PC185_kMed_only <- setdiff(PC185_kMed$Gene, PC185_kMEs$Gene)
# Add a CR_type column
PC185_kME_only = data.frame(Gene = PC185_kME_only, CR_type = "kME")
PC185_kMed_only = data.frame(Gene = PC185_kMed_only, CR_type = "kMed")
# Check that the sum of these two and the intersection add up 
PC185_shared = intersect(PC185_kMEs$Gene, PC185_kMed$Gene)
# Save the lists for later
write.csv(PC185_kME_only, file = "PC185/PC185_kMEsOnly.csv", row.names = FALSE)
write.csv(PC185_kMed_only, file = "PC185/PC185_kMedOnly.csv", row.names = FALSE)
write.csv(PC185_shared, file = "PC185/PC185_sharedkMEkMed.csv", row.names = FALSE)
rm(PC185_kME_only, PC185_kMed_only, PC185_kMEs, PC185_kMed, PC185_shared)


##### ##### ##### ##### ##### ##### 
A03_kME_only= read.csv("A03/A03_kMEsOnly.csv")
A03_kMed_only= read.csv("A03/A03_kMedOnly.csv")
A03_kME_only = A03_kME_only %>% mutate(Geno = "A03")
A03_kMed_only = A03_kMed_only %>% mutate(Geno = "A03")

HN53_kME_only= read.csv("HN53/HN53_kMEsOnly.csv")
HN53_kMed_only= read.csv("HN53/HN53_kMedOnly.csv")
HN53_kME_only = HN53_kME_only %>% mutate(Geno = "HN53")
HN53_kMed_only = HN53_kMed_only %>% mutate(Geno = "HN53")

CC168_kME_only= read.csv("CC168/CC168_kMEsOnly.csv")
CC168_kMed_only= read.csv("CC168/CC168_kMedOnly.csv")
CC168_kME_only = CC168_kME_only %>% mutate(Geno = "CC168")
CC168_kMed_only = CC168_kMed_only %>% mutate(Geno = "CC168")

R500_kME_only= read.csv("R500/R500_kMEsOnly.csv")
R500_kMed_only= read.csv("R500/R500_kMedOnly.csv")
R500_kME_only = R500_kME_only %>% mutate(Geno = "R500")
R500_kMed_only = R500_kMed_only %>% mutate(Geno = "R500")

ACC28_kME_only= read.csv("ACC28/ACC28_kMEsOnly.csv")
ACC28_kMed_only= read.csv("ACC28/ACC28_kMedOnly.csv")
ACC28_kME_only = ACC28_kME_only %>% mutate(Geno = "ACC28")
ACC28_kMed_only = ACC28_kMed_only %>% mutate(Geno = "ACC28")

ACC50_kME_only= read.csv("ACC50/ACC50_kMEsOnly.csv")
ACC50_kMed_only= read.csv("ACC50/ACC50_kMedOnly.csv")
ACC50_kME_only = ACC50_kME_only %>% mutate(Geno = "ACC50")
ACC50_kMed_only = ACC50_kMed_only %>% mutate(Geno = "ACC50")

PCGlu_kME_only= read.csv("PCGlu/PCGlu_kMEsOnly.csv")
PCGlu_kMed_only= read.csv("PCGlu/PCGlu_kMedOnly.csv")
PCGlu_kME_only = PCGlu_kME_only %>% mutate(Geno = "PCGlu")
PCGlu_kMed_only = PCGlu_kMed_only %>% mutate(Geno = "PCGlu")

ZCT_kME_only= read.csv("ZCT/ZCT_kMEsOnly.csv")
ZCT_kMed_only= read.csv("ZCT/ZCT_kMedOnly.csv")
ZCT_kME_only = ZCT_kME_only %>% mutate(Geno = "ZCT")
ZCT_kMed_only = ZCT_kMed_only %>% mutate(Geno = "ZCT")

L58_kME_only= read.csv("L58/L58_kMEsOnly.csv")
L58_kMed_only= read.csv("L58/L58_kMedOnly.csv")
L58_kME_only = L58_kME_only %>% mutate(Geno = "L58")
L58_kMed_only = L58_kMed_only %>% mutate(Geno = "L58")

PC185_kME_only= read.csv("PC185/PC185_kMEsOnly.csv")
PC185_kMed_only= read.csv("PC185/PC185_kMedOnly.csv")
PC185_kME_only = PC185_kME_only %>% mutate(Geno = "PC185")
PC185_kMed_only = PC185_kMed_only %>% mutate(Geno = "PC185")

VT123_kME_only= read.csv("VT123/VT123_kMEsOnly.csv")
VT123_kMed_only= read.csv("VT123/VT123_kMedOnly.csv")
VT123_kME_only = VT123_kME_only %>% mutate(Geno = "VT123")
VT123_kMed_only = VT123_kMed_only %>% mutate(Geno = "VT123")

FT005_kME_only= read.csv("FT005/FT005_kMEsOnly.csv")
FT005_kMed_only= read.csv("FT005/FT005_kMedOnly.csv")
FT005_kME_only = FT005_kME_only %>% mutate(Geno = "FT005")
FT005_kMed_only = FT005_kMed_only %>% mutate(Geno = "FT005")

WO83_kME_only= read.csv("WO83/WO83_kMEsOnly.csv")
WO83_kMed_only= read.csv("WO83/WO83_kMedOnly.csv")
WO83_kME_only = WO83_kME_only %>% mutate(Geno = "WO83")
WO83_kMed_only = WO83_kMed_only %>% mutate(Geno = "WO83")

##### 
A03_DiPALM = rbind(A03_kME_only, A03_kMed_only)
which(duplicated(A03_DiPALM$Gene))

HN53_DiPALM = rbind(HN53_kME_only, HN53_kMed_only)
which(duplicated(HN53_DiPALM$Gene))

CC168_DiPALM = rbind(CC168_kME_only, CC168_kMed_only)
which(duplicated(CC168_DiPALM$Gene))

R500_DiPALM = rbind(R500_kME_only, R500_kMed_only)
which(duplicated(R500_DiPALM$Gene))

ACC28_DiPALM = rbind(ACC28_kME_only, ACC28_kMed_only)
which(duplicated(ACC28_DiPALM$Gene))

ACC50_DiPALM = rbind(ACC50_kME_only, ACC50_kMed_only)
which(duplicated(ACC50_DiPALM$Gene))

PCGlu_DiPALM = rbind(PCGlu_kME_only, PCGlu_kMed_only)
which(duplicated(PCGlu_DiPALM$Gene))

ZCT_DiPALM = rbind(ZCT_kME_only, ZCT_kMed_only)
which(duplicated(ZCT_DiPALM$Gene))

VT123_DiPALM = rbind(VT123_kME_only, VT123_kMed_only)
which(duplicated(VT123_DiPALM$Gene))

FT005_DiPALM = rbind(FT005_kME_only, FT005_kMed_only)
which(duplicated(FT005_DiPALM$Gene))

L58_DiPALM = rbind(L58_kME_only, L58_kMed_only)
which(duplicated(L58_DiPALM$Gene))

PC185_DiPALM = rbind(PC185_kME_only, PC185_kMed_only)
which(duplicated(PC185_DiPALM$Gene))

WO83_DiPALM = rbind(WO83_kME_only, WO83_kMed_only)
which(duplicated(WO83_DiPALM$Gene))


## Save as an Robject for downstream analyses
Brapa_DiPALMgenes = rbind(A03_DiPALM, ACC28_DiPALM, ACC50_DiPALM, CC168_DiPALM, FT005_DiPALM, HN53_DiPALM, 
                          L58_DiPALM, PC185_DiPALM, PCGlu_DiPALM, R500_DiPALM, VT123_DiPALM, WO83_DiPALM, ZCT_DiPALM)
save(Brapa_DiPALMgenes, file = "20250919_BrapaDiPALMgenes_AllLines_withCRtype.RData")
