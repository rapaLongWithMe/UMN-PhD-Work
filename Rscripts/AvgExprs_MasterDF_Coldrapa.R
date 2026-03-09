### Creating a MasterDF to use for plotting and pulling genes of interest

# Compiling the clean DFs per genotype

setwd("~/Desktop/Git/BrapaPangenome_ColdAcclimation/DiPALM_PermResults/ExpressedGenes_FilterAcrossTreats")

# Pull in all DFs at once
files = list.files()

df = lapply(files, function(x) read_csv(x) %>% bind_rows()) 

# Clean up 
df = lapply(df, function(x) {colnames(x)[1] = "Gene"; x})

files = gsub(".csv", "", files)
files = sub("_.*", "", files)
names(df) = files

lapply(df, colnames)
lapply(df, length)


# Define the full processing function
process_expression_df <- function(df_element, name) {
  df_element %>%
    pivot_longer(cols = -Gene, names_to = "sample", values_to = "Expression") %>%
    separate(sample, into = c("Treatment", "Rep", "Timepoint"), sep = "\\.") %>%
    mutate(Genotype = sub("_.*", "", name)) %>%
    
    # Average just the replicate expression at each ZT
    group_by(Genotype, Treatment, Timepoint, Gene) %>%
    mutate(AvgExprs = mean(Expression), 
           sdExprs = sd(Expression)) %>%
    ungroup() %>%
    
    # Calculate Avg and sd over the whole TC by gene
    group_by(Genotype, Treatment, Gene) %>%
    mutate(TCavg = mean(Expression), 
           TCsd = sd(Expression)) %>%
    ungroup() %>%
    
    # Now calculate the zScore to standardize gene-level expression
    mutate(zScore = (AvgExprs - TCavg) / TCsd) %>%
    as.data.frame()  # keep as base data.frame
}

## Apply the function
df_processed <- Map(process_expression_df, df, names(df))

## Save for now 
save(df_processed, file = "20250417_df_processed.RData")
rm(df)

## Make a dataframe by binding rows
exprsDF = do.call(rbind, df_processed)
  # This should be faster than the above
combined_df <- rbindlist(df_processed, idcol = "Genotype")
exprsDF$Gene = gsub("\\.v[0-9].[0-9]", "", exprsDF$Gene)

exprsDF$ZT = exprsDF$Timepoint %>% 
  str_replace_all(c("1" = "17", "2" = "21", "3" = "1", 
                    "5" = "9", "4" = "5", "6" = "13"))

gc()

## Check some things
levels(factor(exprsDF$Genotype))
levels(factor(exprsDF$ZT))
levels(factor(exprsDF$Treatment))
levels(factor(exprsDF$Rep))

  # CCA1
"BraR500.05G351900"
"BraR500.05G351800"
"BrVT123.05G347600"
"BrWO_83.05G350200"
"BrA03.05G351900"
"BrL58.05G348400"
"BrPCGlu.05G349500"

  # TOC1
"BrPCGlu.03G155200"
"BrPCGlu.09G066000"
"BraR500.03G173600"
"BraR500.09G068800"
"BrVT123.03G155300"
"BrVT123.09G066900"
"BrWO_83.03G153200"
"BrA03.03G155400"
"BrA03.09G065800"
"BrL58.03G154100"
"BrL58.09G065900"

exprsDF %>% 
  filter(Gene == "BrL58.09G065900") %>%
  ggplot(aes(x = factor(ZT, levels = c("17", "21", "1", "5", "9", "13")), 
             y = zScore, group = interaction(Treatment, Rep))) +
  geom_line(aes(color = Treatment)) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  facet_wrap(~Genotype) +
  xlab("Timepoint") +
  ylab("zScore") +
  #ggtitle("CCA1")
  ggtitle("TOC1")

## Looks good so far...
load("~/Desktop/Git/BrapaPangenome_ColdAcclimation/Brapa_MasterDF.RData")
Brapa_MasterDF = exprsDF
rm(exprsDF, df_processed)
save(Brapa_MasterDF, file = "~/Desktop/Git/BrapaPangenome_ColdAcclimation/Brapa_MasterDF.RData")
  # Also nice to have just a list of the expressed genes as a simple vector 
exprsGenes_AllGenos = Brapa_MasterDF %>% pull(Gene) %>% unique()
write.csv(exprsGenes_AllGenos, file = "~/Desktop/Git/BrapaPangenome_ColdAcclimation/AllGenos_ExprsGenes_AsVector.csv", row.names = FALSE)

## Now pull out just the average expression to make individual (and a joined) table for each genotype
Brapa_MasterLst = split(Brapa_MasterDF, Brapa_MasterDF$Genotype)
Brapa_AvgLst = lapply(Brapa_MasterLst, function(x) x = x[, -c(3,5)] %>% unique())
  # To save space, write them all out as a compiled list that I can bring in and pull out each geno
save(Brapa_AvgLst, file = "~/Desktop/Git/BrapaPangenome_ColdAcclimation/Brapa_CompiledAvgList.RData")

