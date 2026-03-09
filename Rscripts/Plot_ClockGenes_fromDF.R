## Plotting gene expression of genes of interest

# Set up the environment
setwd("~/Desktop/Git/Brapa_Pangenome_ColdTC/ExpressedGenes/CleanedAndChecked")
require(tidyverse)
require(ggthemes)
require(readxl)

# Load the basic DF previously made
##### CHANGE THE NAME OF THE GENOTYPE #####
Geno = "R500"
##### ##### ##### ##### ##### ##### ##### 

# Load
file = paste0(Geno, "_DF.RData")
load(file)

## Set up the expression data; here I am using a basic DF object that Benjamin Sebastian created in his DiPALM work
rownames(DF)

## Pivot the DF to create plottable columns
  # Rownames get lost when pivoting so bring those in first
exprsDF = rownames_to_column(DF, var = "Gene")
exprsDF = exprsDF %>% pivot_longer(!(Gene), names_to = "ZTs", values_to = "Expression")

# Separate ZT column into multiple
exprsDF = exprsDF %>% separate(ZTs, c("Geno", "Rep", "Treatment", "ZT"))
gc()
# Make and IDs column to pull on later
exprsDF = exprsDF %>% unite("IDs", Geno, Gene, remove = FALSE)

# Will have to reorder ZTs
levels(factor(exprsDF$ZT))
levels(factor(exprsDF$Treatment))
levels(factor(exprsDF$Rep))

## Add average, sd, and zscore (standardized) expression columns for plotting
exprsDF = exprsDF %>% 
  # Average just the replicate expression at each ZT
  group_by(Geno, Treatment, ZT, Gene) %>%
  mutate(AvgExprs = mean(Expression), 
         sdExprs = sd(Expression)) %>%
  ungroup() %>% 
  # Calculate Avg and sd over the whole TC by gene
  group_by(Geno, Treatment, Gene) %>% 
  mutate(TCavg = mean(Expression), 
         TCsd = sd(Expression)) %>%
  ungroup() %>%
  # Now calculate the zScore to standardize gene level (not rep) expression
  mutate(zScore = (AvgExprs - TCavg)/TCsd)

## Note the ZTs are out of order; fix
exprsDF = exprsDF %>% mutate(ZT = factor(ZT, levels = c("1", "5", "9", "13", "17", "21")))
factor(levels(exprsDF$ZT))

# Will have to reorder ZTs
exprsDF = exprsDF %>% mutate(ZT = factor(ZT, levels = c("17", "21", "1", "5", "9", "13")))
levels(exprsDF$ZT)
# Clean up the IDs column
exprsDF <- exprsDF %>%
  mutate(
    # Extract characters after 'Br' and before the first period to get the reference line the genes were mapped to
    Reference = str_extract(IDs, "(?<=Br)[^\\.]+"), 
    # Extract up to the second period to get the gene annotation without the splice variant
    Gene = str_extract(IDs, "[^\\.]+\\.[^\\.]+"))
# Remove the Geno ID in the Gene annotation
exprsDF$Gene = gsub("R500_", "", exprsDF$Gene)

cca1 = "BraR500.05G019300"

exprsDF %>% filter(Gene == cca1) %>%
  ggplot(aes(x = ZT, y = AvgExprs)) +
  geom_line(aes(group = interaction(Rep, Treatment), 
                color = Treatment)) +
  geom_ribbon(aes(ymin = AvgExprs - sdExprs, ymax = AvgExprs + sdExprs, 
                  group = interaction(Rep, Treatment)), alpha = .1) +
  scale_color_manual(values = c("blue", "red")) +
  ylab("Average Log2 TPM")




pdf(file.path("CNSgoi", file = paste(Geno, "_CNSgoi_clockParalogs.pdf")), width = 6, height = 4)
for(i in 1:length(clock_list)){
  print(clock_list[[i]] %>%
          ggplot(aes(x = ZT, y = AvgExprs)) +
          geom_line(aes(group = interaction(Rep, Treatment), 
                        color = Treatment)) +
          geom_ribbon(aes(ymin = AvgExprs - sdExprs, ymax = AvgExprs + sdExprs, 
                          group = interaction(Rep, Treatment)), alpha = .1) +
          scale_color_manual(values = c("blue", "red")) +
          ylab("Average Log2 TPM") +
          facet_wrap(~ClockParalog) +
          ggtitle(paste0(Geno, ":", names(clock_list[i]))) +
          theme_tufte() +
          theme(axis.line = element_line(colour = "grey50"), 
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                strip.text = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
                axis.title.y = element_text(size = 16), 
                axis.title.x = element_text(size = 16), 
                title = element_text(size = 16)))
}
dev.off()