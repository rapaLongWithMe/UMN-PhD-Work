
### Some troubleshooting/package fixes
# install old packages with the repmis package (example is the dependency 'Matrix' for the WGCNA package)
repmis::InstallOldPackages(pkg='Matrix', version='1.6-5')


### Collective reference of Code snippets that I use often 

# remove rows with non Bra/Bo genes (ENSRNA genes)
kME_Sig.01 = kME_Sig.01[which(!grepl("Bra|Bo", kME_Sig.01$Gene) == FALSE), ]
## Or
# Remove any annotation (ie ENSRNA genes) that are not Bra or Bo
TCs_List = lapply(TCs_List, function(x) x[which(str_detect(rownames(x), "Bra|Bo")== TRUE),])


### Chop up a gene name (or character string) based on delimiters
exprsDF <- exprsDF %>%
  mutate(
    # Extract characters after 'Br' and before the first period to get the reference line the genes were mapped to
    Reference = str_extract(IDs, "(?<=Br)[^\\.]+"), 
    # Extract up to the second period to get the gene annotation without the splice variant
    Gene = str_extract(IDs, "[^\\.]+\\.[^\\.]+"))

## Similar result but outside of tidyverse and don't have to make new columns 
exprsDF$IDs = gsub("\\.v[0-9].[0-9]", "", exprsDF$IDs)


# Create a theme for all plots 
my_theme = theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = .8, hjust = .8), 
                              text = element_text(family = "Times"), 
                              axis.title = element_text(size = 14), legend.title = element_blank(),
                              legend.position = "bottom", plot.title = element_text(size = 18))



## Plotting gene expression in a nested list
  ## Great for when you want to plot many genes across many genos at once
load("~/Desktop/Bnapus_drought/Bnapus_TranscriptomeAnalyses/DistributionPlots/Figures/nestedLst.RData")

## Here's how I created nestedLst
# Make a list to loop through and plot
load(file = "~/Desktop/Bnapus_drought/Bnapus_TranscriptomeAnalyses/DistributionPlots/Figures/20240208_OverlapGenes_dat.RData")
datLst = split(dat, dat$InvolvedIn)
nestedLst = lapply(datLst, function(x) split(x, x$Geno))

##### THIS CURRENTLY PRINTS NA AFTER THE FIRST PLOT IN THE TOP LAYER ELEMENT....NEED TO FIX ######
pdf("DistributionPlots/Figures/20240208_WL_ParalogOverlap.pdf")
for(i in 1:length(nestedLst)){
  for(j in 1:length(nestedLst[[i]])){
    plot(nestedLst[[i]][[j]] %>%
           ggplot(aes(x = factor(ZTs, levels = c("1", "7", "13", "19")), y = AvgExprs, 
                      color = Treatment, fill = Treatment, group = interaction(IDs, Treatment))) +
           geom_line(aes(linetype = IDs)) +
           geom_ribbon(aes(ymin = AvgExprs - sdExprs, ymax = AvgExprs + sdExprs), alpha = .1
                       # use this if want to turn off the outlines of the ribbon
                       , linetype = 0) +
           scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
           scale_fill_manual(values = c("darkgoldenrod2", "chartreuse4")) +
           theme_bw() +
           xlab("ZT (hours after lights on)") +
           ylab("Average Expression") +
           ggtitle(paste(names(nestedLst[i][j]))) +
           # Remove legend entirely
           # theme(legend.position = "bottom")
           theme(legend.position = "bottom"))
  }
}
dev.off()


#### More plotting tricks
# In mixed grouped bar plots, maintain the size of the bar width across singles and groups

  #geom_col(aes(color = Geno, fill = Geno, group = interaction(Description, Geno)), 
   #        position = position_dodge2(preserve = c("single"
#, "total")
#),width = .6) +
 

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
    ## Take a character string of gene annotations that currently aren't individual strings (one string per many annotations) and separate them into unique rows
## Load the updated excel sheet with the simplified Descriptors of interest
Freq = read_xlsx("AllModsOverEnrich_BroadDescript.xlsx", sheet = "CarbonNitrogenSulfurPhoto")

# Prepare the data
# Creates a set of functions to apply to this list to do the clean up 
masterFunctions = function(x) {
  #clean up and remove the pipe
  x = data.frame(x %>%
                   mutate(GenesWithAnn = str_replace_all(GenesWithAnn, '\\|', "")))
}


# Split the DF on Broad descriptor and apply the function to all elements
master = split(Freq, Freq$BroadDescriptor)
master = lapply(master, masterFunctions)
DF = rbindlist(master)

  ##### MAY HAVE BEEN SUPERCEEDED WITH SEPARATE_LONGER_DELIM #####
DF = DF %>%
  separate_rows(GenesWithAnn, sep = "  ") 

# Calculate the number of enriched genes by Term in a given module, for each genotype
DF = DF %>% group_by(Geno, Module, GO_Term) %>% mutate(NumGenes = n()) %>% ungroup()

# Clean the data. Want a figure that shows number of genes enriched in a given process in a certain module
DF = DF[, c(4, 5, 6, 10)]
# Remove the Amide metabolic processes
DF = DF %>% filter(BroadDescriptor != "Amide metabolic process")
# Change the module labels from colors to numbers 
levels(factor(DF$Module))

DF$Module = DF$Module %>% str_replace_all(c("black_dudd" = "M5A", "black_uuud" = "M5B", "blue" = "M2", "brown" = "M1", "green" = "M4", 
                                            "greenyellow" = "M10", "magenta" = "M11", "pink" = "M7", "purple" = "M6",
                                            "red_dddu" = "M9A","red_duuu" = "M9B", "turquoise" = "M8", "yellow" = "M3"))

# NOTE that the str_replace function labels 'greenyellow' as 'green' plus 'yellow' which makes M10 = 43
levels(factor(DF$Module))
DF$Module = gsub("M4M3", "M10", DF$Module)
# Reorder
DF = DF %>% mutate(Module = factor(Module, levels = c("M1", "M2", "M3", "M4", "M5A", "M5B", "M6", 
                                                      "M7", "M8", "M9A", "M9B", "M10", "M11")))




##### ##### ##### ##### ##### ##### ##### ##### 
# Add a subgenome identifier to color by 
WL_geneLst = lapply(WL_geneLst, function(x) x %>% mutate(Subgenome = case_when(str_detect(Gene, "Bra") ~ "rapa", 
                                                                               TRUE ~ "ole")))
##### ##### ##### ##### ##### ##### 

##### ##### ##### ##### ##### ##### 
# Adjusting themes to fit within the width of the plot

#axis.text.x changes the font of the xticks
theme(axis.text.x = element_text(size = 10),
      legend.position = "bottom",legend.title = element_blank(),  # Removes legend title for compactness
      legend.text = element_text(size = 10),  # Adjust legend text size
      legend.key.width = unit(1, "cm"),  # Adjust key width if necessary
      legend.box = "horizontal") +  # Arrange legend items horizontally
  guides(color = guide_legend(nrow = 2))  # Arrange legend items in two rows

# Adjusting aesthetics of the color palette 
# Create a named palette to maintain coloring consistency across (ex) genes
paralogPalette <- c("Bo2g161590" = "#456789","Bo7g098590" = "#0A9396","Bo9g014610" = "#94D2BD",
                    "BraA02g44060R" = "#FEAC72", "BraA03g42560R" = "#E47E62","BraA09g07610R" = "#D9A0A0")
# This named list goes in the scale_fill_manaul(values = paralogPalette)
##### ##### ##### ##### ##### ##### 

# Bring in the table from KT with the genes/genos of interest
overlap = read_xlsx("DistributionPlots/KT_Overlap-plotting-candidates.xlsx")
# Simplify to make it easier to match with expression
SmOverlap = overlap[, c(5,6,7,10,11)]
colnames(SmOverlap)[2] = "Geno"
colnames(SmOverlap)[3] = "Gene"
# Pull out the Gene function (name) with the corresponding ATG and matching paralogs
ogs = overlap[, c(1, 2, 5)]
# Fill in the NAs; the function stands for 'last observation carried forward'
ogs$...1 = na.locf(ogs$...1)


##### ##### ##### ##### ##### #####

## Make a list a dataframe again
test = do.call(rbind.data.frame, All_kMEs_ExprsShift)


## Pull out just the numbers for ZT
AvgExprs$Drought_Max = gsub(".*?(\\d.*)", "\\1", AvgExprs$Drought_Max)
AvgExprs$Watered_Max = gsub(".*?(\\d.*)", "\\1", AvgExprs$Watered_Max)

##### Different dataset but still list manipulation
## gsub a named vector within a large list
test = lapply(kMEs, function(x) {names(x) = gsub("gene:", "", x); x})


##### REGEX and gsub #####
  # Get rid of everything after the last underscore

gsub('_[^_]*$', '', colnms)
