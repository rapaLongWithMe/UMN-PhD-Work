### Creating a MasterDF to use for plotting and pulling genes of interest

# Compiling the clean DFs per genotype

setwd("~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses")

# Pull in all DFs at once
files = list.files(file.path("ExpressedGenes/DFs/"))

exprsDF <- lapply(files, function(x) get(load(file.path("ExpressedGenes/DFs/", file = x))))
names(exprsDF) <- files 

# Note they all say drought in the rownames so get rid of that
lapply(exprsDF, rownames)
exprsDF = lapply(exprsDF, function(x) rownames_to_column(x, var = "IDs"))
exprsDF = lapply(exprsDF, function(x) {x$IDs = gsub("Drought_", "", x$IDs); x})


## Remove the Geno ID from colnames to be able to bind columns (names must match in all elements)
exprsDF = lapply(exprsDF, function(x) {colnames(x) = gsub("^[^_]*\\_", "", colnames(x));x})

# Make a dataframe by binding rows
exprsDF = do.call(rbind, exprsDF)

# Pivot to make plottable columns
exprsDF = exprsDF %>% pivot_longer(!(IDs), 
                                   names_to = "ZTs", values_to = "Expression")

# Separate ZT column into multiple. Takes a few minutes.
exprsDF = exprsDF %>% separate(ZTs, c("ZTs", "Treatment", "Rep"))
gc()
  # Will have to reorder ZTs
levels(factor(exprsDF$ZTs))
levels(factor(exprsDF$Treatment))
levels(factor(exprsDF$Rep))

# split IDs to pull out geno and gene; keep IDs for an index
exprsDF = exprsDF %>% separate(IDs, c("Geno", "Gene"), remove = FALSE)
gc()
levels(factor(exprsDF$ZTs))
levels(factor(exprsDF$Treatment))
levels(factor(exprsDF$Rep))
levels(factor(exprsDF$Geno))
tail(exprsDF)

# Keep only Bra* or Bo* genes (ex: ENSRA genes); should be good but might as well
exprsDF = exprsDF[str_detect(exprsDF$Gene, "^Bra|Bo"),]
which(is.na(exprsDF))


## Add average, sd, and zscore (standardized) expression columns for plotting
exprsDF = exprsDF %>% 
    # Average just the replicate expression at each ZT
  group_by(Geno, Treatment, ZTs, Gene) %>%
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


# Check by ploting CCA1 average expression for each genotype
r_cca1.1 = "BraA05g01930R"

exprsDF %>% 
  filter(Gene == "BraA05g01930R") %>%
  ggplot(aes(x = factor(ZTs, levels = c("1", "7", "13", "19")), y = AvgExprs)) +
  geom_line(aes(color = Treatment, group = interaction(Rep, Treatment))) +
  geom_ribbon(aes(ymin = AvgExprs - sdExprs, ymax = AvgExprs + sdExprs, group = interaction(Rep, Treatment), fill = Treatment), alpha = .1) +
  scale_color_manual(values = c("darkgoldenrod1", "chartreuse4")) +
  scale_fill_manual(values = c("darkgoldenrod1", "chartreuse4")) +
  theme_bw() +
  facet_wrap(~Geno) +
  xlab("Zts (hours after lights on)") +
  ylab("Log2 FPKM") +
  ggtitle("CCA1")


# Try it with Standardized expression
exprsDF %>% 
  filter(Gene == "BraA05g01930R") %>%
  ggplot(aes(x = factor(ZTs, levels = c("1", "7", "13", "19")), y = zScore, group = interaction(Treatment, Rep))) +
  geom_line(aes(color = Treatment)) +
  scale_color_manual(values = c("darkgoldenrod1", "chartreuse4")) +
  scale_fill_manual(values = c("darkgoldenrod1", "chartreuse4")) +
  theme_bw() +
  facet_wrap(~Geno) +
  xlab("Zts (hours after lights on)") +
  ylab("zScore") +
  ggtitle("CCA1")



# Bring in the responseScores object with module information to bind together with expression data
load("~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses/ResponseScore/20230326_ResponseScoresAndMods_Only.RData")

tail(responseScores)
levels(factor(responseScores$Module))
which(str_detect(responseScores$IDs, "^Bra|^Bo"))

# check that the numbers match; might have more responseScore rows since I have removed some DR genes based on singleTP analyses
nrow(responseScores) * 2 * 4
nrow(exprsDF) - (nrow(responseScores) * 2 * 3)

test = unique(left_join(exprsDF, responseScores[, c(1, 6)]))


# Save as a RData object to keep working on
save(test, file = "Draft_MasterDF.RData")

