### Creating average expression and master data sheets for Bnapus drought


setwd("~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses")
inPath = "./ExpressedGenes_Final/Avg_DFs"

require(tidyverse)


## Bring in the average expression by geno that has already been checked
all.files = list.files(file.path(inPath), pattern = ".RData")
all.files

exprs <- lapply(all.files, function(x) get(load(file.path(inPath, file = x))))
names(exprs) <- all.files 


# AvgExprs List
AvgExprs <- list(Ab = exprs$Ab_Avg_DF.RData, Al = exprs$Al_Avg_DF.RData, Av = exprs$Av_Avg_DF.RData,
              Br = exprs$Br_Avg_DF.RData, Ca = exprs$Ca_Avg_DF.RData, Da = exprs$Da_Avg_DF.RData, 
              DH12 = exprs$DH12_Avg_DF.RData, DH20 = exprs$DH20_Avg_DF.RData, Gr = exprs$Gr_Avg_DF.RData, 
              Mu = exprs$Mu_Avg_DF.RData, Ne = exprs$Ne_Avg_DF.RData, Se = exprs$Se_Avg_DF.RData, 
              St = exprs$St_Avg_DF.RData, Yu = exprs$Yu_Avg_DF.RData)

# Rownames has the geno info already, so rbind into one DF
AvgExprs = lapply(AvgExprs, function(x) rownames_to_column(x, var = "IDs"))
AvgExprs = do.call(rbind, AvgExprs)
AvgExprs = AvgExprs %>% separate(IDs, c("Geno", "Gene"), remove = FALSE)

which(is.na(AvgExprs))
levels(factor(AvgExprs$Geno))

# Save
getwd()
save(AvgExprs, file = "AllGenos_AvgExprs_DF.RData")



### Make one large expression DF that is not averaged
inPath = "./ExpressedGenes/DFs"

## Bring in the expression DFs by geno 
all.files = list.files(file.path(inPath), pattern = ".RData")
all.files

exprs <- lapply(all.files, function(x) get(load(file.path(inPath, file = x))))
names(exprs) <- all.files 

## Make a list
exprs <- list(Ab = exprs$Ab_DF.RData, Al = exprs$Al_DF.RData, Av = exprs$Av_DF.RData,
                 Br = exprs$Br_DF.RData, Ca = exprs$Ca_DF.RData, Da = exprs$Da_DF.RData, 
                 DH12 = exprs$DH12_DF.RData, DH20 = exprs$DH20_DF.RData, Gr = exprs$Gr_DF.RData, 
                 Mu = exprs$Mu_DF.RData, Ne = exprs$Ne_DF.RData, Se = exprs$Se_DF.RData, 
                 St = exprs$St_DF.RData, Yu = exprs$Yu_DF.RData)

## Remove the Geno ID from colnames to bind
exprs = lapply(exprs, function(x) {colnames(x) = gsub("^[^_]*\\_", "", colnames(x));x})


# Rownames has the Geno Treatment Gene info
exprsDF = lapply(exprs, function(x) rownames_to_column(x, var = "IDs"))
  # Make a DF
exprsDF = do.call(rbind, exprsDF)

# Pivot to make plottable columns
exprsDF = exprsDF %>% pivot_longer(!(IDs), 
                               names_to = "ZTs", values_to = "Expression")
  # separate ZT column into multiple
exprsDF = exprsDF %>% separate(ZTs, c("ZTs", "Treatment", "Rep"))
  # split IDs to pull out geno and gene, remove this treatment cause it no longer matches after pivoting
exprsDF = exprsDF %>% separate(IDs, c("Geno", "Remove", "Gene"))
  # remove the extraneous column now
exprsDF = exprsDF[, -2]
 
# check that the number of obvs makes sense
obs = unlist(lapply(exprs, function(x) nrow(x)))
sum(obs*8*3)
  # Make sure levels don't get lost in the process
levels(factor(exprsDF$Treatment))


# Save for future use
save(exprsDF, file = "checked_exprsDF.RData")


## Add average and sd columns to plot
test = exprsDF %>% 
  group_by(Geno, Treatment, ZTs, Gene) %>%
  rowwise() %>%
  mutate(AvgExprs = mean(Expression), 
         sdExprs = sd(Expression)) %>%
  ungroup() %>% 
  group_by(Geno, Treatment,Gene) %>%
  rowwise() %>%
  mutate(NormExprs = Expression/AvgExprs, 
         sdNorm = sd(NormExprs)) %>%
  ungroup()


# Check by ploting CCA1 average expression for each genotype
r_cca1.1 = "BraA05g01930R"
  
test %>% 
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

# Try it with Normalized expression
test %>% 
  filter(Gene == "BraA05g01930R") %>%
  ggplot(aes(x = factor(ZTs, levels = c("1", "7", "13", "19")), y = NormExprs)) +
  geom_line(aes(color = Treatment, group = interaction(Rep, Treatment))) +
  #geom_ribbon(aes(ymin = NormExprs - sdNorm, ymax = NormExprs + sdNorm, group = interaction(Rep, Treatment), fill = Treatment), alpha = .1) +
  scale_color_manual(values = c("darkgoldenrod1", "chartreuse4")) +
  scale_fill_manual(values = c("darkgoldenrod1", "chartreuse4")) +
  theme_bw() +
  facet_wrap(~Geno) +
  xlab("Zts (hours after lights on)") +
  ylab("Log2 FPKM") +
  ggtitle("CCA1")


blah = test %>% filter(Gene == "BraA05g01930R")


  
