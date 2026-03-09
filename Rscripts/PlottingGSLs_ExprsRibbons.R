## Want to plot standardized expression with ribbons
  # For r500, hn53 and Arabidopsis 
  # select gsl genes (already have written out)


## GLUCOSINOLATE AND SULFUR GENES ## 
myb28_A02 = "BraA02g44060R"
myb28_A03 = "BraA03g42560R"
myb28_A09 = "BraA09g07610R"

myb29_A03 = "BraA03g03990R"
myb29_A10 = "BraA10g28030R"

cyp83B1 = "BraA08g07180R"

cyp83A1 = "BraA04g09110R"

mam1_A02 = "BraA02g27580R"	
mam1_A03 = "BraA03g42970R"

mam3 = "BraA02g27590R"

paps_A08 = "BraA08g28310R"
paps_A06 = "BraA06g38620R"
paps_A09 = "BraA09g55040R"

sur1_A07 = "BraA07g01210R"
sur1_A07_Bo = "Bo7g003330"
sur1_A09 = "BraA09g12830R"

sdi1_2g = "Bo2g150270"
sdi1_9g = "Bo9g009620"

cyp81F1_A10 = "BraA10g14510R"
cyp81F1_A02 = "BraA02g13430R"
cyp81F1_A03 = "BraA03g12530R"

gtr1_A01 = "BraA01g15510R"
gtr1_A06 = "BraA06g18650R"

gtr2_A06 = "BraA06g24760R"
gtr2_A09 = "BraA09g08070R"

pen3 = "BraA07g23920R"

esp_A06_40890 = "BraA06g40890R"
esp_A05 = "BraA05g16440R"
esp_A06_1400 = "BraA06g01400R"

tgg1_A01 = "BraA01g24160R"
tgg1_A02 = "BraA02g41930R"
tgg1_A03 = "BraA03g27140R"
tgg1_A07 = "BraA07g28860R"
tgg1_A08 = "BraA08g29400R"
tgg1_A09 = "BraA09g10900R"


#put into a list to loop through while plotting
gsl = list(myb28_A02,myb28_A03,myb28_A09, myb29_A03,myb29_A10,
           cyp83B1,cyp83A1,mam1_A02,mam1_A03,mam3,paps_A08,
           paps_A06,paps_A09,sur1_A07,sur1_A09,cyp81F1_A10,
           cyp81F1_A02,cyp81F1_A03,gtr1_A01,gtr1_A06,gtr2_A06,
           gtr2_A09,pen3,esp_A06_40890,esp_A05,esp_A06_1400,
           tgg1_A01,tgg1_A02,tgg1_A03,tgg1_A07,tgg1_A08,tgg1_A09)

names(gsl) = c("MYB28_A02","MYB28_A03","MYB28_A09", "MYB29_A03","MYB29_A10",
               "CYP83B1","CYP83A1","MAM1_A02","MAM1_A03","MAM3","PAPS_A08",
               "PAPS_A06","PAPS_A09","SUR1_A07","SUR1_A09","CYP81F1_A10",
               "CYP81F1_A02","CYP81F1_A03","GTR1_A01","GTR1_A06","GTR2_A06",
               "GTR2_A09","PEN3","ESP_A06_40890","ESP_A05","ESP_A06_1400",
               "TGG1_A01", "TGG1_A02","TGG1_A03","TGG1_A07","TGG1_A08","TGG1_A09")


# Just need the expression dfs
  # R500
load("/Users/ricon001/Desktop/Git/RNAseq/Brapa/R500/Pre_2022/R_output/r500_samps.RData")
  # Bnapus
load("/Users/ricon001/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses/ResponseScore/PrevData/20221118_MasterDF_plusPvalues.RData")



# Bring in rownames
r500_exprs = lapply(r500_samps, function(x) {x = rownames_to_column(data.frame(x), var = "Gene");x })

r500_exprs = do.call(cbind, r500_exprs)
  # Keep one gene column since the rest are duplicated
r500_exprs =r500_exprs %>% select(!(c(Cold_R2.Gene, Cold_R3.Gene,Cold_R4.Gene, 
                             Warm_R1.Gene, Warm_R2.Gene, Warm_R3.Gene, Warm_R4.Gene)))
  # Remove the extra info before the period
colnames(r500_exprs) = gsub("Cold_R[0-9].", "", colnames(r500_exprs))
colnames(r500_exprs) = gsub("Warm_R[0-9].", "", colnames(r500_exprs))

# Pivot to have columns to plot on
r500_exprs = r500_exprs %>% pivot_longer(cols = !Gene, names_to = "ColID", values_to = "Exprs")
  # Separate the ColID
r500_exprs = r500_exprs %>% separate(ColID, c("Rep", "Treatment", "Zt"))
  # Rename tps
r500_exprs$Zt = str_replace_all(r500_exprs$Zt, c("1" = "ZT_17", "2" = "ZT_21", "3" = "ZT_1", 
                                                 "4" = "ZT_5", "5" = "ZT_9", "6" = "ZT_13"))
levels(factor(r500_exprs$Zt))

# Remove the extra tp
r500_exprs = r500_exprs %>% filter(Zt != "7")

#the 5 gets replaced so fix that 
r500_exprs$Zt = str_replace_all(r500_exprs$Zt, c("ZT_ZT_9" = "ZT_5"))

levels(factor(r500_exprs$Zt))
r500_exprs$Zt <- factor(r500_exprs$Zt, levels=c("ZT_17", "ZT_21", "ZT_1", 
                                                "ZT_5", "ZT_9", "ZT_13"))


## Calculate average and sd to plot
r500_exprs = r500_exprs %>% group_by(Treatment, Zt, Gene) %>%
  mutate(AvgExprs = mean(Exprs), 
         sd = sd(Exprs)) %>% ungroup()


# Check
r500_exprs %>% 
  filter(Gene == gsl[[1]] & Treatment == "Warm") %>%
  ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
  geom_line() +
  geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .1)

# save the df to plot in the future
save(r500_exprs, file = "r500_exprs_Ribbons.RData")


setwd("~/Desktop/Git/RNAseq")
pdf("GSL_Control_r500Alone.pdf",width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(gsl)) {
  print(r500_exprs %>% 
          filter(Treatment == "Warm" & Gene == gsl[[i]]) %>%
          ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
          geom_line(color = "steelblue", linewidth = 1.2) +
          geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .1, fill = "steelblue4") +
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average expression (control)") +
          ggtitle(paste(names(gsl)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()


### REPEAT WITH NEXT GENOTYPE
# Just need the expression dfs
load("/Users/ricon001/Desktop/Git/RNAseq/Brapa/Hn53/R_output/hn53_samps.RData")

# Bring in rownames
hn53_exprs = lapply(hn53_samps, function(x) {x = rownames_to_column(data.frame(x), var = "Gene");x })

hn53_exprs = do.call(cbind, hn53_exprs)
# Keep one gene column since the rest are duplicated
hn53_exprs =hn53_exprs %>% select(!(c(Cold_R2.Gene, Cold_R3.Gene,Cold_R4.Gene, 
                                      Warm_R1.Gene, Warm_R2.Gene, Warm_R3.Gene, Warm_R4.Gene)))
# Remove the extra info before the period
colnames(hn53_exprs) = gsub("Cold_R[0-9].", "", colnames(hn53_exprs))
colnames(hn53_exprs) = gsub("Warm_R[0-9].", "", colnames(hn53_exprs))

# Pivot to have columns to plot on
hn53_exprs = hn53_exprs %>% pivot_longer(cols = !Gene, names_to = "ColID", values_to = "Exprs")
  # Remove the ZT characters
hn53_exprs$ColID = gsub("_ZT", "", hn53_exprs$ColID)
  # Separate the ColID
hn53_exprs = hn53_exprs %>% separate(ColID, c("Treatment", "Zt", "Rep"))
  # Rename tps
hn53_exprs$Zt = str_replace_all(hn53_exprs$Zt, c("17" = "ZT_17", "21" = "ZT_21", "25" = "ZT_1", 
                                                 "29" = "ZT_5", "33" = "ZT_9", "37" = "ZT_13"))
levels(factor(hn53_exprs$Zt))

# Remove the extra tp
hn53_exprs = hn53_exprs %>% filter(Zt != "41")

hn53_exprs$Zt <- factor(hn53_exprs$Zt, levels=c("ZT_17", "ZT_21", "ZT_1", 
                                                "ZT_5", "ZT_9", "ZT_13"))

## Calculate average and sd to plot
hn53_exprs = hn53_exprs %>% group_by(Treatment, Zt, Gene) %>%
  mutate(AvgExprs = mean(Exprs), 
         sd = sd(Exprs)) %>% ungroup()


# Check
hn53_exprs %>% 
  filter(Gene == gsl[[1]] & Treatment == "Warm") %>%
  ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
  geom_line() +
  geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .1)


# save the df to plot in the future
save(hn53_exprs, file = "hn53_exprs_Ribbons.RData")



setwd("~/Desktop/Git/RNAseq")
pdf("GSL_Control_hn53Alone.pdf",width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(gsl)) {
  print(hn53_exprs %>% 
          filter(Treatment == "Warm" & Gene == gsl[[i]]) %>%
          ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
          geom_line(color = "chocolate", linewidth = 1.2) +
          geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .3, fill = "chocolate4") +
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average expression (control)") +
          ggtitle(paste(names(gsl)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()


### AND AGAIN WITH THE TWO ARABIDOPSIS LINES 
# Just need the expression dfs
load("/Users/ricon001/Desktop/Git/RNAseq/Arabidopsis/R_output/C24/C24_TCs_List.RData")

# Bring in rownames
c24_exprs = lapply(TCs_List, function(x) {x = rownames_to_column(data.frame(x), var = "Gene");x })

c24_exprs = do.call(cbind, c24_exprs)
# Keep one gene column since the rest are duplicated
c24_exprs =c24_exprs %>% select(!(c(Freeze_2.Gene, Freeze_3.Gene, Freeze_4.Gene, 
                                    Control_1.Gene, Control_2.Gene, Control_3.Gene, Control_4.Gene)))
# Remove the extra info before the period
colnames(c24_exprs) = gsub("Freeze_[0-9].", "", colnames(c24_exprs))
colnames(c24_exprs) = gsub("Control_[0-9].", "", colnames(c24_exprs))

# Pivot to have columns to plot on
c24_exprs = c24_exprs %>% pivot_longer(cols = !Gene, names_to = "ColID", values_to = "Exprs")
# Separate the ColID
c24_exprs = c24_exprs %>% separate(ColID, c("Geno", "Rep", "Treatment", "Zt"))
# Rename tps
c24_exprs$Zt = str_replace_all(c24_exprs$Zt, c("ZT17" = "ZT_17", "ZT21" = "ZT_21", "ZT1" = "ZT_1", 
                                                 "ZT5" = "ZT_5", "ZT9" = "ZT_9", "ZT13" = "ZT_13"))
levels(factor(c24_exprs$Zt))


levels(factor(c24_exprs$Zt))
c24_exprs$Zt <- factor(c24_exprs$Zt, levels=c("ZT_17", "ZT_21", "ZT_1", 
                                                "ZT_5", "ZT_9", "ZT_13"))

## Calculate average and sd to plot
c24_exprs = c24_exprs %>% group_by(Treatment, Zt, Gene) %>%
  mutate(AvgExprs = mean(Exprs), 
         sd = sd(Exprs)) %>% ungroup()

# Remove the Treatment ID 
c24_exprs$Gene = gsub("Freeze_", "", c24_exprs$Gene)
c24_exprs$Gene = gsub("Control_", "", c24_exprs$Gene)
  # and anything after the . for now
c24_exprs$Gene = gsub("\\.[0-9]", "", c24_exprs$Gene)


# Check
# Check
c24_exprs %>% 
  filter(Gene == "AT5G61420" & Treatment == "Control") %>%
  ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
  geom_line() +
  geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .1)

# save the df for later plotting
save(c24_exprs, file = "c24_exprs_Ribbons.RData")



### REPEAT WITH TSU
load("/Users/ricon001/Desktop/Git/RNAseq/Arabidopsis/R_output/TSU/TSU_TCs_List.RData")

# Bring in rownames
tsu_exprs = lapply(TCs_List, function(x) {x = rownames_to_column(data.frame(x), var = "Gene");x })

tsu_exprs = do.call(cbind, tsu_exprs)
# Keep one gene column since the rest are duplicated
tsu_exprs =tsu_exprs %>% select(!(c(Freeze_2.Gene, Freeze_3.Gene, Freeze_4.Gene, 
                                    Control_1.Gene, Control_2.Gene, Control_3.Gene, Control_4.Gene)))
# Remove the extra info before the period
colnames(tsu_exprs) = gsub("Freeze_[0-9].", "", colnames(tsu_exprs))
colnames(tsu_exprs) = gsub("Control_[0-9].", "", colnames(tsu_exprs))

# Pivot to have columns to plot on
tsu_exprs = tsu_exprs %>% pivot_longer(cols = !Gene, names_to = "ColID", values_to = "Exprs")
# Separate the ColID
tsu_exprs = tsu_exprs %>% separate(ColID, c("Geno", "Rep", "Treatment", "Zt"))
# Rename tps
tsu_exprs$Zt = str_replace_all(tsu_exprs$Zt, c("ZT17" = "ZT_17", "ZT21" = "ZT_21", "ZT1" = "ZT_1", 
                                               "ZT5" = "ZT_5", "ZT9" = "ZT_9", "ZT13" = "ZT_13"))
levels(factor(tsu_exprs$Zt))

tsu_exprs$Zt <- factor(tsu_exprs$Zt, levels=c("ZT_17", "ZT_21", "ZT_1", 
                                              "ZT_5", "ZT_9", "ZT_13"))

## Calculate average and sd to plot
tsu_exprs = tsu_exprs %>% group_by(Treatment, Zt, Gene) %>%
  mutate(AvgExprs = mean(Exprs), 
         sd = sd(Exprs)) %>% ungroup()

# Remove the Treatment ID 
tsu_exprs$Gene = gsub("Freeze_", "", tsu_exprs$Gene)
tsu_exprs$Gene = gsub("Control_", "", tsu_exprs$Gene)
# and anything after the . for now
tsu_exprs$Gene = gsub("\\.[0-9]", "", tsu_exprs$Gene)


# Check
tsu_exprs %>% 
  filter(Gene == "AT5G61420" & Treatment == "Control") %>%
  ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
  geom_line() +
  geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .1)

# save the df for later plotting
save(tsu_exprs, file = "tsu_exprs_Ribbons.RData")


## Just a subset of genes
  # Katie wants 
AT_selectGenes = list(MYB28 = "AT5G61420", CYP83A1 = "AT4G13770", SUR1 = "AT2G20610", 
                      GTR2 =	"AT5G62680", PEN3 = "AT1G59870")


setwd("~/Desktop/Git/RNAseq")
load("c24_exprs_Ribbons.RData")

pdf("GSL_Control_C24Alone.pdf",width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(AT_selectGenes)) {
  print(c24_exprs %>% 
          filter(Treatment == "Control" & Gene == AT_selectGenes[[i]]) %>%
          ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
          geom_line(color = "olivedrab3", linewidth = 1.2) +
          geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .3, fill = "olivedrab") +
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average expression (control)") +
          ggtitle(paste(names(AT_selectGenes)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()


load("tsu_exprs_Ribbons.RData")

pdf("GSL_Control_tsuAlone.pdf",width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(AT_selectGenes)) {
  print(tsu_exprs %>% 
          filter(Treatment == "Control" & Gene == AT_selectGenes[[i]]) %>%
          ggplot(aes(x = Zt, y = AvgExprs, group = Gene)) +
          geom_line(color = "darkorchid3", linewidth = 1.2) +
          geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd), alpha = .3, fill = "darkorchid4") +
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average expression (control)") +
          ggtitle(paste(names(AT_selectGenes)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()



## Bnapus

#Calculate average exprs and sd
MasterDF = MasterDF %>% group_by(Geno, Treatment, Zts, Gene) %>%
  mutate(AvgExprs = mean(Expression), 
         sd = sd(Expression)) %>% ungroup()

setwd("~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses/Figures")

pdf("EnrichedGSLs.pdf",width = 10,height = 5)
for(i in 1:length(gsl)) {
  print(MasterDF %>% 
          filter(Gene == gsl[[i]] & Geno == St) %>%
          ggplot(aes(x = Zts, y = AvgExprs, group = Gene)) +
          geom_line(aes(color = Treatment, group = interaction(Replicate, Treatment))) +
          geom_ribbon(aes(ymin = AvgExprs - sd, ymax = AvgExprs + sd, group = interaction(Replicate, Treatment), fill = Treatment), alpha = .1) +
          scale_color_manual(values = c("darkgoldenrod1", "chartreuse4")) +
          scale_fill_manual(values = c("darkgoldenrod1", "chartreuse4")) +
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average expression (control)") +
          ggtitle(paste(names(gsl)[i])))
}
dev.off()


## Removing outliers
  # Starting by just removing any data point outside of the whiskers, or the 1.4*IQR
