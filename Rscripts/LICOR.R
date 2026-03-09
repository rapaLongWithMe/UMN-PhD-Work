### First look at LICOR data

## Set up the enviroment
require(tidyverse)
require(readxl)

setwd("~/Desktop/Git/Bnapus_drought/Physiology/RoundThree/LiCOR")
FigPath = "./Figures"


## Bring in and clean the data
  # Most of the cleaning was done in excel 
licor = read_xlsx("LicorData_ForR.xlsx", sheet = 1)

levels(factor(licor$Geno))
levels(factor(licor$Treatment))
levels(factor(licor$ExperimentalDay))
levels(factor(licor$ZT))


## Potential colors for palettes; darker is the next line and the lighter version below that 
##595B48 #C55E44 #5E0B15 #5296A5 #04007A
##A3A58D #D99382 #F6B6BD #9FC7D0 #B0ADFF


## A quick look
# All genos as panels; change the variable of interest
licor %>% 
  filter(Treatment == "WW" | Treatment == "WL") %>%
  ggplot(aes(x = ExperimentalDay, y = gsw, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
 # turn off outlier shape to remove the extra points produced by the default
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("gsw") +
  xlab("") +
  theme(legend.position = "bottom") +
  facet_wrap(~Geno)
###

#Plot each geno's phys to a separate pdf; combine after
pdf(file.path(FigPath, file = "20230811_StPhys_UnderWL.pdf"), width = 6, height = 4)
# A
licor %>% 
  filter(Geno == "St") %>%
  ggplot(aes(x = ExperimentalDay, y = A, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("A (umols/m2/s)") +
  ggtitle("St:CO2 Assimilation") +
  xlab("") +
  theme(legend.position = "bottom")

# gsw 
licor %>% 
  filter(Geno == "St") %>%
  ggplot(aes(x = ExperimentalDay, y = gsw, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("gsw (mols/m2/s)") +
  ggtitle("St: Stomatal conductance") +
  xlab("") +
  theme(legend.position = "bottom")


# E
licor %>% 
  filter(Geno == "St") %>%
  ggplot(aes(x = ExperimentalDay, y = E, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("E (mols/m2/s)") +
  ggtitle("St: Transpiration rate (E) ") +
  xlab("") +
  theme(legend.position = "bottom")
dev.off()
###


###
pdf(file.path(FigPath, file = "20230811_MuPhys_UnderWL.pdf"), width = 6, height = 4)
# A
licor %>% 
  filter(Geno == "Mu") %>%
  ggplot(aes(x = ExperimentalDay, y = A, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("A (umols/m2/s)") +
  ggtitle("Mu:CO2 Assimilation") +
  xlab("") +
  theme(legend.position = "bottom")

# gsw 
licor %>% 
  filter(Geno == "Mu") %>%
  ggplot(aes(x = ExperimentalDay, y = gsw, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("gsw (mols/m2/s)") +
  ggtitle("Mu: Stomatal conductance") +
  xlab("") +
  theme(legend.position = "bottom")


# E
licor %>% 
  filter(Geno == "Mu") %>%
  ggplot(aes(x = ExperimentalDay, y = E, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("E (mols/m2/s)") +
  ggtitle("Mu: Transpiration rate (E) ") +
  xlab("") +
  theme(legend.position = "bottom")
dev.off()
###

###
pdf(file.path(FigPath, file = "20230811_AbPhys_UnderWL.pdf"), width = 6, height = 4)
# A
licor %>% 
  filter(Geno == "Ab") %>%
  ggplot(aes(x = ExperimentalDay, y = A, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("A (umols/m2/s)") +
  ggtitle("Ab:CO2 Assimilation") +
  xlab("") +
  theme(legend.position = "bottom")

# gsw 
licor %>% 
  filter(Geno == "Ab") %>%
  ggplot(aes(x = ExperimentalDay, y = gsw, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("gsw (mols/m2/s)") +
  ggtitle("Ab: Stomatal conductance") +
  xlab("") +
  theme(legend.position = "bottom")


# E
licor %>% 
  filter(Geno == "Ab") %>%
  ggplot(aes(x = ExperimentalDay, y = E, color = Treatment, group = interaction(Treatment, ExperimentalDay, ZT))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ylab("E (mols/m2/s)") +
  ggtitle("Ab: Transpiration rate (E) ") +
  xlab("") +
  theme(legend.position = "bottom")
dev.off()
###



# Connected replicates
pdf(file.path(FigPath, file = "PhotoAssim_TwoWksAfter.pdf"), width = 6, height = 4)
licor %>% 
  filter(Geno == "St" & date >= 20230522) %>%
  ggplot(aes(x = factor(Zt), y = A, color = Treatment, group = interaction(Treatment, LicorID))) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("St: Two weeks after 20% FC") +
## Add this if want to see which rep is doing what
  # require(ggrepel)
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") 

licor %>% 
  filter(Geno == "Mu" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = A, color = Treatment, group = interaction(Treatment, Replicate))) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Mu: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") 


licor %>% 
  filter(Geno == "Ab" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = A, color = Treatment, group = interaction(Treatment, Replicate))) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Ab: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") 

dev.off()


pdf(file.path(FigPath, file = "gsw_TwoWksAfter.pdf"), width = 6, height = 4)
licor %>% 
  filter(Geno == "St" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = gsw, color = Treatment, group = interaction(Treatment, Replicate))) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("St: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") 

licor %>% 
  filter(Geno == "Mu" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = gsw, color = Treatment, group = interaction(Treatment, Replicate))) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Mu: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") 


licor %>% 
  filter(Geno == "Ab" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = gsw, color = Treatment, group = interaction(Treatment, Replicate))) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Ab: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") 

dev.off()


## Plot a little more seriously...

# Boxplots; check out the spread
pdf(file.path(FigPath, file = "BoxPlots_PhotoAssim_TwoWksAfter.pdf"), width = 6, height = 4)
licor %>% 
  filter(Geno == "St" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = A, color = Treatment)) +
  geom_boxplot() + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("St: Two weeks after 20% FC") +
  xlab("") +
  theme(legend.position = "bottom")

licor %>% 
  filter(Geno == "Mu" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = A, color = Treatment)) +
  geom_boxplot() + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Mu: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") +
  theme(legend.position = "bottom")


licor %>% 
  filter(Geno == "Ab" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = A, color = Treatment)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Ab: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") +
  theme(legend.position = "bottom")

dev.off()


pdf(file.path(FigPath, file = "BoxPlots_gsw_TwoWksAfter.pdf"), width = 6, height = 4)
licor %>% 
  filter(Geno == "St" & date >= 20230522) %>%
  ggplot(aes(x = factor(Zt), y = gsw, color = Treatment)) +
  geom_boxplot() + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("St: Two weeks after 20% FC") +
  xlab("") +
  theme(legend.position = "bottom") +
  facet_wrap(~ExperimentalDay)

licor %>% 
  filter(Geno == "Mu" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = gsw, color = Treatment)) +
  geom_boxplot() + 
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Mu: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") +
  theme(legend.position = "bottom")


licor %>% 
  filter(Geno == "Ab" & date >= 20230522) %>%
  ggplot(aes(x = factor(ZT), y = gsw, color = Treatment)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = .2), alpha = .6) +
  scale_color_manual(values = c("darkgoldenrod2", "chartreuse4")) +
  theme_bw() +
  ggtitle("Ab: Two weeks after 20% FC") +
  ## Add this if want to see which rep is doing what
  #geom_label_repel(aes(label = Replicate), max.overlaps = 40) +
  xlab("") +
  theme(legend.position = "bottom")

dev.off()



  

