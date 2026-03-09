## Analyzing soil moisture data (Volumetric Water Content) and water usage for Brassica napus


## Set up Environment
setwd("~/Desktop/Git/Bnapus_drought/Physiology")

require(tidyverse)
require(readxl)


## Bring in and clean the data
dat = read_xlsx("CompiledSoilMoisture.xlsx")


levels(factor(dat$Treatment))
levels(factor(dat$Genotype))
# Removing "extra" pots
dat = dat %>% filter(Treatment != "EXTRA")


##### ##### ##### ##### 
# Not working? 
test = dat %>% 
  mutate(
    across(everything(), str_replace_na(.))
  )
##### ##### ##### ##### 



## Look at the relationship between VWC and Pre water weight

### Becomes more linear as the drought progresses, which is good 
dat %>% 
  ggplot(aes(y = VWC_PreWater, x = PreWaterWeight, color = Treatment)) +
  geom_point(alpha = .9) +
  scale_color_manual(values = c("goldenrod2", "chartreuse4")) +
  ylim(10, 40) +
  facet_grid(Genotype ~ TreatmentDay) +
  theme_bw() +
  ylab("Volumeric Water Content") +
  xlab("Pre-watering weight") +
  ggtitle("Days after drought")

dat %>% 
  ggplot(aes(y = VWC_PreWater, x = PreWaterWeight, fill = Treatment)) +
  geom_boxplot(alpha = .6) +
  geom_point(aes(color = Treatment), alpha = .3) +
  scale_fill_manual(values = c("goldenrod2", "chartreuse4")) +
  scale_color_manual(values = c("goldenrod2", "chartreuse4")) +
  ylim(10, 40) +
  facet_grid(Genotype ~ TreatmentDay) +
  theme_bw() +
  ylab("Volumeric Water Content") +
  xlab("Pre-watering weight") +
  ggtitle("Days after drought")


## Separate the genos into individual plots
dat %>% 
  filter(Genotype == "St") %>% 
  ggplot(aes(y = VWC_PreWater, x = PreWaterWeight, fill = Treatment)) +
  geom_boxplot(alpha = .6) +
  geom_point(aes(color = Treatment), alpha = .3) +
  scale_fill_manual(values = c("goldenrod2", "chartreuse4")) +
  scale_color_manual(values = c("goldenrod2", "chartreuse4")) +
  ylim(10, 40) +
  facet_grid(~TreatmentDay) +
  theme_bw() +
  ylab("Volumeric Water Content") +
  xlab("Pre-watering weight") +
  ggtitle("Days after drought")

## Why does St only have one point (not a boxplot) on day7 in drought? 
dat %>% filter(Genotype == "St" & TreatmentDay == 7 & Treatment == "Drought") %>% pull(VWC_PreWater)

dat %>% filter(Genotype == "St" & TreatmentDay == 7 & Treatment == "Drought") %>%
  ggplot(aes(y = VWC_PreWater, x = PreWaterWeight, fill = Treatment)) +
  geom_boxplot(alpha = .6)




## Quick look at water to add and VWC
  ### Will have to pivot the table 


