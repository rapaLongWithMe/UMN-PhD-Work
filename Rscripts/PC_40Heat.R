### Visualizing heat stress data on pennycress mutants
### March 2026 Greenham Lab

## Set up the environment
setwd("~/Desktop/PC_HeatStress/")
require(readxl)
require(tidyverse)


##### ##### ##### ##### ##### ##### ##### 
##### COMPILE ALL DATA AT ONCE ##### ####
##### ##### ##### ##### ##### ##### ##### 

### ~~~ CONTROL DATA 
controlFiles = list.files(path = "ControlData",
                          # Allows you to keep other files in there but will read in ANY csv file
                          pattern = "\\.csv$",
                          # Need this to get a proper path for the next step
                          full.names = TRUE
)

# Read them into R and do some cleaning 
controlData <- lapply(controlFiles, function(x){
  dat <- read_csv(x, skip = 13, show_col_types = FALSE)
  dat = dat[-1, ]
  dat$file <- basename(x) %>% str_remove("\\.csv$")
    # Probs want to change TimeOfDay to ZT
  dat = dat %>% separate(file, into = c("Day", "TimeOfDay", "Treatment"), sep = "_")
  dat
}) %>% bind_rows()

# Check and simplify
colnames(controlData)
ControlDat_short = controlData[, c(1,7:9,11,13,18, 299:301)]

unique(ControlDat_short$Genotype)
unique(ControlDat_short$TimeOfDay)
unique(ControlDat_short$Day)
unique(ControlDat_short$Rep)


### ~~~ HEAT DATA 
heatFiles = list.files(path = "HeatData",
                          # Allows you to keep other files in there but will read in ANY excel file
                          pattern = "\\.csv$",
                          # Need this to get a proper path for the next step
                          full.names = TRUE
)

# Read them into R and do some cleaning 
heatData <- lapply(heatFiles, function(x){
  dat <- read_csv(x, skip = 13, show_col_types = FALSE)
  dat = dat[-1, ]
  dat$file <- basename(x) %>% str_remove("\\.csv$")
  # Probs want to change TimeOfDay to ZT
  dat = dat %>% separate(file, into = c("Day", "TimeOfDay", "Treatment"), sep = "_")
  dat
}) %>% bind_rows()

# Check and simplify
colnames(heatData)
heatDat_short = heatData[, c(1,7:9,11,13,18, 299:301)]


unique(heatDat_short$Genotype)
unique(heatDat_short$TimeOfDay)
unique(heatDat_short$Day)
unique(heatDat_short$Rep)


## Combine the two datasets 
dat = rbind(ControlDat_short, heatDat_short)
# Clean up the labels 
dat = dat %>% unite("Label", Day:Treatment, remove = FALSE)
unique(dat$Label)

dat = dat %>% mutate(Genotype = factor(Genotype, levels = c("mn106 wt", "elf3 1683", "elf3 7", "toc1", "sp32 wt", "hag1hag3", "myc3", "aop2")))
dat = dat %>% mutate(Label = factor(Label, levels = c("Day1_4pm_Control", "Day1_4pm_40C", "Day2_10am_Control", "Day2_10am_40C", 
                                              "Day2_4pm_Control", "Day2_4pm_40C", "Day2_10pm_Control", "Day2_10pm_40C", 
                                              "Day3_4am_Control", "Day3_4am_40C", "Day3_10am_Control", "Day3_10am_40C", "Day3_4pm_Control", "Day3_4pm_40C")))


# Separate into genetic background and type of mutant
dat %>%
  filter(Genotype == "sp32 wt" | Genotype == "hag1hag3") %>%
  ggplot(aes(x = interaction(Genotype, Label), y = as.numeric(gsw), fill = interaction(Genotype, Treatment), color = interaction(Genotype, Treatment))) +
  geom_boxplot(alpha = .5, position = position_dodge(width = .8)) +
  geom_point(position = position_dodge(width = .8)) +
  theme_bw() +
  #facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ## Add vertical lines to indicate heat/no heat
  # on: 10am; off: 10pm 
  geom_vline(xintercept=(8.5), linetype = "dashed", color = "red") +
  geom_vline(xintercept=(16.5), linetype = "dashed", color = "red") +
  geom_vline(xintercept=(24.5), linetype = "dashed", color = "red")


## Add vertical lines to indicate heat/no heat
    # on: 10am; off: 10pm 
# get the coordinates
dat %>%
  filter(Genotype %in% c("sp32 wt", "hag1hag3")) %>%
  with(levels(interaction(Genotype, Label)))

## Print all mutant plots at once: 

pdf("Pennycress_HeatStressPlots.pdf")

## Split by mutant
dat %>%
  filter(Genotype == "sp32 wt" | Genotype == "aop2" | Genotype == "hag1hag3" | Genotype == "myc3") %>%
  ggplot(aes(x = Label, y = as.numeric(gsw), fill = interaction(Genotype, Treatment), color = interaction(Genotype, Treatment))) +
  geom_boxplot(alpha = .5, position = position_dodge(width = .8)) +
  geom_point(position = position_dodge(width = .8)) +
  theme_bw() +
  #facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~Genotype) +
  ggtitle("Stomatal Conductance")

dat %>%
  filter(Genotype == "sp32 wt" | Genotype == "aop2" | Genotype == "hag1hag3" | Genotype == "myc3") %>%
  ggplot(aes(x = Label, y = as.numeric(A), fill = interaction(Genotype, Treatment), color = interaction(Genotype, Treatment))) +
  geom_boxplot(alpha = .5, position = position_dodge(width = .8)) +
  geom_point(position = position_dodge(width = .8)) +
  theme_bw() +
  #facet_grid(~Treatment) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~Genotype) +
  ggtitle("Carbon Assimilation")

dev.off()




##### ##### ##### ##### ##### ##### ##### 
##### #####  ONE AT A TIME ##### ##### ##
##### ##### ##### ##### ##### ##### ##### 

# Bring in the first timepoint
heat1 = read_xlsx("2026-03-10-TP3Heat.xlsx", skip = 14)
  # Get rid of the row with units
heat1 = heat1[-1,]
  # Check 
unique(heat1$Genotype)
unique(heat1$Rep)


# Bring in the second timepoint
heat2 = read_xlsx("2026-03-10-TP2Heat.xlsx", skip = 14)
# Get rid of the row with units
heat2 = heat2[-1,]
# Check 
unique(heat2$Genotype)
unique(heat2$Rep)


# Quick plot
heat1 %>%
  ggplot(aes(x = Genotype, y = as.numeric(gsw))) +
  geom_boxplot()

heat2 %>%
  ggplot(aes(x = Genotype, y = as.numeric(gsw))) +
  geom_boxplot()

# Combine the datasets 
  # Add the timepoint 
heat1$Timepoint = "2"
heat2$Timepoint = "3"
heatDat = rbind(heat1, heat2)

heatDat %>%
  ggplot(aes(x = interaction(factor(Timepoint), factor(Genotype)), y = as.numeric(gsw))) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~Timepoint)



