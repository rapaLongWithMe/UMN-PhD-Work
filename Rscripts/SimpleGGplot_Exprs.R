## Plotting one or two genes of interest (expression)

## Set up the environment
require(tidyverse)
require(ggthemes)

setwd("~/Desktop/Git/RNAseq/Brapa/GRNs/2022_YellowSarsons")
FigPath = "~/Desktop/Git/RNAseq/Brapa/GRNs/2022_YellowSarsons/Figures"

## Bring in the master sheet, or at least expression data
load("/Users/ricon001/Desktop/Git/Clean_Rscripts/ACC28_PlottingDF.RData")
  # Here I am renaming the 'test' DF I was working on at the time
ACC28_PlottingDF = ACC28_PlottingDF %>% filter(Geno == "ACC28")


## Plotting single genes
#geneName_chromosome

#Cold responsive genes
## COR27: AT5G42900
cor27_a06 = "BraA06g41820R"
cor27_a09 = "BraA09g19760R"
## cor15B:AT2G42530
cor15b_a03_50 = "BraA03g21950R"
cor15b_a03_60 = "BraA03g21960R"
## cor28: AT4G33980
cor28_a08 = "BraA08g16300R"
## cbf1: AT4G25490
cbf1_a08 = "BraA08g19970R"
cbf1_a03 = "BraA03g15690R"
## cbf2: AT4G25470
cbf2_a03 = "BraA03g51570R"
## dreb1: AT2G38340
dreb1_a05 = "BraA05g07750R"
## dreb2: AT4G25470
dreb2_a03 = "BraA03g51570R"
## ice1: AT3G26744
ice1_a06 = "BraA06g36240R"
ice1_a09 = "BraA09g03540R"
ice1_a02 = "BraA02g38130R"
## camta5: AT4G16150
camta5_a08 = "BraA08g11820R"
## camta3: AT2G22300
camta3_a09_60 = "BraA09g52860R"
camta3_a09_70 = "BraA09g52870R	"
camta3_a04 = "BraA04g15660R"
#toc1 is cold and clock
## toc1: AT5G61380	
toc1_a03 = "BraA03g42600R"	
toc1_a09 = "BraA09g07570R"



# Make sure Zt is factored appropriately (in the right order)
levels(ACC28_PlottingDF$ZT)
  # You may or may not need this depending on if your ZT column is organized properly or not
  # This column needs to be set as a factor; you can either do it here or in the ggplot call itself
#ACC28_PlottingDF$ZT <- factor(ACC28_PlottingDF$ZT, levels=c(1, 5, 9, 13, 17, 21))

# Print a single gene
  # The next line ('pdf.....') will print your plot directly into a new pdf file. 
  # If you just want to do a quick look in R you can delete that line and also delete 'dev.off()' 
pdf(file.path(FigPath,"acc28_cor27_a06.pdf"),width = 6,height = 4)
ACC28_PlottingDF %>%
  # Here you replace "Bra...." with either what you named your gene or the annotation in quotes
    # I am filtering my MasterDF to plot this gene for a single genotype....you may need to do the same if you have multiple genos in the same DF
  filter(Gene == cor27_a06 & Geno == "ACC28") %>%
  ggplot(aes(x = ZT, y = AvgExprs)) +
  # If you don't have reps you will have to change the group, probably to Treatment (depends on how your DF is set up)
  geom_line(aes(x = ZT, group = interaction(Rep, Treatment), color = Treatment)) +
  geom_ribbon(aes(ymin = Expression - sdExprs, ymax = Expression + sdExprs, 
                  group = interaction(Rep, Treatment)), alpha = .2) +
  scale_color_manual(values = c("blue", "red")) +
  ylab("Average Log2 FPKM") +
  # Change the title to something informative for your gene; the Bra is usually helpful
  ggtitle("ACC28 (COR27; A06)") +
  # you can delete the rest of this, it's just customization 
  theme_tufte() +
  theme(axis.line = element_line(colour = "grey50"), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        title = element_text(size = 16))
dev.off()

