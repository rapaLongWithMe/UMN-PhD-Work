## Plotting any genes of interest for Katie 
require(tidyverse)
require(ggthemes)

## Set up environment
setwd("~/Desktop/Git/RNAseq/Brapa/For_KT")
FigPath = "./Figures"
Clock_FigPath = "./ClockGene_Comparisons/"
Cold_FigPath = "./ColdGene_Comparisons/"
GSL_FigPath = "~/Desktop/Git/RNAseq/Brapa/GSLgenes_Plots/"

## Set up the dataframes

                ## THIS SECTION USES AVERAGE EXPRESSION ## 
## Bring in master sheet with avg expression for all YS
load("/Users/ricon001/Desktop/Git/RNAseq/Brapa/ExpressedGenes/AverageExprs/ys_expressionAvg_AR.RData")


#### MAKE IT PLOTTABLE ####
ys_expressionAvg = ys_expressionAvg %>% 
  pivot_longer(cols = Cold_1:Warm_6, names_to = "Treat_Zt", values_to = "Avg_Exprs")

levels(factor(ys_expressionAvg$Geno))

#split the pivoted column to make Treatment and Zt columns
ys_expressionAvg = ys_expressionAvg %>% separate(Treat_Zt, c("Treatment", "Zt"))

levels(factor(ys_expressionAvg$Geno))
levels(factor(ys_expressionAvg$Treatment))
levels(factor(ys_expressionAvg$Zt))

#Change the time point numbers to actual Zts
  #have to add the Zt here otherwise some numbers get replaced bc the function is greedy
test = ys_expressionAvg
test$Zt = str_replace_all(ys_expressionAvg$Zt, c("1" = "Zt_17", "2" = "Zt_21", "3" = "Zt_1", 
                        "4" = "Zt_5", "5" = "Zt_9", "6" = "Zt_13"))
levels(factor(test$Zt))
#the 5 gets replaced so fix that 
test$Zt = str_replace_all(test$Zt, c("Zt_Zt_9" = "Zt_5"))

#make sure Zt is factored appropriately (in the right order)
levels(factor(test$Zt))
test$Zt <- factor(test$Zt, levels=c("Zt_17", "Zt_21", "Zt_1", 
                                    "Zt_5", "Zt_9", "Zt_13"))
levels(factor(test$Zt))

#add sd to plot as well



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

            ## THIS SECTION USES ALL REPLICATES ## 
#Make a combined df of all three genos if you haven't already


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 


                  ## DETOXIFICATION OF ELECTROPHILES ## 
# Bring in the data if you haven't already. Just using R500 for now
  ## Bring in the master sheet, or at least expression data
load("/Users/ricon001/Desktop/Git/RNAseq/Brapa/WorkingMaster_DF_long_June15.RData")
Master_DF = test
rm(test)

#make sure Zt is factored appropriately (in the right order)
levels(Master_DF$Zt)

## Make sure you're using AvgExprs for ribbons and not Expression

#print a single gene
pdf(file.path(FigPath,"r500_SuccHydrolase.pdf"),width = 6,height = 4)
Master_DF %>%
  filter(Gene == "BraA05g34240R" & Geno == "R500") %>%
  ggplot(aes(x = Zt, y = AvgExprs)) +
  geom_line(aes(group = interaction(Rep, Treatment), 
                color = Treatment)) +
  geom_ribbon(aes(ymin = AvgExprs - sdExprs, ymax = AvgExprs + sdExprs, 
                  group = interaction(Rep, Treatment)), alpha = .1) +
  scale_color_manual(values = c("blue", "red")) +
  ylab("Average Log2 FPKM") +
  ggtitle("Succinyl hydrolase (R500)") +
  theme_tufte() +
  theme(axis.line = element_line(colour = "grey50"), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_blank(), legend.position = "bottom", legend.title = element_blank(), 
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        title = element_text(size = 16))
dev.off()



                    ## COLD RESPONSIVE GENES ## 

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


#print a single gene to test
test %>%
  filter(Gene == "BraA06g41820R") %>%
  ggplot(aes(x = Zt, y = Avg_Exprs)) +
  geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment))) +
  ylab("Average Log2 FPKM") +
  ggtitle("COR27 (BraA06)") +
  theme_bw() +
  facet_grid(~Treatment)


#need to scale gene expression to compare
test = test %>%
  group_by(Treatment, Geno, Zt) %>% 
  mutate(ScaledCentered = as.numeric(scale(Avg_Exprs)))

test %>%
  filter(Gene == "BraA06g41820R") %>%
  ggplot(aes(x = Zt, y = ScaledCentered)) +
  geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment)), size = 1.2) +
  scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
  ylab("Average Log2 FPKM") +
  ggtitle("COR27 (BraA06)") +
  theme_bw() +
  facet_grid(~Treatment)


## Plotting multiple genes into a single pdf
  #make genes a list to loop through and print
cold = list(cor27_a06, cor27_a09, cor15b_a03_50, cor15b_a03_60, cor28_a08, 
            cbf1_a03, cbf1_a08, cbf2_a03, dreb1_a05, dreb2_a03, ice1_a02, 
            ice1_a06, ice1_a09, camta3_a04, camta3_a09_60, camta3_a09_70, 
            camta5_a08, toc1_a03, toc1_a09)
names(cold) = c("COR27_A06", "COR27_A09","COR15B_A03_50", "COR15B_A03_60", 
                "COR28_A08", "CBF1_A03", "CBF1_A08", "CBF2_A03", "DREB1_A05", 
                "DREB2_A03", "ICE1_A02", "ICE1_A06", "ICE1_A09", "CAMTA3_A04", 
                "CAMTA3_A09_60", "CAMTA3_A09_70", "CAMTA5_A08", "TOC1_A03", "TOC1_A09")

## Plot the 
pdf(file.path(Cold_FigPath,"ColdGenesOnly_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(cold)) {
  print(test %>%
          filter(Gene == cold[[i]] & Treatment == "Cold") %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment)), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
    #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste("Cold", names(cold)[i])) +
    #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
    #move the legend around, = "none" to remove
          legend.position = "bottom"))
}
dev.off()


pdf(file.path(Cold_FigPath,"ControlVersionOnly_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(cold)) {
  print(test %>%
          filter(Gene == cold[[i]] & Treatment == "Warm") %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment)), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste("Warm", names(cold)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()


#try one with all six lines on one plot
pdf(file.path(Cold_FigPath,"ColdAndControl_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(cold)) {
  print(test %>%
          filter(Gene == cold[[i]]) %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment), linetype = Treatment), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste(names(cold)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()


                      ## CLOCK GENES ## 
  ## pulling from my "core clock gene" list in the RNAseq folder
elf3_A04 = "BraA04g18480R"
elf3_A09_20 = "BraA09g50820R"
elf3_A09_30 = "BraA09g50830R"

gi = "BraA09g38670R"
lhy = "BraA10g01800R"
lnk1 = "BraA06g25990R"
lnk2_A04 = "BraA04g06240R"
lnk2_A07 = "BraA07g20200R"
lnk2_A09 = "BraA09g43520"

cca1 = "BraA05g01930R"
elf4_A03 = "BraA03g20890R"
elf4_A04 = "BraA04g27180R"
elf4_A05 = "BraA05g06590R"

lux = "BraA06g19700R"

prr1_A03 = "BraA03g42600R"
prr1_A09 = "BraA09g07570R"

prr5_A02 = "BraA02g43100R"
prr5_A06 = "BraA06g29990R"
prr5_A09 = "BraA09g06770R"

prr7_A02 = "BraA02g01670R"
prr7_A10 = "BraA10g31240R" 

prr9_A02 = "BraA04g31570R"
prr9_A10 = "BraA05000188R"

rve4_A03 = "BraA03g01770R"
rve4_A09 = "BraA09g45590R"
rve4_A10 = "BraA10g31210R"

rve8_A05 = "BraA05g35220R"
rve8_A09 = "BraA09g56570R"

#put into a list to loop through while plotting
clock = list(elf3_A04,elf3_A09_20,elf3_A09_30,gi,lhy,lnk1,lnk2_A04,
             lnk2_A07,lnk2_A09,cca1,elf4_A03,elf4_A04,elf4_A05,lux,
             prr1_A03,prr1_A09,prr5_A02,prr5_A06,prr5_A09,prr7_A02,
             prr7_A10,prr9_A02,prr9_A10,rve4_A03,rve4_A09,rve4_A10,
             rve8_A05,rve8_A09)
names(clock) = c("ELF3_A04","ELF3_A09_20","ELF3_A09_30","GI", "LHY","LNK1","LNK2_A04",
             "LNK2_A07","LNK2_A09","CCA1","ELF4_A03","ELF4_A04","ELF4_A05","LUX",
             "PRR1_A03","PRR1_A09","PRR5_A02","PRR5_A06","PRR5_A09","PRR7_A02",
             "PRR7_A10","PRR9_A04","PRR9_A05","RVE4_A03","RVE4_A09","RVE4_A10",
             "RVE8_A05","RVE8_A09")


pdf(file.path(Clock_FigPath,"ColdAndControl_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(clock)) {
  print(test %>%
          filter(Gene == clock[[i]]) %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment), linetype = Treatment), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste(names(clock)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()

## Plot individual panels per clock gene by each genotype
pdf("~/Desktop/Git/RNAseq/Brapa/ClockGene_Plots/r500_AllCoreClock.pdf",width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(clock)) {
  print(Master_DF %>%
          filter(Gene == clock[[i]] & Geno == "R500") %>%
          ggplot(aes(x = Zt, y = AvgExprs, fill = Treatment)) +
          geom_line(aes(x = Zt, color = Treatment, group = interaction(Geno, Treatment))) +
          geom_ribbon(aes(ymax = AvgExprs + sdExprs, ymin = AvgExprs - sdExprs, 
                          group = interaction(Rep, Treatment)), alpha = .1) +
          scale_color_manual(values = c("blue", "red")) +
          scale_fill_manual(values = c("blue", "red")) +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste(names(clock)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()

pdf("~/Desktop/Git/RNAseq/Brapa/ClockGene_Plots/acc28_AllCoreClock.pdf",width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(clock)) {
  print(Master_DF %>%
          filter(Gene == clock[[i]] & Geno == "ACC28") %>%
          ggplot(aes(x = Zt, y = AvgExprs, fill = Treatment)) +
          geom_line(aes(x = Zt, color = Treatment, group = interaction(Geno, Treatment))) +
          geom_ribbon(aes(ymax = AvgExprs + sdExprs, ymin = AvgExprs - sdExprs, 
                          group = interaction(Rep, Treatment)), alpha = .1) +
          scale_color_manual(values = c("blue", "red")) +
          scale_fill_manual(values = c("blue", "red")) +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste(names(clock)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()


pdf(file.path(Clock_FigPath,"Cold_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(clock)) {
  print(test %>%
          filter(Gene == clock[[i]] & Treatment == "Cold") %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment), linetype = Treatment), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste("Cold", names(clock)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()

pdf(file.path(Clock_FigPath,"Control_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(clock)) {
  print(test %>%
          filter(Gene == clock[[i]] & Treatment == "Warm") %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment), linetype = Treatment), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste("Control", names(clock)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()



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
sur1_A09 = "BraA09g12830R"

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


pdf(file.path(GSL_FigPath,"GSL_ColdAndControl_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(gsl)) {
  print(test %>%
          filter(Gene == gsl[[i]]) %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment), linetype = Treatment), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste(names(gsl)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()

pdf(file.path(GSL_FigPath,"GSL_Coldl_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(gsl)) {
  print(test %>%
          filter(Gene == gsl[[i]] & Treatment == "Cold") %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment), linetype = Treatment), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste("Cold",names(gsl)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()


pdf(file.path(GSL_FigPath,"GSL_Control_YS.pdf"),width = 10,height = 5)
### ### ### ### ### ### ### ###
for(i in 1:length(gsl)) {
  print(test %>%
          filter(Gene == gsl[[i]] & Treatment == "Warm") %>%
          ggplot(aes(x = Zt, y = ScaledCentered)) +
          geom_line(aes(x = Zt, color = Geno, group = interaction(Geno, Treatment), linetype = Treatment), size = 1.2) +
          scale_color_manual(values = c("darkolivegreen", "burlywood3", "deepskyblue3"), name = "Genotype") +
          #base size is the font size
          theme_classic(base_size = 18) +
          xlab("Zts (hours after lights on)") +
          ylab("Average Log2 FPKM") +
          ggtitle(paste("Warm", names(gsl)[i])) +
          #axis.text.x changes the font of the xticks
          theme(axis.text.x=element_text(size=10), 
                #move the legend around, = "none" to remove
                legend.position = "bottom"))
}
dev.off()

