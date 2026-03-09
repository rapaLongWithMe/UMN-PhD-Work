#plotting gene expression for genes of interest
require(DiPALM)


#example code to start with#
  #napus data, drought palette and three reps
pdf(file.path(FigPath,"ACC28_CCA1_rapa.pdf"),width = 4,height = 3)
PlotTCs(TClst = TCs_List,tgene = r_cca1.1, 
        scale = T, xlab="ZT Time (hours)", xAxsLabs = c(seq(1,23,4),seq(1,23,4)), 
        ledgeX="topright",
        tcols = c("darkgoldenrod1","darkgoldenrod1", "darkgoldenrod1",
                  "chartreuse4", "chartreuse4", "chartreuse4"), 
        tltys = c(1,2, 3,1,2, 3), 
        main = "B. rapa: cca1")
dev.off()

##name your genes of interest

#clock genes
#rapa 
r_cca1.1 = "BraA05g01930R"
r_toc1.1 = "BraA03g42600R"
r_toc1.2 = "BraA09g07570R"
r_prr5.1 = "BraA02g43100R"
r_prr5.2 = "BraA06g29990R"
r_prr5.3 = "BraA09g06770R"
r_lhy1.1 = "BraA10g01800R"


#r_cor27.1 = "BraA06g41820R"
#r_cor27.2 = "BraA09g19760R"

o_cca1.1 = "gene:Bo4g006930"
o_toc1.1 = "gene:Bo7g098570"
o_toc1.2 = "gene:Bo9g014570"
o_prr5.1 = "gene:Bo7g096790"
o_prr5.2 = "gene:Bo9g012630"
# o_prr5.3 lost


#light harvesting PSI
lhca = "BraA08g27640R"

#napus drought genes
r_dreb1.1 = "BraA09g09840R"
r_dreb1.2 = "BraA07g15580R"
o_dreb1.1 = "Bo7g063980"
o_dreb1.2 = "Bo9g020380"


## Compile your genes of interest into a named list

rapa_clockGenes = list(r_cca1.1, r_toc1.1, r_toc1.2)
names(rapa_clockGenes) = c("cca1.1", "toc1.1", "toc1.2")
o_clockGenes = list(o_cca1.1, o_toc1.1, o_toc1.2)
names(o_clockGenes) = c("cca1.1", "toc1.1", "toc1.2")

rapa_prrs = list(r_prr5.1, r_prr5.2, r_prr5.3)
o_prrs = list(o_prr5.1, o_prr5.2)

#napus drought 
rapa_drought = list(r_dreb1.1, r_dreb1.2)
o_drought = list(o_dreb1.1, o_dreb1.2)


## Load the expression data in list form. Need each named element to represent the entire timeseries for a replicate (by treatment)
load(file.path("~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses/ExpressedGenes", file = "AllGenos_exprs_List.RData"))

#or for a single geno
st = exprs$St



##Plot to a pre-definited folder
#clock genes
ClockPath = "~/Desktop/Git/RNAseq/Bnapus/ExpressedGenes/ClockGenes"
  #work around to print out the subset of genos that have some gene expressed when others don't
toc1.2Genos = list(exprs[[1]], exprs[[10]], exprs[[12]])
names(toc1.2Genos) = c("Ab", "Mu", "Se")

#this will print the first element of the gene list but for all genos
pdf(file.path(ClockPath,"checking_toc1.2_oleracea.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
          tgenes = o_toc1.2,
          scale = T, 
          tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                    "chartreuse4","chartreuse4", "chartreuse4"), 
          tltys = c(1,2, 3, 1,2, 3), 
          main = paste("TOC1.2:", 
                       names(exprs)[i], sep=" "))
}
dev.off()

#drebs
FigPath = "~/Desktop/Git/Bnapus_drought/Bnapus_TranscriptomeAnalyses"
pdf(file.path(FigPath,"DREB: B. rapa Copy 1.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = rapa_drought[[1]],
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "B. rapa Copy 1")
}
dev.off()

pdf(file.path(FigPath,"DREB: B. oleracea Copy 1.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = o_drought[[1]],
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "B. oleracea Copy 1")
}
dev.off()

pdf(file.path(FigPath,"DREB: B. rapa Copy 2.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = rapa_drought[[2]],
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "B. rapa Copy 2")
}
dev.off()

pdf(file.path(FigPath,"DREB: B. oleracea Copy 2.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = o_drought[[2]],
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "B. rapa Copy 2")
}
dev.off()

#prrs
pdf(file.path(FigPath,"prr5: B. rapa Copy 2.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = r_prr5.2,
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "prr5: B. rapa Copy 2")
}
dev.off()

pdf(file.path(FigPath,"prr5: B. rapa Copy 3.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = r_prr5.3,
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "prr5: B. rapa Copy 3")
}
dev.off()


pdf(file.path(FigPath,"prr5: B. oleracea Copy 1.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = o_prr5.1,
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "prr5: B. oleracea Copy 1")
}
dev.off()

pdf(file.path(FigPath,"prr5: B. oleracea Copy 2.pdf"),width = 10,height = 5)
for(i in 1:length(exprs)) {
  PlotTCsRibbon(TClst = exprs[[i]], 
                tgenes = o_prr5.2,
                scale = T, 
                tcols = c("darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1",
                          "chartreuse4","chartreuse4", "chartreuse4"), 
                tltys = c(1,2, 3, 1,2, 3), 
                main = "prr5: B. oleracea Copy 2")
}
dev.off()







#oleracea
o_cca1.1 = "gene:Bo4g006930"
o_toc1.1 = "gene:Bo7g098570"

#####
#GSL biosynthesis genes
#rapa
r_myb28.1 = "BraA09g07610R"
r_myb28.2 = "BraA02g44060R"
r_myb28.3 = "BraA03g42560R"

o_myb28.1 = "gene:Bo7g098590"
o_myb28.2 = "gene:Bo2g161590"
o_myb28.3 = "gene:Bo9g014610"

r_myb29.1 = "BraA10g28030R"   #subscript out of bounds? 
r_myb29.2 = "BraA03g03990R"   #subscript out of bounds? 
r_myb29.3 = "BraA03g04020R"   #subscript out of bounds? 
r_myb29.4 = "BraA03g04050R"   #subscript out of bounds? 

o_myb29.1 = "gene:Bo9g175680"   #subscript out of bounds? 
o_myb29.2 = "gene:Bo3g004500"   #subscript out of bounds? 

#####
#sulfur related
#rapa
r_paps.1 = "BraA10g30900R"
r_paps.2 = "BraA10g30910R"
r_paps.3 = "BraA02g01920R"
r_paps.4 = "BraA03g01990R"
r_paps.5 = "BraA03g02030R"

#oleracea
o_paps.1 = "gene:Bo3g002340"
o_paps.2 = "gene:Bo2g004070"

#####
#drought 
#rapa
r_dreb19 = "BraA05g07750R"  #subscript out of bounds? 

#oleracea
o_dreb19 = "gene:Bo4g027470"   #subscript out of bounds? 




##### CHANGE THE FILE NAME BASED ON GENOTYPE #####
#pdf(file.path(FigPath,"DH12_CCA1_rapa.pdf"),width = 4,height = 3)

r_elf4 = Master_DF %>%
  filter(Clock == "ELF4" & Geno == "St") %>%
  select(Clock)

PlotTCs(TClst = TCs_List,tgene = r_myb29.1, 
        scale = T, xlab="ZT Time (hours)", xAxsLabs = c(seq(1,23,4),seq(1,23,4)), 
        ledgeX="topright",
        tcols = c("darkgoldenrod1","darkgoldenrod1", "darkgoldenrod1",
                  "chartreuse4", "chartreuse4", "chartreuse4"), 
        tltys = c(1,2, 3,1,2, 3), 
        main = "B. rapa: myb29")

#dev.off()


