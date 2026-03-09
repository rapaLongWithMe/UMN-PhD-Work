## Compile all GO enrichment files into a single master list

## Set Up
setwd("~/Desktop/Git/RNAseq/Brapa/EileenKosola/GO_Enrichment/acc28")
require(tidyverse)

## Bring in the files and combine into a single dataframe
  #create a vector that contains the name of files to loop over
filenames <- list.files(pattern = "*.csv", full.names = FALSE)
filenames

##### CHANGE BASED ON GENOTYPE #####
acc28_go = lapply(filenames, function(x) read_csv(x) %>% 
  #if you don't bind you'll end up with a tibble per csv file
  bind_rows %>%
  mutate(file = x, 
         Geno = "acc28")) 
  
acc28_go = do.call("rbind", lapply(acc28_go, as.data.frame))
colnames(acc28_go)[1] = "GO_term"

levels(as.factor(acc28_go$Geno))
levels(as.factor(acc28_go$file))


## Also interested in compiling a yellow sarson (acc28, acc50 and r500) GO sheet 
  ## repeat for the other genotypes
filenames <- list.files(pattern = "*.csv", full.names = FALSE)
filenames

##### CHANGE BASED ON GENOTYPE #####
acc50_go = lapply(filenames, function(x) read_csv(x) %>% 
                  #if you don't bind you'll end up with a tibble per csv file
                  bind_rows %>%
                  mutate(file = x, 
                         Geno = "acc50")) 

acc50_go = do.call("rbind", lapply(acc50_go, as.data.frame))
colnames(acc50_go)[1] = "GO_term"

levels(as.factor(acc50_go$Geno))
levels(as.factor(acc50_go$file))


##### CHANGE BASED ON GENOTYPE #####
r500_go = lapply(filenames, function(x) read_csv(x) %>% 
                    #if you don't bind you'll end up with a tibble per csv file
                    bind_rows %>%
                    mutate(file = x, 
                           Geno = "r500")) 

r500_go = do.call("rbind", lapply(r500_go, as.data.frame))
colnames(r500_go)[1] = "GO_term"

levels(as.factor(r500_go$Geno))
levels(as.factor(r500_go$file))


## Bind the dfs together into one ys Master sheet
ys_GO = rbind(acc28_go, acc50_go, r500_go)
levels(factor(ys_GO$Geno))
levels(factor(ys_GO$file))


## save the master file
write.csv(ys_GO, file.path("~/Desktop/Git/RNAseq/Brapa/EileenKosola/GO_Enrichment", file = "ys_GO.csv"))

## Now want a df that contains just the overlapping terms between all three genos
  #First pull out only sig. enriched pathways that are over enriched 

sig_ys_GO = ys_GO %>%
  filter(`Adj_P-Value` <= 0.01 & `Enrichment Direction` == "Over")

#note we're now down to t phase groups
levels(factor(sig_ys_GO$file))
levels(factor(sig_ys_GO$Geno))

#split the GO_Term column to keep the descriptions separate
terms = str_split(sig_ys_GO$GO_term, " - ", simplify = TRUE)
#bind in these terms to the master sheet
sig_ys_GO = cbind(sig_ys_GO, terms)
#peel off the original column
sig_ys_GO = sig_ys_GO[, -1]
#rename the new columns
glimpse(sig_ys_GO)
colnames(sig_ys_GO)[12] = "GO_Term"
colnames(sig_ys_GO)[13] = "GO_Description"

#save this as well 
getwd()
write.csv(sig_ys_GO, file.path("~/Desktop/Git/RNAseq/Brapa/EileenKosola/GO_Enrichment", file = "OverSig_ys_GO.csv"))


## We can now explore the master data sheet for pathways of interest
  #read it in, if starting here
MasterSheet = "~/Desktop/Git/RNAseq/Brapa/EileenKosola/GO_Enrichment/Master_ys_GO.csv"
go = read.csv(MasterSheet)
  #remove the column of rownumbers
go = go[, -1]

#the files are named by the phase group that they were enriched in
levels(factor(go$file))

#split the GO_Term column to keep the descriptions separate
terms = str_split(go$GO_term, " - ", simplify = TRUE)
#bind in these terms to the master sheet
go = cbind(go, terms)
#peel off the original column
go = go[, -1]
#rename the new columns
glimpse(go)
colnames(go)[12] = "GO_Term"
colnames(go)[13] = "GO_Description"

#filter for just sig terms
go = go %>% filter(Adj_P.Value <= 0.1)


gsl = go %>% filter(str_detect(go$GO_Description, "glucosinolate"))
glycosin = go %>% filter(str_detect(go$GO_Description, "glycosinolate"))
#bind together
gsl = rbind(gsl, glycosin)

write.csv(gsl, file = "SigenrichedGSLs_YSphaseGroups.csv", row.names = FALSE)
write.csv(go, file = "YS_SigGOterms.1.csv", row.names = FALSE)
