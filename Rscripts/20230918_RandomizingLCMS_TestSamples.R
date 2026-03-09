#Randomizing samples for LC-MS Runs

### PREVIOUSLY 
# simplest way is to randomly pick numbers that represent the tube number 
# sample(78:156,78, replace=FALSE)
# can also add pools and concatenate both into a string to randomly sample from
samps = seq(5, 28)
pools = c("Pool_Inj.5", "Pool_Inj.1", "Pool_Inj.25", "GP.1", "GP.25", "GP.5")
template = sample(c(samps, pools), replace=FALSE)

# decided to add another dilution to the 1:2 and 1:10 dilutions, about halfway between these 
addOns = seq(26, 34)
template2 = sample(c(template, addOns), replace = FALSE)
  
  
# Write out as a csv to copy and paste into sequence template
setwd("~/Desktop/Git/Metabolomics/2022_TimeSeries/CircadianTC_Summer2022")
write.csv(template, "20230918_RandomizedTestSamples.csv", row.names = FALSE)


