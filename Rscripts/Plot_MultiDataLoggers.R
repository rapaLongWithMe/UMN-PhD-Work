### Plotting datalogger data from many files at once 

## Set up your environment
# The require function is less verbose; the first time you install packages you will have to use install.packages and the library functions
require(tidyverse)
require(readr)
require(lubridate)
require(ggtext)

setwd("~/Desktop/Git/Bnapus_drought/Physiology/RoundThree/DataLoggers")

## Bring in the data and add the file name to separate
all.files = list.files(pattern = ".txt")
all.files

# I've decided to ignore my first few which were just getting the chamber set up originally 
all.files = all.files[6:17]


# Make a function to process each file
processFile <- function(f) {
  # Turn off the check.names function so it ignores duplicate column names
  df <- read.csv(paste0(f), check.names = FALSE) %>% bind_rows %>% mutate(file = f)
  # Clip off the file extension
  df$file = gsub(".txt", "", df$file)
  # If you name your files with metadata seprated by `_` you can now pull that information into new columns
    # add `remove = FALSE` if you want to keep the original file ID
  df = df %>% separate(file, c("WeekOfExperiment", "Date", "Logger"))
  # Simplfy the df
  df = df[,-c(1, 5,6)]
}


dat = lapply(all.files, processFile)
# Bind together into one large df
all_dat = do.call("rbind", dat)
# Rename columns with UTF8 characters otherwise you'll get warning messages the entire session
colnames(all_dat)[2] = "Celcius"
colnames(all_dat)[3] = "Humidity"

# Make the time column an actual date
all_dat = all_dat %>% separate(Time, c("Date", "Time"), sep = " ")
all_dat = all_dat %>% mutate(Date = as.Date(Date))

## Average the environmental values by hour to smooth out the plot
all_dat = all_dat %>% separate(Time, c("Hour", "Minute", "Second"), remove = FALSE)
all_dat = all_dat %>% group_by(Date, Logger, Hour) %>% mutate(AvgTemp = mean(Celcius)) %>% ungroup()

# Make a new xaxis to plot one point every hour over the entire time course
all_dat = all_dat %>% separate(as.numeric(Date), c("Year", "Month", "Day"))


## Add the new x-axis to the df
xvar = paste(all_dat$Month, all_dat$Day, all_dat$Hour, sep = "_") 
xvar = paste(all_dat$Month, all_dat$Day, sep = "_") 
all_dat$xvar = xvar


## Make yourself a color palette; type `help(colorRampPalette` in the console for more info
#Cols<-c(rgb(.2,0,0.6,1),"white","darkorange")
  
# Divergent palette, not color blind friendly. better with larger numbers of loggers
divergeCols = c( "cornflowerblue","burlywood1", "orangered2")
Func<-colorRampPalette(colors = divergeCols, bias = 0.7)
blueRed<-Func(5)
barplot(rep(1,5),col=blueRed, main="Temperature palette")
# A pretty one if you're not too worried about individual loggers and just want broad patterns
warmColors = c("orangered2", "sienna2", "tan2", "firebrick", "tomato3")


## Set up a theme(s)
shortTC =  theme_classic() + theme(axis.text.x = element_text(angle = 45)) 

longTC = theme_classic() + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
                                 axis.ticks.x = element_blank(), text = element_text(family = "Times"), 
                                 axis.title = element_text(size = 14), legend.title = element_blank(),
                                 legend.position = "bottom", plot.title = element_text(size = 18))

pptTheme = theme_classic() + theme(axis.text = element_text(size = 16), axis.title.x = element_blank(), 
                                   axis.ticks.x = element_blank(), text = element_text(family = "Times"), 
                                   axis.title = element_text(size = 20), 
                                   plot.title = element_text(size = 36))

## You are likely more interested in knowing the chamber position instead of the logger name
  # Create a labelling function to plot the position instead
levels(factor(all_dat$Logger))

LoggerLabels = as_labeller(c(`CCR2` = "St4", `GI` = "St1", `LUX` = "St2", `PIF4` = "Ab2", `ZTL` = "Mu2"))



## Color by datalogger (i.e. 'chamber position')
all_dat %>%
  ggplot(aes(x = xvar, group = Logger)) +
  geom_line(aes(y = AvgTemp, color = Logger)) +
  scale_color_manual(values = blueRed) +
  shortTC

## For just a few loggers, might be easier to just pick individual values
all_dat %>%
  ggplot(aes(x = xvar, group = Logger)) +
  geom_line(aes(y = AvgTemp, color = Logger)) +
  scale_color_manual(values = c("#1F487E", "#B03707", "#2BB6B6", "#450920", "#45305F"), labels = LoggerLabels) +
  theme_classic() +
  ylab("Average air temperature (C)") +
  ggtitle("*Arabidopsis freeze trials*") +
  longTC +
  theme(plot.title = element_markdown()) 
# Might be easier to see the patterns if they're separate
  # + facet_wrap(~Logger, labeller = LoggerLabels)




