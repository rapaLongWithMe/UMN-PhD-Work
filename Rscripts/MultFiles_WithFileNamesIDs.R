## Compile multiple files and name with filename
setwd("~/Desktop/Git/RNAseq/Bnapus/Enriched_GO_Terms/20220307/ByGeno")

# Make a function to process each file
processFile <- function(f) {
  df <- read.csv(paste0(f))
}

# Find all .csv files within different folders
  #need recursive = TRUE to look within folders, don't need if wd in a single folder with lists of (ex) csv's
files <- dir(".", recursive=TRUE, full.names=TRUE, pattern="\\.csv$")
files

# Apply the function to all files.
enrich_list <- lapply(files, processFile)

#chop up "files" object to create a list of names as IDs
mods = gsub(".*/", "", files)
mods = gsub(".csv", "", mods)
mods
genos = gsub("_.*", "", files)
genos = gsub("./", "", genos)
genos
#name the elements of the list for clarity
names(enrich_list) = paste(genos, mods, sep = "_")

#Add an id column to each element of the list
enrich_master = list()
for(i in 1:length(enrich_list)) {
  enrich_master[[i]] = enrich_list[[i]] %>%
    mutate(ID = names(enrich_list[i]))
}  



## Bring in the files and combine into a single df
  # set pattern = something common in the filename for files of interest (ex: .txt)
filenames = list.files(pattern = "L02")
filenames
#set up a vector of column names if you only want to keep a subset of the columns. Have to name them all first though
headers = c("Horizontal", "Vertical", "File")

traces = lapply(filenames, function(x) read_csv(x) %>% 
                  #you'll end up with a tibble per csv file
                  bind_rows %>%
                  #add the file name as an identifier
                  mutate(File = paste(x)) %>%
                  data.frame() %>%
                  #set the column names so can pull the ones you want
                  setNames(nm = headers) %>%
                  select(Vertical, File))

#just want the vertical movement data (col = 2) and the filename
traces = do.call("rbind", lapply(traces, function(x) as.data.frame(x)))
#chop up the file name into more useful information ('vector', 'start', 'stop')
traces$Line = substr(traces$File, 1, 3)
traces$Camera = substr(traces$File, 4, 6)
traces$Row = substr(traces$File, 7, 8)
traces$Position = substr(traces$File, 9, 11)
