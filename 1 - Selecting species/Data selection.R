##### Choosing a dataset from the RAM Legacy database #####

### Packages
library(lme4)
library(brms)
library(lmerTest)
library(lmtest)

### Load in raw data
load("~/GitHub/Extended-Independent-Project/RAM Legacy R Data/DBdata[asmt][v4.66].RData")

# Reduce stock table to variables I care about
cs <- as.data.frame(cbind(stock$stockid, stock$scientificname, stock$commonname, stock$areaid, stock$region))
colnames(cs) <- c('stockid', 'scientificname', 'commonname', 'areaid', 'region')

# Reduce areas to only ones that are repeated more than twice (have more than 2 species in)
count <- table(cs$areaid)
repeated <- names(count[count > 2])
cs <- cs[cs$areaid %in% repeated,]

# Add area names
areas <- as.data.frame(cbind(area$areaid, area$areaname))
colnames(areas) <- c('areaid', 'areaname')
cs$areaname <- areas$areaname[match(cs$areaid, areas$areaid)]

# Reorder and remove area id
cs <- cs[,c('region','areaname', 'commonname', 'scientificname', 'stockid')]

# How many different areas are there?
n_distinct(cs$areaname) # 132 different areas
