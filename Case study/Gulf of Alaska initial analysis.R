library(dplyr)

##### Gulf of Alaska initial analysis #####

# Read in data
alaska1 <- read.csv("~/GitHub/extended-independent-project/NOAA Stock Smart data/Assessment_TimeSeries_Data_Part_1.csv")
alaska2 <- read.csv("~/GitHub/extended-independent-project/NOAA Stock Smart data/Assessment_TimeSeries_Data_Part_2.csv")
alaska <- cbind(alaska1, alaska2) # Combine separate files into 1

# Tidying up
alaska <- as.data.frame(t(as.matrix(alaska)))
colnames(alaska)[1:7] <- alaska[1,1:7]
colnames(alaska)[8:73] <- alaska[2,8:73]
alaska <- alaska[-c(1,2),]
