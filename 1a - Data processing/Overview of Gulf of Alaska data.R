##### Overview of Gulf of Alaska data #####

### Packages
library(ggplot2)
library(ggokabeito)
library(ggpubr)

### Load in data
load("~/GitHub/Extended-Independent-Project/Gulf_Of_Alaska.RData")

### Extract SSB and recruitment data
ssb <- GoA.data[,1:3]
recruits <- GoA.data[,c(1,2,4)]

# SSB plot

p1 <- ggplot(ssb, aes(x = year, y = log(ssb), col = species, group = species)) + 
      geom_line(linewidth = 1.5) + 
      scale_color_okabe_ito(labels = GoA.species$commonname, name = "Species") + 
      scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12)) +
      theme(axis.title = element_text(size = 12)) +
      xlab("Year") +
      ylab("log(SSB [tonnes] )") +
      ggtitle("Spawning stock biomass of 6 species in the Gulf of Alaska 1980-2017")

# Recruitment plot

p2 <- ggplot(recruits, aes(x = year, y = log(recruits), col = species, group = species)) + 
      geom_line(linewidth = 1.5) + 
      scale_color_okabe_ito(labels = GoA.species$commonname, name = "Species") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12)) +
      theme(axis.title = element_text(size = 12)) +
      xlab("Year") +
      ylab("log(R)") +
      ggtitle("Recruitment of 6 species in the Gulf of Alaska 1980-2017")

ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = T, legend = "bottom")

