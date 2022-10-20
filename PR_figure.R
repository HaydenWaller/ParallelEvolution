library(ggplot2)
library(dplyr)
library(hrbrthemes)

songs<-read.csv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/4_species_RNAseq/4_species_songs.csv", header=TRUE)
ggplot(songs, aes(x=PR_values, fill=Species, colour=Species))+
  geom_histogram(alpha = 0.5, bins=30)+
  scale_x_continuous(breaks=seq(0,4.5, by=0.1))
  theme_minimal()