library(vegan)
library(ecodist)
library(ggplot2)
library(dplyr)
library(tidyr) 

setwd('~/OneDrive - Flinders/HONOURS 2021/Combined Paper 2023/')

## Figure on Susie's data  
## PCO ordination at the species level 
data <- read.csv('Susie Data/level3_readcounts.tsv', sep = '\t', row.names = 1)
data[is.na(data)] <- 0
data <- subset(data, rowSums(data != 0) > 0)
meta <- read.table(file = 'Susie Data/metadata.csv',sep = ',', header = TRUE)

# calculate the bray curtis distance - can take some time 
bray_curtis_dist <- vegdist(t(sqrt(data)), method = 'bray')

# compute the pcoa 
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

#get the the variance 
per_variance = 100*(bray_curtis_pcoa$values/sum(bray_curtis_pcoa$values))

# Create a plot
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = meta$location)) +
  geom_point(  aes(shape = meta$sewage_date), size = 3) +
  labs(x = paste( "PCo1 (" ,toString(round(per_variance[1],2)) ,'%)'),
       y = paste( "PCo2 (" ,toString(round(per_variance[2],2)) ,'%)') ) + stat_ellipse(size=0.75, show.legend = F) +
  theme_bw()+
  labs(shape="Sewage sample", colour="Tube position")+
  scale_color_brewer(palette = 'Set1') 
bray_curtis_plot 

## PERMANOVA
## difference between the beginning and the end 
permanova_location <- adonis2(bray_curtis_dist ~ meta$location, permutations=9999)
## difference between the sewage date 
permanova_sewagedate <- adonis2(bray_curtis_dist ~ meta$sewage_date, permutations=9999)

## PERMDISP 
permdisp_location <- betadisper(bray_curtis_dist, meta$location) 
permdisp_results <- permutest(permdisp_location, permutations = 9999) 


### Repeat for Jims functions 

## PCO ordination at the species level 
  
data <- read.csv('Jim Data/superfocus_normed/level3_abundance.tsv', sep = '\t', row.names = 1)
data[is.na(data)] <- 0
data <- subset(data, rowSums(data != 0) > 0)
metadata <- read.table(file='Jim Data/Jim_metadata .csv', fill=TRUE, sep=',', row.names=1, header=TRUE)

# order the samples by their order in the tube 
data <- data[rownames(metadata)]

# calculate the bray curtis distance - can take some time 
bray_curtis_dist <- vegdist(t(sqrt(data)), method = 'bray')

# compute the pcoa 
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

# Calculate lead values for x and y
bray_curtis_pcoa_df$lead_pcoa1 <- c(bray_curtis_pcoa_df$pcoa1[-1], NA)
bray_curtis_pcoa_df$lead_pcoa2 <- c(bray_curtis_pcoa_df$pcoa2[-1], NA)

bray_curtis_pcoa_df[is.na(bray_curtis_pcoa_df)] <- 0

#replace the cells with nan 
bray_curtis_pcoa_df[28,"lead_pcoa1"] <- NA
bray_curtis_pcoa_df[28, "lead_pcoa2"] <- NA

#get the the variance 
per_variance = 100*(bray_curtis_pcoa$values/sum(bray_curtis_pcoa$values))

# Create a plot with arrows between points
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df[1:28,], aes(x = pcoa1, y = pcoa2)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = lead_pcoa1, yend = lead_pcoa2),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),  # Adjust arrow size as needed
               size = 0.5) +
  geom_point(data = bray_curtis_pcoa_df[29,], aes(color = "Band", shape = "Residual Community"), size = 5) +  
  labs(x = paste("PCoA1 (", toString(round(per_variance[1], 2)), '%)'),
       y = paste("PCoA2 (", toString(round(per_variance[2], 2)), '%)')) +
  theme_bw() +
  labs(shape = "Sewage sample", colour = "Tube position")
bray_curtis_plot 
