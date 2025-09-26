library(vegan)
library(ecodist)
library(ggplot2)

setwd('../')


# read in the data 
data <- read.table(file = 'data/16S/community_formation_genus_abund.csv', header=TRUE, sep=',', row.names = 1)

#overwrite the column names to be consistent with the metadata 
metadata <- read.table(file='metadata/community_formation_metadata.csv', fill=TRUE, sep=',', row.names=1, header=TRUE)
metadata$Sewage_sample <- as.character(metadata$Sewage_sample)
metadata$Sewage_label <- paste0('Sewage sample ', metadata$Sewage_sample)
metadata$Time <- as.character(metadata$Time)


# calculate the bray curtis distance - can take some time 
bray_curtis_dist <- vegdist(t(sqrt(data)), method = 'bray')

# compute the pcoa 
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

#get the the variance 
per_variance = 100*(bray_curtis_pcoa$values/sum(bray_curtis_pcoa$values))


# Visualise the pcoa 
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x = pcoa1, y = pcoa2, color = metadata$Time)) + 
  geom_point(aes(shape = metadata$Sewage_label), size = 3) +
  theme_bw() + 
  labs(x = paste("PCo1 (", toString(round(per_variance[1], 2)), "%)"),
       y = paste("PCo2 (", toString(round(per_variance[2], 2)), "%)"),
       color = "Time (hours)",   # Change color legend title
       shape = "Sewage replicate" # Change shape legend title
  ) + 
  scale_color_brewer(palette = 'Set1') 

ggsave("figures/community_formation_genus_PCoA.png", plot = bray_curtis_plot + theme(aspect.ratio = 1), width = 5, height = 5, units = "in")

# do permanova and permdisp on this data 
## 1. Perform PERMANOVA
permanova_result_sewage_date <- adonis2(bray_curtis_dist ~ Sewage_sample, data = metadata, method = "bray", permutations=9999)
print(permanova_result_sewage_date)
permanova_result_location <- adonis2(bray_curtis_dist ~ Time, data = metadata, method = "bray", permutations=9999)
print(permanova_result_location)

## 2. Perform PERMDISP
# Test for homogeneity of dispersions for 'location'
disper_location <- betadisper(bray_curtis_dist, metadata$Time)
permdisp_location <- permutest(disper_location, permutations = 9999)
print(permdisp_location)

