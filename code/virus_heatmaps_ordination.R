## ordination and heatmaps of the viruses

library(pheatmap)
library(RColorBrewer)
library(vegan)
library(ecodist)
library(ggplot2)
library(dplyr)

setwd('../')

### Make heatmap of long term migration viruses 
# read in the data 
data = read.table(file='output_files/long_term_migration_virus_hecatomb_species_0.05.tsv', fill=TRUE, sep = '\t',header=TRUE, row.names=1)
colnames(data) <- lapply(colnames(data), function(x) substr(x, 1, 10))

# read in the metadata
metadata <- read.table(file='metadata/long_term_migration_metadata.csv', sep = ',', row.names=2, header=TRUE)

# data  
#TODO reorganise the data so that it is in the right order 
colnames(data) <- metadata$location

# Modify ordering of the clustering clustering callback option
callback = function(hc, mat){
  mat[is.na(mat)] <- 0
  sv = svd(mat)$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

#change the names to italics 
newnames <- gsub("\\s+", " ", rownames(data))
newnames <- lapply(
  newnames,
  function(x) bquote(italic(.(x))))


data <- data[, order(names(data), decreasing=TRUE)]
df <- as.matrix(data)
df[df == 0] <- NA


# Check which rows have a maximum value less than 0.05
rows_to_remove <- apply(df, 1, function(row) max(row) < 0.10)

# Remove rows where the maximum value is less than 0.05
df_filtered <- df[!rows_to_remove, ]

# Remove rows with NaN values
df_filtered <- df_filtered[complete.cases(df_filtered), ]

heatmap_plot <- pheatmap(sqrt(df_filtered), colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "PuOr")))(100),cluster_rows = T,  cluster_cols = F,  clustering_callback = callback, fontsize_row =8.5, fontsize_column=8.5,show_colnames = T,  cellheight=12, cellwidth =8,labels_row = as.expression(newnames), gaps_col=c(12, 12) )
heatmap_plot$legend$bottom <- FALSE  # Remove legend from the bottom
heatmap_plot$legend$right <- FALSE  # Remove legend from the right
heatmap_plot$legend$top <- TRUE  # Show legend on the top
heatmap_plot$legend$labels$size <- 8  # Adjust label size
heatmap_plot$legend$labels$padding <- 5  # Adjust padding between labels
heatmap_plot

ggsave("figures/long_term_migration_virus_species_heatmap.png", plot = heatmap_plot, width = 10, height = 3, units = "in", 
       limitsize = FALSE, 
       antialias = "cleartype", ) 


## PCoA on the residual community 
data <- read.table(file = 'output_files/residual_community_virus_hecatomb_species.tsv', header = TRUE, row.names = 1, sep = '\t')
data[is.na(data)] <- 0
data <- subset(data, rowSums(data != 0) > 0)
data <- data[complete.cases(data), ]
data < data.frame(data)

#overwrite the column names to be consistent with the metadata 
colnames(data) <- substr(colnames(data), start = 1, stop = 10)

# read in the metadata
metadata <- read.table(file='metadata/residual_community_metadata.csv', fill=TRUE, sep=',', row.names=1, header=TRUE)

# order the samples by their order in the tube 
data <- data[rownames(metadata)]

# calculate the bray curtis distance - can take some time 
bray_curtis_dist <- vegdist(t(sqrt(data)), method = 'bray')

# compute the pcoa 
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
bray_curtis_pcoa_df$Order <- metadata$Order
#get the the variance 
per_variance = 100*(bray_curtis_pcoa$values/sum(bray_curtis_pcoa$values))

# Assuming you have already executed the code to convert "Order" to a factor with numerical order levels.

# Calculate lead values for x and y
bray_curtis_pcoa_df$lead_pcoa1 <- c(bray_curtis_pcoa_df$pcoa1[-1], NA)
bray_curtis_pcoa_df$lead_pcoa2 <- c(bray_curtis_pcoa_df$pcoa2[-1], NA)

bray_curtis_pcoa_df[is.na(bray_curtis_pcoa_df)] <- 0

#replace the cells with nan 
bray_curtis_pcoa_df[28,"lead_pcoa1"] <- NA
bray_curtis_pcoa_df[28, "lead_pcoa2"] <- NA

# Create a plot with arrows between points
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df[1:28,], aes(x = pcoa1, y = pcoa2)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = lead_pcoa1, yend = lead_pcoa2),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),  # Adjust arrow size as needed
               size = 0.5) +
  geom_point(data = bray_curtis_pcoa_df[29,], aes(color = "Band", shape = "Residual Community"), size = 3) +  
  labs(x = paste("PCoA1 (", toString(round(per_variance[1], 2)), '%)'),
       y = paste("PCoA2 (", toString(round(per_variance[2], 2)), '%)')) +
  theme_bw() +
  labs(shape = "Sewage sample", colour = "Tube position")


# generate ordination showing the correlation 
bray_curtis_plot_with_order <- bray_curtis_plot + geom_segment(aes(x = 0, y = 0, xend = cor(bray_curtis_pcoa_df$pcoa1, bray_curtis_pcoa_df$Order), 
                                                                   yend = cor(bray_curtis_pcoa_df$pcoa2, bray_curtis_pcoa_df$Order)),
                                                               arrow = arrow(type = "closed", length = unit(0.15, "inches")), 
                                                               color = "red", size = 1.5)
 

# Display the plot
bray_curtis_plot

ggsave("figures/residual_community_virus_species_PCoA.png", plot = bray_curtis_plot + theme(aspect.ratio = 1), width = 5, height = 5, units = "in")


### Repeat ordination for long term migration tube samples 
data <- read.table(file = 'output_files/long_term_migration_virus_hecatomb_species.tsv', header = TRUE, row.names = 1, sep = '\t')
data[is.na(data)] <- 0
data <- subset(data, rowSums(data != 0) > 0)
data <- data[complete.cases(data), ]
data < data.frame(data)

#overwrite the column names to be consistent with the metadata 
colnames(data) <- substr(colnames(data), start = 1, stop = 10)

# read in the metadata
metadata <- read.table(file='metadata/long_term_migration_metadata.csv', fill=TRUE, sep=',', row.names=1, header=TRUE)

# order the samples by their order in the tube 
data <- data[metadata$FAME]

# calculate the bray curtis distance - can take some time 
bray_curtis_dist <- vegdist(t(sqrt(data)), method = 'bray')

# compute the pcoa 
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

#get the the variance 
per_variance = 100*(bray_curtis_pcoa$values/sum(bray_curtis_pcoa$values))

# Create a plot
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = metadata$location)) +
  geom_point(  aes(shape = metadata$sewage_date), size = 3) +
  labs(x = paste( "PCo1 (" ,toString(round(per_variance[1],2)) ,'%)'),
       y = paste( "PCo2 (" ,toString(round(per_variance[2],2)) ,'%)') ) + stat_ellipse(size=0.75, show.legend = F) +
  theme_bw()+
  labs(shape="Sewage sample", colour="Tube position")+
  scale_color_brewer(palette = 'Set1') 
bray_curtis_plot 
ggsave("figures/long_term_migration_virus_species_PCoA.png", plot = bray_curtis_plot + theme(aspect.ratio = 1), width = 5, height = 5, units = "in")


# Do a PERMANOVA to see if there is a significant difference in viruses 
## difference between the beginning and the end 
permanova_location <- adonis2(bray_curtis_dist ~ meta$location, permutations=9999)
## difference between the sewage date 
permanova_sewagedate <- adonis2(bray_curtis_dist ~ meta$sewage_date, permutations=9999)

# Do a PERMDISP to see if there is a signficant difference in dispersion 
permdisp_location <- betadisper(bray_curtis_dist, meta$location)  not sure what is up with this 
permdisp_results <- permutest(permdisp_location, permutations = 9999) 