library(pheatmap)
library(RColorBrewer)
library(vegan)
library(ecodist)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd('~/OneDrive - Flinders/HONOURS 2021/Combined Paper 2023/')


####  FIGURE ON JIMS TUBE DATA #### 

# read in the data 
data = read.table(file='Jim Data/kraken_prinseq_paired_genus_0.002.tsv', fill=TRUE, sep = '\t',header=TRUE, row.names=1)
colnames(data) <- lapply(colnames(data), function(x) substr(x, 1, 10))

# read in the metadata
metadata <- read.table(file='Jim Data/Jim_metadata .csv', fill=TRUE, sep=',', row.names=1, header=TRUE)

# data  
data <- data[rownames(metadata)]

# Modify ordering of the clustersusing clustering callback option
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
newnames <- gsub("\\s+", "", rownames(data))
newnames <- lapply(
  newnames,
  function(x) bquote(italic(.(x))))

#data <- data[order(data$FAME000132, decreasing=TRUE), ]
df <- as.matrix(data)
colnames(df) <- c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11" ,"12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24","25", "26", "27", "28", "Band")
df[df == 0] <- NA
pheatmap(sqrt(df), colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "Spectral")))(100),cluster_rows = T,  cluster_cols = F,  clustering_callback = callback, fontsize_row =8.5, fontsize_column=8.5,show_colnames = T,  cellheight=10, cellwidth =15,labels_row = as.expression(newnames), gaps_col=c(28, 28, 28, 28))

#### FIGURE ON SUSIES DATA
# read in the data 
data = read.table(file='Susie Data/kraken_prinseq_paired_genus_0.005.tsv', fill=TRUE, sep = '\t',header=TRUE, row.names=1)
colnames(data) <- lapply(colnames(data), function(x) substr(x, 1, 10))

# read in the metadata
metadata <- read.table(file='Susie Data/metadata.csv', sep = ',', row.names=2, header=TRUE)

# data  
colnames(data) <- metadata$location

# Modify ordering of the clustersusing clustering callback option
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
newnames <- gsub("\\s+", "", rownames(data))
newnames <- lapply(
  newnames,
  function(x) bquote(italic(.(x))))


data <- data[, order(names(data), decreasing=TRUE)]
df <- as.matrix(data)
df[df == 0] <- NA
heatmap_plot <- pheatmap(sqrt(df), colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "Spectral")))(100),cluster_rows = T,  cluster_cols = F,  clustering_callback = callback, fontsize_row =8.5, fontsize_column=8.5,show_colnames = T,  cellheight=12, cellwidth =12,labels_row = as.expression(newnames), gaps_col=c(12, 12))

ggsave("~/Downloads/taxa_heatmap.png", plot = heatmap_plot, width = 10, height = 5, units = "in", 
      limitsize = FALSE, 
       antialias = "cleartype", )

### Repeat at the Species level 
data = read.table(file='Susie Data/kraken_prinseq_paired_genus.tsv', fill=TRUE, sep = '\t',header=TRUE, row.names=1)
colnames(data) <- lapply(colnames(data), function(x) substr(x, 1, 10))

# read in the metadata
metadata <- read.table(file='Susie Data/metadata.csv', sep = ',', row.names=2, header=TRUE)

# data  
colnames(data) <- metadata$location


# Modify ordering of the clustering callback option
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
newnames <- gsub("\\s+", "", rownames(data))
newnames <- lapply(
  newnames,
  function(x) bquote(italic(.(x))))


data <- data[, order(names(data), decreasing=TRUE)]
df <- as.matrix(data)
df[df == 0] <- NA


rownames(df) <- gsub("\\s+", "", rownames(df))

# Check which rows have a maximum value less than 0.05
rows_to_remove <- apply(df, 1, function(row) max(row) < 0.01)

# Remove rows where the maximum value is less than 0.05
df_filtered <- df[!rows_to_remove, ]

# Remove rows with NaN values
df_filtered <- df_filtered[complete.cases(df_filtered), ]

pheatmap(sqrt(df_filtered), colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "Spectral")))(100),cluster_rows = T,  cluster_cols = F,  clustering_callback = callback, fontsize_row =8.5, fontsize_column=8.5,show_colnames = T,  cellheight=10, cellwidth =15,labels_row = as.expression(newnames), gaps_col=c(12, 12))


## DO the PCO ordination

## PCO ordination at the genus level 
data <- read.table(file = 'Susie Data/kraken_prinseq_paired_genus.tsv', fill = TRUE, sep = '\t', row.names = 1, header = TRUE)
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
## difference between the beginning and the end 
permdisp_location <- betadisper(bray_curtis_dist, meta$location) 
permdisp_results <- permutest(permdisp_location, permutations = 9999) 

## PCO ordination at the species level - Jims tube 
data <- read.table(file = 'Jim Data/residual_community_bracken_genus_confidence_0.1.tsv', fill = TRUE, sep = '\t', row.names = 1, header = TRUE)
data[is.na(data)] <- 0
data <- subset(data, rowSums(data != 0) > 0)
data <- data[complete.cases(data), ]

#overwrite the column names to be consistent with the metadata 
colnames(data) <- substr(colnames(data), start = 1, stop = 10)

# read in the metadata
metadata <- read.table(file='Jim Data/Jim_metadata .csv', fill=TRUE, sep=',', row.names=1, header=TRUE)

# order the samples by their order in the tube 
data <- data[rownames(metadata)]

# Transpose the data so that features are in columns and samples are in rows
transposed_data <- t(data)

# calculate the bray curtis distance - can take some time 
bray_curtis_dist <- vegdist(t(sqrt(data)), method = 'bray')

# compute the pcoa 
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
bray_curtis_pcoa_df$Order <- metadata$Order

#get the the variance 
per_variance = 100*(bray_curtis_pcoa$values/sum(bray_curtis_pcoa$values))

# Calculate lead values for x and y
bray_curtis_pcoa_df$lead_pcoa1 <- c(bray_curtis_pcoa_df$pcoa1[-1], NA)
bray_curtis_pcoa_df$lead_pcoa2 <- c(bray_curtis_pcoa_df$pcoa2[-1], NA)

# Calculate the correlation of each feature with both PCoA1 and PCoA2
feature_correlations <- data.frame(
  feature = colnames(transposed_data),
  pcoa1_corr = apply(transposed_data, 2, function(feature) cor(feature, bray_curtis_pcoa_df$pcoa1)),
  pcoa2_corr = apply(transposed_data, 2, function(feature) cor(feature, bray_curtis_pcoa_df$pcoa2))
)

# Calculate the combined correlation (Euclidean norm) for each feature
feature_correlations$combined_corr <- sqrt(feature_correlations$pcoa1_corr^2 + feature_correlations$pcoa2_corr^2)

# Select the top 5 features based on the combined correlation
top_features <- feature_correlations[order(-feature_correlations$combined_corr), ][1:5, ]

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

  
# Add arrows to represent the top 5 features driving the PCoA1 and PCoA2 axes
bray_curtis_plot_with_taxa_correlations <- bray_curtis_plot +
  geom_segment(data = top_features, aes(x = 0, y = 0, 
                                        xend = pcoa1_corr, 
                                        yend = pcoa2_corr),
               arrow = arrow(type = "closed", length = unit(0.15, "inches")),
               color = "blue", size = 1) +
  geom_text(data = top_features, aes(x = pcoa1_corr * 1.1, 
                                     y = pcoa2_corr * 1.1, 
                                     label = feature),
            color = "blue", size = 4, hjust = 0)


# Adding the arrow for 'Order'
bray_curtis_plot_with_order <- bray_curtis_plot + geom_segment(aes(x = 0, y = 0, xend = cor(bray_curtis_pcoa_df$pcoa1, bray_curtis_pcoa_df$Order), 
                                                                     yend = cor(bray_curtis_pcoa_df$pcoa2, bray_curtis_pcoa_df$Order)),
                                                                 arrow = arrow(type = "closed", length = unit(0.15, "inches")), 
                                                                 color = "red", size = 1.5)


# Display the plot
bray_curtis_plot
ggsave("~/Downloads/residual_community_taxa.png", plot = bray_curtis_plot + theme(aspect.ratio = 1), width = 5, height = 5, units = "in")


## PCO ordination at the species level - Susie Data 
data <- read.table(file = 'Susie Data/kraken_prinseq_paired_species.tsv', fill = TRUE, sep = '\t', row.names = 1, header = TRUE)
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

# maybe pick some different colours for this figure as not to be confused with the colors in the heatmaps 

# Create a plot
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2, color = meta$location)) +
  geom_point(  aes(shape = meta$sewage_date), size = 3) +
  labs(x = paste( "PCo1 (" ,toString(round(per_variance[1],2)) ,'%)'),
       y = paste( "PCo2 (" ,toString(round(per_variance[2],2)) ,'%)') ) + stat_ellipse(size=0.75, show.legend = F) +
  theme_bw()+
  labs(shape="Sewage sample", colour="Tube position")+
  scale_color_brewer(palette = 'Set1') 

bray_curtis_plot




### plot lineplot comparing motile and non-motile taxa 
taxaplot_df <- read.table('Jim Data/motility_kraken_prinseq_paired_genus_0.002.txt', sep = '\t', header = TRUE, row.names = 1)
# subset the non-motile species 
motile_df <- taxaplot_df[taxaplot_df$motile == 'N', ]  
motile_df <- select(motile_df, -c('motile'))
motile_df$Genera <- rownames(motile_df)

# reshape the data so that it is in 'long' format 
plot_motile_df <- gather(motile_df, Sample, Abundance, -Genera)
plot_motile_df[is.na(plot_motile_df)] <- 0

# build this plot 
ggplot(plot_motile_df, aes(x = factor(Sample, levels = colnames(motile_df)), y = Abundance, color = Genera, group = Genera)) +
  geom_point() +
  geom_line() +  # Add a line connecting points
  labs(
    title = "Line Plot by Timepoint",
    x = "Sample in Tube",
    y = "Relative Abundance"
  ) +
  theme_minimal() +  # Change the overall theme to minimal
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_line(color = "white", size = 0.2),  # White major grid lines
    panel.grid.minor = element_line(color = "white", size = 0.1),  # White minor grid lines
    panel.background = element_rect(fill = "white")  # Set the background to white
  ) +
  scale_color_hue() +
  scale_x_discrete(
    limits = colnames(motile_df),
    labels = seq_along(colnames(motile_df))  # Set sequential numbers as labels
  ) +
  scale_y_continuous(limits = c(0, max(plot_motile_df$Abundance)) )  # Set y-axis to start at 0


### Next run correlations between the taxa to try and find motile species which correlate with the motile species 

# subset the non-motile species 
non_motile_df = taxaplot_df[taxaplot_df$motile == 'Y', ] 
non_motile_df <- select(non_motile_df, -c('motile'))
non_motile_df$Genera <- rownames(non_motile_df)

# reshape the data so that it is in 'long' format 
plot_non_motile_df <- gather(non_motile_df, Sample, Abundance, -Genera)
plot_non_motile_df[is.na(plot_non_motile_df)] <- 0
