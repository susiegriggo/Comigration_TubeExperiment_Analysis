# libraries 
library(pheatmap)
library(RColorBrewer)
library(vegan)
library(ecodist)
library(ggplot2)
library(dplyr)
library(tidyr) 

setwd('../')

# Residual community PCoA 

# read in the data 
data <- read.table(file = 'output_files/residual_community_bracken_genus_abund_confidence_0.1.tsv', fill = TRUE, sep = '\t', row.names = 1, header = TRUE)
data[is.na(data)] <- 0
data <- subset(data, rowSums(data != 0) > 0)
data <- data[complete.cases(data), ]

#overwrite the column names to be consistent with the metadata 
colnames(data) <- substr(colnames(data), start = 1, stop = 10)
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

# Display the plot
bray_curtis_plot
ggsave("figures/residual_community_genus_PCoA.png", plot = bray_curtis_plot + theme(aspect.ratio = 1), width = 5, height = 5, units = "in")


# correlation of order with the first principal component 
# Ensure both 'Order' and 'PCoA1' are numeric
bray_curtis_pcoa_df$Order <- as.numeric(bray_curtis_pcoa_df$Order)
bray_curtis_pcoa_df$pcoa1 <- as.numeric(bray_curtis_pcoa_df$pcoa1)

# Calculate the Pearson correlation coefficient between 'Order' and 'PCoA1'
correlation_order_pcoa1 <- cor(bray_curtis_pcoa_df$Order, bray_curtis_pcoa_df$pcoa1, method = "pearson")

# Print the correlation coefficient
print(correlation_order_pcoa1)

## PcoA of long-term migration 
# start vs end PCoA 
data <- read.table(file = 'output_files/long_term_migration_bracken_genus_abund_confidence_0.1.tsv', fill = TRUE, sep = '\t', row.names = 1, header = TRUE)
data[is.na(data)] <- 0
data <- subset(data, rowSums(data != 0) > 0)
meta <- read.table(file = 'metadata/long_term_migration_metadata.csv',sep = ',', header = TRUE)

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

# start vs end heatmap 
  # read in the data 
  data = read.table(file='output_files/long_term_migration_bracken_species_abund_confidence_0.1.tsv', fill=TRUE, sep = '\t',header=TRUE, row.names=1)
  colnames(data) <- lapply(colnames(data), function(x) substr(x, 1, 10))
  
  # read in the metadata
  metadata <- read.table(file='metadata/long_term_migration_metadata.csv', sep = ',', row.names=2, header=TRUE)
  
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
  
  data <- data[, order(names(data), decreasing=TRUE)]
  df <- as.matrix(data)
  df[df == 0] <- NA
  
  # filter to include with max above cutoff 
  # Check which rows have a maximum value less than 0.05
  rows_to_remove <- apply(df, 1, function(row) max(row) < 0.05)
  
  # Remove rows where the maximum value is less than 0.05
  df_filtered <- df[!rows_to_remove, ]
  
  # Remove rows with NaN values
  df_filtered <- df_filtered[complete.cases(df_filtered), ]
  
  
  #change the names to italics 
  newnames <- gsub("\\s+", " ", rownames(df_filtered))
  newnames <- lapply(
    newnames,
    function(x) bquote(italic(.(x))))
  heatmap_plot <- pheatmap(sqrt(df_filtered), colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                       "Spectral")))(100),cluster_rows = T,  cluster_cols = F,  clustering_callback = callback, fontsize_row =8.5, fontsize_column=8.5,show_colnames = T,  cellheight=12, cellwidth =8,labels_row = as.expression(newnames), gaps_col=c(12, 12))
  
  ggsave("figures/long_term_migration_species_heatmap.png", plot = heatmap_plot, width = 10, height = 5, units = "in", 
         limitsize = FALSE, 
         antialias = "cleartype", )

# Jim Tube figures 

# read in the data 
data = read.table(file='output_files/residual_community_bracken_species_abund_confidence_0.1.tsv', fill=TRUE, sep = '\t',header=TRUE, row.names=1)
colnames(data) <- lapply(colnames(data), function(x) substr(x, 1, 10))

# read in the metadata
metadata <- read.table(file='metadata/residual_community_metadata.csv', fill=TRUE, sep=',', row.names=1, header=TRUE)

# data  
data <- data[rownames(metadata)]

# Modify ordering of the clusters using clustering callback option
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

#data <- data[order(data$FAME000132, decreasing=TRUE), ]
df <- as.matrix(data)
df[is.na(df)] <- 0
colnames(df) <- c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11" ,"12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24","25", "26", "27", "28", "Band")


# Check which rows have a maximum value less than a cutoff 
rows_to_remove <- apply(df, 1, function(row) max(row) < 0.02)

# Remove rows where the maximum value is less than 0.05
df_filtered <- df[!rows_to_remove, ]


#change the names to italics 
newnames <- gsub("\\s+", "", rownames(df_filtered))
newnames <- lapply(
  newnames,
  function(x) bquote(italic(.(x))))


plot <- pheatmap(sqrt(df_filtered), colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "Spectral")))(100),cluster_rows = T,  cluster_cols = F,  clustering_callback = callback, fontsize_row =8.5, fontsize_column=8.5,show_colnames = T,  cellheight=10, cellwidth =15,labels_row = as.expression(newnames), gaps_col=c(28, 28, 28, 28))

plot 