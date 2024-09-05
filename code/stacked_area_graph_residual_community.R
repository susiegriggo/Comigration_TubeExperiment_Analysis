#### for species 1 to 30 

library(ggplot2)
library(dplyr)
library(scales)
library(viridis)
library(tidyr)
library(RColorBrewer)

setwd('../')

bacteria <- read.csv("output_files/residual_community_bracken_species_abund_confidence_0.1.tsv", stringsAsFactors = FALSE, sep='\t')
metadata <- read.csv("metadata/residual_community_metadata.csv")

# Step 1: Reorder columns in bacteria based on metadata
bacteriaReorder <- bacteria[, metadata$FAME.ID]

# Subset bacteriaNEW based on metadata$FAME.ID..leave.blank.
selected_columns <- c("name", metadata$FAME.ID)  # Include "name" as the first column

# Subset bacteriaNEW
bacteriaReorder <- bacteriaNEW[, selected_columns]

# Add a new column 'TotalAbundance' that sums the values of the FAME columns for each row
bacteriaReorder <- bacteriaReorder %>%
  rowwise() %>%
  mutate(TotalAbundance = sum(c_across(starts_with("FAME")), na.rm = TRUE)) %>%
  ungroup()

# Name the first column
colnames(bacteriaReorder)[1] <- "Species"

# Generate new column names
new_column_names <- paste( 1:30, sep = "_")

# Rename columns 2 to 30
colnames(bacteriaReorder)[2:30] <- new_column_names

# Arrange the dataset by total_abundance in descending order
bacteriaReorder <- bacteriaReorder %>%
  arrange(desc(TotalAbundance))

# Step 4: Select top 30 rows
top30_bacteriaReorder <- bacteriaReorder[1:30, ]

# Step 5: Reshape data to long format
long_bacteriaReorder <- gather(top30_bacteriaReorder, key = "Sample", value = "Abundance", -Species)

# Ensure 'Time' is numeric and 'Value' is numeric
long_bacteriaReorder$Sample <- as.numeric(as.character(long_bacteriaReorder$Sample))
long_bacteriaReorder$Abundance <- as.numeric(as.character(long_bacteriaReorder$Abundance))

View(long_bacteriaReorder)

# Assuming your dataframe is named df

long_bacteriaReorder <- long_bacteriaReorder %>%
  filter(!is.na(Abundance) & is.finite(Abundance))

# Step 2: Create desired order based on the current sequence
desired_order <- unique(long_bacteriaReorder$Species)  # Extract unique species names in their current order

# Step 3: Convert species to factor with desired order
long_bacteriaReorder$Species <- factor(long_bacteriaReorder$Species, levels = desired_order)

## Create custom color palette using the Spectral palette for 30 groups
spectral_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(37)
print(spectral_palette)

# Create a copy of your data to manipulate
long_bacteriaReorder_modified <- long_bacteriaReorder

# Identify the species to highlight in black
species_to_highlight <- "Enterobacter sp. RHBSTW-00994"


##### reorder colours for new order

custom_palette3 <- c(
  "Bacteroides graminisolvens" = "black",
  "uncultured Bacteroides sp." = "black", 
  "Aeromonas caviae" = "#5E4FA2",
  "Aeromonas hydrophila" = "#515EA9", 
  "Aeromonas media" = "#456EB1", 
  "Morganella morganii" = "#397EB8", 
  "Lactococcus raffinolactis" = "black", 
  "Raoultella ornithinolytica" = "black", 
  "Citrobacter freundii" = "#378EBA",
  "Aeromonas veronii" = "#469EB3", 
  "Kluyvera cryocrescens" = "#54AEAD",
  "Escherichia coli" = "#63BEA6", 
  "Microvirgula aerodenitrificans" = "#75C8A4", 
  "Klebsiella pneumoniae" = "black", 
  "Aeromonas allosaccharophila" = "#9BD7A4",
  "Enterobacter hormaechei" = "#AEDEA3", 
  "Leclercia adecarboxylata" = "#BEE4A0", 
  "Lelliottia amnigena" = "#EAF69E", 
  "Citrobacter braakii" = "#F1F9A9", 
  "Enterobacter roggenkampii" = "#F8FCB4", 
  "Enterobacter asburiae" = "#FEF6B0", 
  "Salmonella enterica" = "#FEEDA2", 
  "Hafnia paralvei" = "#FEE593", 
  "Leclercia sp. 119287" = "#FDCC7A", 
  "Leclercia sp. G3L" = "#FDBE6F", 
  "Clostridium beijerinckii" = "#FDB063", 
  "Enterobacter sp. RHBSTW-00994" = "#FB9F5A", 
  "Citrobacter sp. RHBSTW-00137" = "#F88D52", 
  "Enterobacter sp. 638" = "#F88D52", 
  "Clostridium sp. MF28" = "#F26A43"
)

custom_palette <-  custom_palette3

# Extract the species and their abundances in the first sample from top30_bacteriaReorder
species_order <- top30_bacteriaReorder %>%
  select(Species, `1`) %>%
  arrange(desc(`1`)) %>%
  pull(Species)

# Convert the Species column to a factor with the desired order
top30_bacteriaReorder$Species <- factor(top30_bacteriaReorder$Species, levels = species_order)

# Reshape the data to long format for ggplot
long_top30_bacteriaReorder <- top30_bacteriaReorder %>%
  pivot_longer(cols = -c(Species, TotalAbundance), names_to = "Sample", values_to = "Abundance")

# Convert Sample to numeric for proper ordering in the plot
long_top30_bacteriaReorder$Sample <- as.numeric(long_top30_bacteriaReorder$Sample)

# Plot with legend 
ggplot(long_top30_bacteriaReorder, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_area(position = "fill", colour = "black", alpha = 0.75) +  # Adjust alpha here for stronger transparency
  labs(x = "Distance (cms)", y = "Abundance", fill = "Species") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),  # Bold and larger x-axis title
    axis.title.y = element_text(face = "bold", size = 14),  # Bold and larger y-axis title
    axis.text = element_text(size = 12),  # Larger axis text
    legend.title = element_text(face = "bold", size = 14),  # Bold and larger legend title
    legend.text = element_text(face = "italic", size = 10),  # Italic legend text
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  scale_fill_manual(values = custom_palette3) +  # Custom color palette
  geom_vline(xintercept = 28, linetype = "solid", color = "white", size = 1.5) + # Add vertical white line
  scale_x_continuous(breaks = c(0, 10, 20, 30),  # Specify the positions of the breaks
                     labels = c("0", "30", "60", "90"))  # Specify the labels


###plot without legend 
ggplot(long_top30_bacteriaReorder, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_area(position = "fill", colour = "black", alpha = 0.75) + 
  labs(x = "Distance (cms)", y = "Abundance", fill = "Species") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),  # Bold and larger x-axis title
    axis.title.y = element_text(face = "bold", size = 20),  # Bold and larger y-axis title
    axis.text = element_text(size = 12),  # Larger axis text
    legend.title = element_text(face = "bold", size = 14),  # Bold and larger legend title
    legend.text = element_text(face = "italic", size = 10),  # Italic legend text
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  scale_fill_manual(values = custom_palette3) +  # Custom color palette
  geom_vline(xintercept = 28, linetype = "solid", color = "white", size = 1.5) + # Add vertical white line
  scale_x_continuous(breaks = c(0, 10, 20, 30),  # Specify the positions of the breaks
                     labels = c("0", "30", "60", "90")) +  # Specify the labels
  guides(fill = FALSE)  # Remove the legend


###### 31-61##########################################################################################

# Step 4: Select top 30 rows
second30_bacteriaReorder <- bacteriaReorder[31:61, ]
View(second30_bacteriaReorder)

# Step 5: Reshape data to long format
library(tidyr)
long_bacteriaReorder2nd30 <- gather(second30_bacteriaReorder, key = "Sample", value = "Abundance", -Species)
View(long_bacteriaReorder2nd30)


# Ensure 'Time' is numeric and 'Value' is numeric
long_bacteriaReorder2nd30$Sample <- as.numeric(as.character(long_bacteriaReorder2nd30$Sample))
long_bacteriaReorder2nd30$Abundance <- as.numeric(as.character(long_bacteriaReorder2nd30$Abundance))

library(dplyr)


# Step 2: Create desired order based on the current sequence
desired_order2 <- unique(long_bacteriaReorder2nd30$Species)  # Extract unique species names in their current order

# Step 3: Convert species to factor with desired order
long_bacteriaReorder2nd30$Species <- factor(long_bacteriaReorder2nd30$Species, levels = desired_order2)



custom_palette4 <- c(
  "Aliarcobacter butzleri" = "#9E0142",
  "uncultured Tolumonas sp." = "black",
  "Kluyvera ascorbata" = "#AD1145",
  "Aeromonas rivipollensis" = "#BC2249",
  "Lactococcus lactis" = "#D8434D",
  "Acinetobacter johnsonii" = "black",
  "Vagococcus fluvialis" = "#E1504A",
  "Aeromonas dhakensis" = "#F26A43",
  "Prevotella herbatica" = "black",
  "Citrobacter pasteurii" = "#FDB063",
  "Kluyvera intermedia" = "#FDBE6F",
  "Aeromonas sp. ASNIH1" = "#FDDA86",
  "Enterobacter kobei" = "#FEE593",
  "Citrobacter portucalensis" = "#FEEDA2",
  "Hafnia alvei" = "#FEF6B0",
  "Comamonas terrigena" = "#FFFFBF",
  "Enterobacter ludwigii" = "#F8FCB4",
  "Serratia marcescens" = "#F1F9A9",
  "Citrobacter sp. RHBSTW-00017" = "#EAF69E",
  "Citrobacter werkmanii" = "#DFF299",
  "Serratia liquefaciens" = "#CFEB9C",
  "Citrobacter sp. RHBSTW-00881" = "#BEE4A0",
  "Leclercia sp. 29361" = "#AEDEA3",
  "Leclercia sp. LSNIH3" = "#9BD7A4",
  "Citrobacter sp. RHBSTW-00944" = "#88CFA4",
  "Clostridium diolis" = "#75C8A4",
  "Leclercia sp. Colony189" = "#63BEA6",
  "Leclercia sp. W17" = "#54AEAD",
  "Leclercia sp. LSNIH1" = "#469EB3",
  "Salmonella sp. SSDFZ54" = "#378EBA",
  "Enterobacter sp. MGH 14" = "#397EB8"
)

######### change the order ######

# Extract the species and their abundances in the first sample from top30_bacteriaReorder
species_order <- second30_bacteriaReorder %>%
  select(Species, `1`) %>%
  arrange(desc(`1`)) %>%
  pull(Species)

# Convert the Species column to a factor with the desired order
second30_bacteriaReorder$Species <- factor(second30_bacteriaReorder$Species, levels = species_order)

# Reshape the data to long format for ggplot
long_second30_bacteriaReorder <- second30_bacteriaReorder %>%
  pivot_longer(cols = -c(Species, TotalAbundance), names_to = "Sample", values_to = "Abundance")

# Convert Sample to numeric for proper ordering in the plot
long_second30_bacteriaReorder$Sample <- as.numeric(long_second30_bacteriaReorder$Sample)

# Plot with legend
ggplot(long_second30_bacteriaReorder, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_area(position = "fill", colour = "black", alpha = 0.75) +  # Adjust alpha here for stronger transparency
  labs(x = "Distance (cms)", y = "Abundance", fill = "Species") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),  # Bold and larger x-axis title
    axis.title.y = element_text(face = "bold", size = 14),  # Bold and larger y-axis title
    axis.text = element_text(size = 12),  # Larger axis text
    legend.title = element_text(face = "bold", size = 14),  # Bold and larger legend title
    legend.text = element_text(face = "italic", size = 10),  # Italic legend text
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  scale_fill_manual(values = custom_palette4) +  # Custom color palette
  geom_vline(xintercept = 28, linetype = "solid", color = "white", size = 1.5) + # Add vertical white line
  scale_x_continuous(breaks = c(0, 10, 20, 30),  # Specify the positions of the breaks
                     labels = c("0", "30", "60", "90"))  # Specify the labels

# Plot with no legend 
species_order2 <- levels(long_second30_bacteriaReorder$Species)
print(species_order2)

ggplot(long_second30_bacteriaReorder, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_area(position = "fill", colour = "black", alpha = 0.75) +  # Adjust alpha here for stronger transparency
  labs(x = "Distance (cms)", y = "Abundance", fill = "Species") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),  # Bold and larger x-axis title
    axis.title.y = element_text(face = "bold", size = 14),  # Bold and larger y-axis title
    axis.text = element_text(size = 12),  # Larger axis text
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none"  # Remove the legend
  ) +
  scale_fill_manual(values = custom_palette4) +  # Custom color palette
  geom_vline(xintercept = 28, linetype = "solid", color = "white", size = 1.5) +  # Add vertical white line
  scale_x_continuous(breaks = c(0, 10, 20, 30),  # Specify the positions of the breaks
                     labels = c("0", "30", "60", "90"))  # Specify the labels




