###
Use DESeq2 to analyse changes in the abundance of the subsystems categories 
###

# import libraries 
library(EnhancedVolcano)
library(DESeq2)


# set the working directory 
setwd('../')

# read in unnormalised counts from long term migration 
level3 <- read.csv('data/superfocus/long_term_migration_level3_readcounts.tsv', sep = '\t')

# replace colnames 
colnames(level3) <- substr(colnames(level3), 1, 10)

# read in the metadata 
metadata <- read.table(file = 'metadata', row.names = NULL, sep = ',', header = TRUE)
colnames(metadata) <- c("Sample", "id","sample_id" ,"tube","sewage_date", "distance_cm","location","all"  )
#rownames(metadata) <- metadata$id
metadata <- metadata[, c("id","sample_id" ,"tube","sewage_date", "distance_cm","location","all"  )]

# create a deseq object 
dds <- DESeqDataSetFromMatrix(countData=level3, 
                              colData=metadata, 
                              design=~tube + location, tidy = TRUE) # is this correct design to pair location and tube 

# run DESEQ on the object 
dds <- DESeq(dds)

# look at the results table 
res <- results(dds, contrast=c("location","end", "start"))
head(results(dds, tidy=TRUE)) 

# sort the results by p-value
res[order(res$padj),]

# example of plotting a count 
plotCounts(dds, gene="Bacterial Chemotaxis", intgroup="location")


keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$pvalue < 10e-15, '#377eb8',
  ifelse(res$log2FoldChange > 0.5 & res$pvalue < 10e-15, '#e41a1c',
         'black'))
         keyvals[is.na(keyvals)] <- 'black'
           names(keyvals)[keyvals == '#e41a1c'] <- 'Enriched after migration'
           names(keyvals)[keyvals == '#377eb8'] <- 'Enriched before migration'

# make an enhanced volcano plot 
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                pCutoff = 10e-15,
                colCustom = keyvals,
                FCcutoff = 0.5,
                pointSize = 1.6,
                labSize = 0) +  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


### Make a plot where the top 10 significantly changed are shown as a bar plot and the top ten decreased as separate plots 


# attempt to plot some bar charts of the log fold change and difference in proportion 
r <- as.data.frame(res)
sig <- as.data.frame(subset(res, padj<.1e-15 & abs(log2FoldChange)>0.5))
sig$Subsystem <- rownames(sig)
sig<- sig[order(sig$log2FoldChange), ]

# read in mapping of the different level 1 subsystems 
level1 <- read.csv('lvl1_lvl3_subsystems.tsv', sep = '\t')
increase_level1 <- level1[level1$Subsystem.Level.3 %in% rownames(sig_increase), ]
increase_level1 <- increase_level1[order(match(increase_level1$Subsystem.Level.3, rownames(sig_increase))), ]
decrease_level1 <- level1[level1$Subsystem.Level.3 %in% rownames(sig_decrease), ]
decrease_level1 <- decrease_level1[order(match(decrease_level1$Subsystem.Level.3, rownames(sig_decrease))), ]


# sig increased 
sig_increase <- as.data.frame(subset(res, padj<.1e-15 & log2FoldChange>0.5))
sig_increase$Subsystem <- rownames(sig_increase)
sig_increase<- sig_increase[order(sig_increase$log2FoldChange), ]


# Create a basic barplot with reordered y-axis (vertical) and gradient fill
ggplot(head(sig_increase,10), aes(y = reorder(Subsystem, log2FoldChange), x = log2FoldChange, fill = log2FoldChange)) +
  geom_bar(stat = "identity", color = "white", position = "dodge", width = 0.7) +  # Remove black outline, add gap
  labs(
    x = "log2FoldChange",
    y = "Subsystem"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")  # Very light grey color
  ) +
  scale_fill_gradient(low = "orange", high = "red") +  # Adjust the colors for the gradient
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



# sig decreased 
sig_decrease <- as.data.frame(subset(sig, padj<.1e-15 & log2FoldChange<0.5))
sig_decrease$Subsystem <- rownames(sig_decrease)
sig_decrease<- sig_decrease[order(sig_decrease$log2FoldChange),]

# sig increased 
sig_increase <- as.data.frame(subset(sig, padj<.1e-15 & log2FoldChange>0.5))
sig_increase$Subsystem <- rownames(sig_increase)
sig_increase<- sig_increase[order(sig_increase$log2FoldChange),]
rownames(sig_increase) <- c("ABC transporter branched-chain amino acid",         
                            "Alkylphosphonate utilization" ,                                    
                            "Bacterial Chemotaxis",                                       
                            "Pyrimidine utilization",                                           
                            "p-Aminobenzoyl-Glutamate Utilization",                             
                           "Autoinducer 2 transport and processing",
                           "Malonate decarboxylase",
                           "Two partner secretion pathway")

# save all of the significant changes to a dataframe 
write.csv(sig, 'signficiant_changed_subsystems.csv', row.names = TRUE)

# Create a basic barplot with reordered y-axis (vertical) and gradient fill
ggplot(head(sig_decrease,30), aes(y = reorder(Subsystem, log2FoldChange), x = log2FoldChange, fill = log2FoldChange)) +
  geom_bar(stat = "identity", color = "white", position = "dodge", width = 0.7) +  # Remove black outline, add gap
  labs(
    x = "log2FoldChange",
    y = "Subsystem"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")  # Very light grey color
  ) +
  scale_fill_gradient(low = "light blue", high = "dark blue") +  # Adjust the colors for the gradient
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  # Rotate x-axis labels for better readability


# Create a basic barplot with reordered x-axis
# would be good to do this analysis so that I can pair the data and have a standard error of the mean 
# Does Deseq2 have an option for paired data 

ggplot(sig_increase, aes(x = log2FoldChange, y = reorder(Subsystem, log2FoldChange))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    x = "log2FoldChange",
    y = "Subsystem"
  ) +
  theme_bw() +  # Set background to white
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines


# Could have a plot which shows the changes in the gene expression data through time and plot chemotaxis and then have functions within chemotaxis as separate lines 

# save the results table to text file 
write.csv(r, file = 'deSeq_results_level3_subsystems.csv', row.names = TRUE)


