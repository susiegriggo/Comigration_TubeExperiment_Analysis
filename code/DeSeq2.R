###
# Use DESeq2 to analyse changes in the abundance of the SEED level 3 Subsystems 
###

# import libraries 
library(EnhancedVolcano)
library(DESeq2)

# set the working directory 
setwd('../')

# read in unnormalised counts from long term migration 
level3 <- read.csv('data/superfocus/long_term_migration_level3_readcounts.tsv', sep = '\t')

# fix the column names 
colnames(level3) <- substr(colnames(level3), 1, 10)

# read in the metadata 
metadata <- read.table(file = 'metadata/long_term_migration_metadata.csv', row.names = NULL, sep = ',', header = TRUE)
colnames(metadata) <- c("Sample", "id","sample_id" ,"tube","sewage_date", "distance_cm","location" )
rownames(metadata) <- metadata$id
metadata <- metadata[, c("id","sample_id" ,"tube","sewage_date", "distance_cm","location")]

# create a deseq object 
dds <- DESeqDataSetFromMatrix(countData=level3, 
                              colData=metadata, 
                              design=~tube + location, tidy = TRUE) # is this correct design to pair location and tube 

# run DESEQ on the object 
dds <- DESeq(dds)

# look at the results table 
res <- results(dds, contrast=c("location","end", "start"))
head(results(dds, tidy=TRUE)) 


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


# save the significantly changed subsystems 
r <- as.data.frame(res)
sig <- as.data.frame(subset(res, padj<.1e-15 & abs(log2FoldChange)>0.5))
sig$Subsystem <- rownames(sig)
sig<- sig[order(sig$log2FoldChange), ]

# save all of the significant changes to a dataframe 
write.csv(sig, 'output_files/long_term_migration_signficiant_changed_subsystems.csv', row.names = TRUE)

# save the results table to text file 
write.csv(r, file = 'output_files/long_term_migration_deSeq_results_level3_subsystems.csv', row.names = TRUE)


