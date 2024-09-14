# Co-migration of hundreds of species over metres drives selection and promotes non-motile hitchhikers

Code used to generate figures for the paper 'Co-migration of hundreds of species over metres drives selection and promotes non-motile hitchhikers' by Grigson et. al 

## How-to-use
Raw sequencing files can be obtained from the NCBI Sequence Read Archive, BioProject number PRJNA788639. Sample accession numbers range from SRR17326388-SRR17326409 and SRR29503425-SRR29503452 and were processed using scripts from the atavide pipeline https://github.com/linsalrob/atavide_lite/tree/main/slurm. </br> 

Intermediate files and output are available [here](https://figshare.com/s/953b050065ca18fae420![image](https://github.com/user-attachments/assets/e9268431-d12c-4874-8a6e-f17a4c6de1f6)
). Files should be unpacked into this  respository into the `data` directory to allow generation of figures and statistics in the `code` directory. The notebook `curves_rank_abundance.ipynb` should be run first to generate the output files needed to generate ordination plots. Similarly, the code in `filter_hecatomb.ipynb` should be run before `virus_heatmaps_ordination.R` and `DeSeq2.R` should be run before `chemotaxis.ipynb`. 

## Dependencies 

**Bioinformatics tools used for processing and analysing metagenomic data included:**
* fastp (version 0.23.4)
* PRINSEQ (1.2.4)
* Kraken2 (version 2.1.3)
* Bracken (version 2.9)
* seqtk (version 1.4)
* SUPERFOCUS (version 1.6)
* Hecatomb (version 1.1.0)
* MEGAHIT (version 1.2.9)
* Koverage (version 0.1.11)
* VAMB (version 4.1.4)
* Checkm2 (version 1.0.1)
* GTDB-Tk (version 2.1.1)
* pyani (version 0.2.12)
  
Where applicable, scripts were utilised from the atavide pipeline https://github.com/linsalrob/atavide_lite. 

**For analysis carried out in R, version 4.4.0 was used in RStudio (version 2024.04.2+764). The following packages were used:**
* ggplot2 (version 3.5.1)
* pheatmap (version 1.0.12)
* ecodist (2.1.3)
* vegan (2.6-8)
* EnhancedVolcano (version 1.22.0)
* DESeq2 (version 1.44.0)
* rstatix (version 0.7.2)

**For analysis carried out in python, version 3.12.0 was used with the packages:**
* seaborn (version 0.13.0)
* statsmodel (0.14.1)
* matplotlib (3.8.2)


