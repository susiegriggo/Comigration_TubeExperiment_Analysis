# Co-migration of hundreds of species over metres drives selection and promotes non-motile hitchhikers

Code used to generate figures for the paper 'Co-migration of hundreds of species over metres drives selection and promotes non-motile hitchhikers' by Grigson et. al 

## How-to-use
Raw sequencing files can be obtained from the NCBI Sequence Read Archive, BioProject number PRJNA788639. Sample accession numbers range from SRR17326388-SRR17326409 and SRR29503425-SRR29503452 and were processed using scripts from the atavide pipeline https://github.com/linsalrob/atavide_lite/tree/main/slurm. </br> 

Intermediate files and output are available [here](https://figshare.com/s/953b050065ca18fae420![image](https://github.com/user-attachments/assets/e9268431-d12c-4874-8a6e-f17a4c6de1f6)
). Files should be unpacked into this  respository into the `data` directory to allow generation of figures and statistics in the `code` directory. The notebook `curves_rank_abundance.ipynb` should be run first to generate the output files needed to generate ordination plots. Similarly, the code in `filter_hecatomb.ipynb` should be run before `virus_heatmaps_ordination.R` and `DeSeq2.R` should be run before `chemotaxis.ipynb`. 



