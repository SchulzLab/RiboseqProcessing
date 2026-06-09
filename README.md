# Ribo-seq Preprocessing scripts for Split-ORF analysis

In order to run the ribo-cov module of the [Split-ORF pipeline](https://github.com/SchulzLab/SplitOrfs) the Ribo-seq datasets need to be supplied in preprocessed format.

The ribo-cov module performs the mapping of the adapter trimmed and filtered reads using STAR, if UMIs are involved the data need to be supplied as BAM files that are mapped with STAR and deduplicated.

The scripts in this directory are tailored for the preprocessing of the Ribo-seq data used in the study and were not created or tested for use on other data or systems.

The New_datasets_riboseq_so_paper folder contains the bash scripts used to preprocess glioblastoma (PRJNA591767) and breast cancer (PRJNA523167) datasets.


The UPF1_deletion folder contains the bash script to align the data of Boehm et al. upon UPF1 degradation (E-MTAB- 13837).


All data are preprocessed and or aligned using the Riboseq_analysis_pipeline/Riboseq_pipeline.sh script, which contains different options to enable the preprocessing of different types of Ribo-seq data.
