# PC-based-analysis-from-gene-list

deposited by Marcus

R-scripts which pull specific datasets from the Lusislab SQL server, then takes user-defined lists to assess variation of gene sets across the HMDP


This example focuses on the atherosclerosis study, where a list of human GWAS candidate genes are used to generate the PCs and correlated against traits

# Packages used
### RODBC

Pull data from SQL server

### reshape2

dcast function used to generate aggregate matrix for correlations

### WGCNA 

bicorAndPvalue function used to construct correlations

### gplots

heatmap.2 function to generate heatmap of gene X gene and gene x trait correlations

### factoextra

fviz_pca_var function used to visualize gene-contributions to PCs

# Plots generated

The following plots are generated from the script:

### Human GWAS list HMDP gene X gene correlation.pdf

Heatmap showing correlation structure of specified genes across the HMDP Aorta expression arrays

### PC contribution - ALL genes.pdf

Contribution of individual genes to PC1 and PC2

### PC contribution - outliers removed.pdf

Contribution of individual genes to PC1 and PC2 after removing the "outliers" shown in plot above

### PC1 and PC2 x trait heatmap.pdf

Heatmap showing correlation of all traits across PC1 and PC2 

# Files

The following are files used to run the scripts

### Human_gwas_gene_list.txt 

lists of genes to be used in constructing Principle component based vectors

### Human_MS_Orths.txt 

Mouse and Human Orthologous genes from JAX

### pre-processed_transcripts.txt pre-processed_traits.txt

Tables of pre-processed files produced from merged SQL data which will be used for downstream analysis

The following files are generated form the scripts

### PC X trait bicor and pvalue.txt

The resulting bicor coefficient and pvalues for PC1 and PC2 across HMDP Ath traits
                                              
