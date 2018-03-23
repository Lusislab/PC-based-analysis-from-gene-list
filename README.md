## PC-based-analysis-from-gene-list
# deposited by Marcus
# R-scripts which pull specific datasets from the Lusislab SQL server, then takes user-defined lists to asses variation of gene sets across the HMDP
# This example focuses on the atherosclerosis study, where a list of human GWAS candidate genes are used to generante the PCs and correlated against traits

## Packages used
RODBC
# Pull data from SQL server

reshape2
# dcast function used to generate aggregate matrix for correlations

WGCNA 
bicorAndPvalue function used to construct correlations

gplots
heatmap.2 function to generate heatmap of gene X gene and gene x trait correlations

factoextra
fviz_pca_var function used to visualize gene-contributions to PCs

## Plots generated
The following plots are generated from the script:

Human GWAS list HMDP gene X gene correlation.pdf
Heatmap showing correlation structure of specified genes accross the HMDP Aorta expression arrays

PC contribution - ALL genes.pdf
Contribution of individual genes to PC1 and PC2

PC contribution - outliers removed.pdf
Contribution of individual genes to PC1 and PC2 after removing the "outliers" shown in plot above


                                              
