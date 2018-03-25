library(odbc)
library(WGCNA)
library(reshape2)
library(tidyverse)
library(magrittr)
library(factoextra)

# connect to the SQL server
con <- dbConnect(odbc::odbc(),
                 .connection_string = 'SERVER=10.14.48.13;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={ODBC Driver 13 for SQL Server}')

# Retrieve female trait data from HMDP Atherosclerosis panel (there is a view that contains this):
trait = dbGetQuery(con, "select * from ClinicalTraits.AtherosclerosisFemale") %>%
  select(-Maternal_strain, -Asp_Time)

# Retrieve female aorta expression data from HMDP Atherosclerosis panel:
expr_q = "select mice.Maternal_strain, expr.*, info.gene_symbol
from TranscriptAbundance.AtherosclerosisAortaFemale expr
	inner join IndividualMetadata.Atherosclerosis mice
		on expr.mouse_number = mice.mouse_number
	inner join TranscriptAnnotation.[Affymetrix_HT_MG-430_PM_v33] info
		on expr.probesetID = info.probesetID"

expr = dbGetQuery(con, expr_q)

# The next goal will be to use a gene list in order to generate Principal Components from transcripts
# We will use a list of candidate human GWAS genes supplied by Hooman Allayee
# To find mouse orthologues for human genes, we use the deposited list "Human_MS_Orths.txt"
human_mouse_orthologs <- read.delim("Human_MS_Orths.txt", stringsAsFactors = F)

# Now we read in the list of genes used to construct the PC - This is supplied as "Human_gwas_gene_list.txt"
human_gwas_genes <- read.delim("Human_gwas_gene_list.txt", stringsAsFactors = F)

# Add column of orthologous mouse genes, keeping only genes with mouse orthologs
human_gwas_genes %<>%
  inner_join(select(human_mouse_orthologs, human_orth, Symbol), by = c('Gene.based.Analysis1' = 'human_orth')) %>%
  rename(mouse_orth = Symbol)

# get expression for genes with human orthologs, average by gene symbol if more than 1 probeset exists for it,
# and rearrange to have 1 column per gene symbol
expr_human_orthologous = expr %>%
  select(mouse_number, gene_symbol, expression_value) %>%
  semi_join(human_gwas_genes, by = c('gene_symbol' = 'mouse_orth')) %>%
  group_by(mouse_number, gene_symbol) %>%
  summarise(expr_avg = mean(expression_value)) %>%
  spread(gene_symbol, expr_avg) %>%
  ungroup()

# identify mice with both expression and trait data
common = intersect(expr_human_orthologous$mouse_number, trait$mouse_number)

# define filtering function
keep_common = . %>%
  filter(mouse_number %in% common) %>%
  arrange(mouse_number) %>%
  select(-mouse_number)

# keep only the common mice
expr_human_orthologous %<>% keep_common
trait %<>% keep_common

# calculate correlations using biweight midcorrelation (robust to outliers)
gg.cor = bicorAndPvalue(expr_human_orthologous, use = 'pairwise.complete.obs')

# cluster results to make better-looking heatmap
dist_bc = dist(gg.cor$bicor)
hc = hclust(dist_bc)
bicor_values = melt(gg.cor$bicor) %>%
  mutate_at(vars(Var1, Var2), ~fct_relevel(., hc$labels[hc$order])) %>%
  setNames(c('gene1', 'gene2', 'bicor'))

# make heatmap
pdf('Human GWAS list HMDP gene X gene correlation.pdf', width=8, height=7)
ggplot(bicor_values, aes(x=gene1, y=gene2, fill=bicor)) + geom_tile() +
  scale_fill_gradient2(low='yellow', mid='grey', high='blue') +
  xlab('') + ylab('') +
  theme(axis.text = element_text(size=4),
        axis.text.x = element_text(angle=90))
dev.off()

# do principal components analysis (PCA)
pca = prcomp(expr_human_orthologous)
summary(pca)

# visualize PCA results
pdf('PC contribution - ALL genes.pdf', width=11, height=8.5)
fviz_pca_var(pca, col.var = "contrib",
             gradient.cols = c("grey", "blue"),
             ggtheme = theme_minimal())
dev.off()

# remove outlying genes and redo PCA
pca = select(expr_human_orthologous, -Myh11, -Mybphl) %>% prcomp
pdf('PC contribution - outliers removed.pdf', width=11, height=8.5)
fviz_pca_var(pca, col.var = "contrib",
             gradient.cols = c("grey", "blue"),
             ggtheme = theme_minimal())
dev.off()

# use the first 2 PCs - note, can use more by speciying length of columns
comps = pca$x[,1:2]

# get correlation between position on PC axis vs trait
trait.pc.cor = bicorAndPvalue(comps, trait, use = "p")

# perform clustering to get better-looking heatmap
dist_bc = dist(t(trait.pc.cor$bicor))
hc = hclust(dist_bc)

bicor_pvalues = cbind(melt(trait.pc.cor$bicor), melt(trait.pc.cor$p)[,3]) %>%
  mutate_at('Var2', ~fct_relevel(., hc$labels[hc$order])) %>%
  setNames(c('PC', 'trait_name', 'bicor', 'pvalue'))

# export correlation values
bicor_pvalues %>% gather(variable, value, -PC, -trait_name) %>%
  mutate(dest = str_glue_data(., '{PC}_{variable}')) %>%
  select(-PC, -variable) %>%
  spread(dest, value) %>%
  write.table(file="PC X trait bicor and pvalue.txt", row.names=F, col.names=T, sep='\t', quote=F)

# The heatmap produced is listed as "PC1 and PC2 x trait heatmap.pdf"
pdf("PC1 and PC2 x trait heatmap.pdf")
ggplot(bicor_pvalues, aes(x=PC, y=trait_name, fill=bicor)) + geom_tile() +
  scale_fill_gradient2(low='red', mid='white', high='green') +
  xlab('') + ylab('') + theme_bw() +
  theme(axis.text = element_text(size=4))
dev.off()
