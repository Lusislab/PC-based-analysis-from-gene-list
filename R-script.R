#set the directory where the files will be written
setwd("C:/Users/mseldin/Desktop/lab files/Hooman Lab/2-12-18 gene-only/outliers removed")

#load packages
library(RODBC)
library(reshape2)
library(WGCNA)
library(gplots)
library(factoextra)

#establish connection with SQL server
ch = odbcDriverConnect('SERVER=10.14.48.13;DATABASE=HMDP;DRIVER={ODBC Driver 13 for SQL Server};Trusted_Connection=Yes')

#pull a specific query from the SQL server - here we are pulling all traits from the HMDP Tables
query1 = "select *
from ClinicalTraits.Atherosclerosis"
traits = sqlQuery(ch, query1, as.is=T, stringsAsFactors=F)

#pull the annotation table for clinical traits
query2 = "select *
from ClinicalTraitAnnotation.Atherosclerosis"
traits.annot = sqlQuery(ch, query2, as.is=T, stringsAsFactors=F)

#pull the metadata for the study
query3 = "select *
from IndividualMetadata.Atherosclerosis"
match = sqlQuery(ch, query3, as.is=T, stringsAsFactors=F)

#add the annotation for each trait, based on a matching "trait_id"
traits$trait_name = traits.annot$trait_name[match(traits$trait_id, traits.annot$trait_id)]
traits$c = traits.annot$trait_category[match(traits$trait_id, traits.annot$trait_id)]
traits = traits[!grepl("lipidomics", traits$c,]
                       
#Add the "Maternal_strain" to each individual based on matching "mouse_number"
traits$strain = match$Maternal_strain[match(traits$mouse_number, match$mouse_number)]

#Make sure the values for each "trait" and "mouse_number" are numeric class
traits$value = as.numeric(traits$value)
traits$mouse_number = as.numeric(traits$mouse_number)

#Here we will only use the female mice - so we remove the males based on the "sex" column
traits$sex = match$Gender[match(traits$mouse_number, match$mouse_number)]
traits =traits[!grepl("M",traits$sex),]

#apply the dcast function to generate a new matrix of traits, where the rows are "mouse_number" and the columns are "trait_name"
trait = dcast(traits, traits$mouse_number ~ traits$trait_name, fun.aggregate = mean, value.var = 'value')

#Now we apply a similar set of functions to the expression data - the goal is to generate a similar matrix of "mouse_number" (rows) and "gene_symbol" (columns)
#Pull the transcript abundance table
query4 = "select *
from TranscriptAbundance.AtherosclerosisAortaFemale"
trx = sqlQuery(ch, query4, as.is=T, stringsAsFactors=F)

#Pull the annotation file containing "gene_symbol" - make sure the arrray matches the Transcript.Abundace table
query5 = "select *
from TranscriptAnnotation.[Affymetrix_HT_MG-430_PM_v33]"
trx.annot = sqlQuery(ch, query5, as.is=T, stringsAsFactors=F)

#Add the gene_symbol based on matching by "probesetID"
trx$gene_symbol = trx.annot$gene_symbol[match(trx$probesetID, trx.annot$probesetID)]

#Make sure "expression_value" and "mouse_number" are numeric class
trx$expression_value = as.numeric(trx$expression_value)
trx$mouse_number = as.numeric(trx$mouse_number)

#use dcast function to generate a table of "mouse_number" as rows and "gene_symbol" as columns
tr = dcast(trx, trx$mouse_number ~  trx$gene_symbol, fun.aggregate = mean, value.var = 'expression_value')

#use the match fucntion to only retain "mouse_number" which match between matrices
tr$m = match(tr$`trx$mouse_number`, trait$`traits$mouse_number`, nomatch = 0)
tr$mm = tr$m >0
tr = tr[!grepl("FALSE",tr$mm),]
tr$m = NULL
tr$mm = NULL
trait$m = match(trait$`traits$mouse_number`, tr$`trx$mouse_number`, nomatch = 0)
trait$mm = trait$m >0
trait = trait[!grepl("FALSE",trait$mm),]
trait$m = NULL
trait$mm = NULL

#Now we order the data so the rows "mouse_number" are in the same order
t = trait[order(trait$`traits$mouse_number`, decreasing=T), , drop = FALSE]
g = tr[order(tr$`trx$mouse_number`, decreasing=T), , drop = FALSE]

#Add the "mouse_number" as row.names and remove the old column
row.names(t) = t$`traits$mouse_number`
t$`traits$mouse_number` = NULL
row.names(g) = g$`trx$mouse_number`
g$`trx$mouse_number` = NULL

#We have generated 2 dataframes which both contain "mouse number" as rows and either "gene_symbol" (g) or "trait_name" (t)
#These data frames contain only "mouse_number" which match between the two
#these preprocessed examples are provided as "pre-processed_traits.txt" for traits or "pre-processed_transcripts.txt" for transcripts

#The next goal will be to use a gene list in order to generate Principle Components from transcripts
#We will use a list of candidate human GWAS genes supplied by Hooman Allayee
#To find mouse orthologues for human genes, we use the deposited list "Human_MS_Orths.txt"
HOM_MouseHumanSequence_rpt <- read.delim("Human_MS_Orths.txt", "\t", header = T)

#Now we read in the list of genes used to construct the PC - This is supplied as "Human_gwas_gene_list.txt"
new_list <- read.delim("Human_gwas_gene_list.txt", "\t", header = T)

#Add a new column which contains the mouse orthologue for each human gene_symbol
new_list$mouse_orth = HOM_MouseHumanSequence_rpt$Symbol[match(new_list$Gene.based.Analysis1, HOM_MouseHumanSequence_rpt$human_orth)]

#Generate a new data frame from "g" (all transcripts) which only contains those from the "Human_gwas_gene_list.txt" 

gg = t(g)
gg = as.data.frame(gg)
gg$m = match(row.names(gg), new_list$mouse_orth, nomatch = 0)
gg$mm = gg$m >0
gg = gg[!grepl("FALSE",gg$mm),]
gg$m = NULL
gg$mm = NULL
gg = t(gg)
gg = as.data.frame(gg)

#look at gene X gene correlation and generate heatmap of the correlation structure of the genes themselves
#This plot is shown as "Human GWAS list HMDP gene X gene correlation.pdf"
gg.cor = bicorAndPvalue(gg, gg, use = 'pairwise.complete.obs')
my_colors <- colorRampPalette(c("yellow", "grey", "blue"))(n = 299)
heatmap.2(as.matrix(gg.cor$bicor), col=my_colors, density.info="none", trace="none", cexRow =  0.1, cexCol = 0.1, dendrogram="none",Colv = T, Rowv = T, symm=F,symkey=T, symbreaks=T, scale="none", na.rm=T)

#Perform Principle component analysis on the list of genes and their correlation within the HMDP
pcca = prcomp(gg)

#The summarry will print features of each principle comonent
summary(pcca)

#The postion of each strain on gene-driven PCs is used to proceed forward
position = pcca$x

#Visualize the contribution of each gene to PC1 and PC2
#This plot is shown as "PC contribution - ALL genes.pdf"
fviz_pca_var(pcca, col.var = "contrib",
             gradient.cols = c("grey", "blue"),
             ggtheme = theme_minimal())
#We can see that PC1 and PC2 are being driven by two genes.  If needed, these can be removed from the gene expression matrix
#Remove the "outlier genes" - Note the column names are specified as the "gene_symbol" where genes "Myh11" and "Mybphl" are removed
gg[,"Myh11"] = NULL
gg[,"Mybphl"] = NULL
pcca = prcomp(gg)

#Repeat the PC analysis on outlier-removed genes
summary(pcca)
position = pcca$x

#Visualize the contribution of each gene to PC1 and PC2 after removing the two outliers.  Note the genes are more evenly spread in their contribution
#This plot is shown as "PC contribution - outliers removed.pdf"
fviz_pca_var(pcca, col.var = "contrib",
             gradient.cols = c("grey", "blue"),
             ggtheme = theme_minimal())


#Now with the "pcca" result representative vectors of PC1 and PC2 will be used for correlation
gg = as.data.frame(gg)

#assign the contribution for wach individual
contr = abs(pcca$rotation)
#use the first 2 PCs - note, can use more by speciying length of columns
comps = as.data.frame(cbind(pcca$x[,1:2]))

#construct correlation between position on PC axis vs trait
trait.pc.cor = bicorAndPvalue(comps, t, use = "p")
comp.trait.cor = as.data.frame(t(trait.pc.cor$bicor))
comp.trait.p = as.data.frame(t(trait.pc.cor$p))

#Assign colnames to eahc of the data.frames
colnames(comp.trait.cor) = c("PC1_bicor", "PC2_bicor")
colnames(comp.trait.p) = c("PC1_pvalue", "PC2_pvalue")
alternate.cols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}

#This produces a dataframe which shows the correlations between PC1 and PC2 vs traits 
#this files is listed as "PC X trait bicor and pvalue.txt"
PC_trait_cor = alternate.cols(comp.trait.cor ,comp.trait.p)
write.table(PC_trait_cor, file="PC X trait bicor and pvalue.txt",row.names=T, col.names=T, sep='\t', quote=F)

#Now generate a heatmap of each PC cxorrelated against all traits
comp.trait.cor = na.omit(comp.trait.cor)

#specify colors
my_colors <- colorRampPalette(c("red", "white", "green"))(n = 299)
heatmap.2(as.matrix(comp.trait.cor), col=my_colors, density.info="none", trace="none", cexRow =  0.1, cexCol = 0.1, dendrogram="none",Colv = T, Rowv = T, symm=F,symkey=T, symbreaks=T, scale="none", na.rm=T)

#The heatmap produced is listed as "PC1 and PC2 x trait heatmap.pdf"

