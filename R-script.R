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

new_list <- read.delim("C:/Users/mseldin/Desktop/lab files/Hooman Lab/Gene List for Marcus-Calvin.txt", header = T)

new_list$mouse_orth = HOM_MouseHumanSequence_rpt$Symbol[match(new_list$Gene.based.Analysis1, HOM_MouseHumanSequence_rpt$human_orth)]


gg = t(g)
gg = as.data.frame(gg)

gg$m = match(row.names(gg), new_list$mouse_orth, nomatch = 0)
gg$mm = gg$m >0
gg = gg[!grepl("FALSE",gg$mm),]
gg$m = NULL
gg$mm = NULL
gg = t(gg)
gg = as.data.frame(gg)
gg$Myh11 = NULL
gg$Mybphl = NULL

#look at gene X gene correlation

gg.cor = bicor(gg, gg, use = 'pairwise.complete.obs')

my_colors <- colorRampPalette(c("yellow", "grey", "blue"))(n = 299)
heatmap.2(as.matrix(gg.cor), col=my_colors, density.info="none", trace="none", cexRow =  0.1, cexCol = 0.1, dendrogram="none",Colv = T, Rowv = T, symm=F,symkey=T, symbreaks=T, scale="none", na.rm=T)

pcca = prcomp(gg)
summary(pcca)
position = pcca$x



fviz_pca_var(pcca, col.var = "contrib",
             gradient.cols = c("grey", "blue"),
             ggtheme = theme_minimal())

#print variance explained


contr = abs(pcca$rotation)

comps = as.data.frame(cbind(pcca$x[,1:3]))
trait.pc.cor = bicorAndPvalue(comps, t, use = "p")
comp.trait.cor = as.data.frame(t(trait.pc.cor$bicor))
comp.trait.p = as.data.frame(t(trait.pc.cor$p))

write.table(comp.trait.cor, file="PC X trait bicor",row.names=T, col.names=T, sep='\t', quote=F)
write.table(comp.trait.p, file="PC X trait pvalue",row.names=T, col.names=T, sep='\t', quote=F)
write.table(contr[,1:3], file="PC gene contributions",row.names=T, col.names=T, sep='\t', quote=F)
