if (!require("BiocManager", quietly = FALSE)) install.packages("BiocManager")
BiocManager::install("AnnotationHub")
BiocManager::install("AnnotationDbi")
BiocManager::install("stats4")
BiocManager::install("IRanges")
BiocManager::install("S4Vectors")
BiocManager::install("AnnotationForge")
BiocManager::install("biomaRt") #required for makeOrgPackageFromNCBI - processing GO data
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("ggplot2")
install.packages("Rtools")
library(dplyr)

# CREATE AN ORGANISM PACKAGE -----------------------------------------------------------
# https://bioconductor.org/packages/devel/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html

library(AnnotationDbi)
library(AnnotationForge)

# SPECIES: Paenarthrobacter nicotinovorans ATCC 49919, NCBI TAXON ID 29320

makeOrgPackageFromNCBI(version = "0.1",
                       author = "Amada El-Sabeh <amadaelsabeh@gmail.com>",
                       maintainer = "Amada El-Sabeh <amadaelsabeh@gmail.com>",
                       outputDir = ".",
                       tax_id = "29320",
                       genus = "Paenarthrobacter",
                       species = "nicotinovorans",
                       rebuildCache=FALSE)

# created a package with the db for Pnic, now we will install this package for the first use
install.packages("D:/Activitati/Doctorat/4_REFERAT2/DESeq/org.Pnicotinovorans.eg.db", repo = NULL, type = "source")

# call the library for subsequent uses of org.Pnicotinovorans.eg.db
library(org.Pnicotinovorans.eg.db)

## ---- Inspect content of OrgPackage -----------------------------------------------------------

# inspect the object (output in console)
org.Pnicotinovorans.eg.db 

# keytype - This is the source of the annotation (gene ids). The options vary for each annotation. 
# check which options are available with the following 

keytypes(org.Pnicotinovorans.eg.db) 

# (output in console):
# [1] "ACCNUM"      "ALIAS"       "ENTREZID"    "EVIDENCE"    "EVIDENCEALL" "GENENAME"    "GID"         "GO"          "GOALL"       "ONTOLOGY"    "ONTOLOGYALL"
# [12] "PMID"        "REFSEQ"      "SYMBOL"  
 
# #extract some sample keys of a particular type using keys and head to actually look at the output in the console
# #ALIAS: the locus tags: JMY_...; GENENAME: actual gene names; GO: Gene Ontology ID etc.
# head(keys(org.Pnicotinovorans.eg.db, keytype="ALIAS"))

# GENE SET ENRICHMENT ANALYSIS (GSEA) FOR INDIVIDUAL DGE RESULTS FILES -----------------------------------------------------------

## ----- Create and set output directory and read input  -----------------------------------------------------------
library(stringr)
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/")
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_single_DE_file/")
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_single_DE_file/nic_T2_vs_nic_T1/") #no sDEGs 
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_single_DE_file/nic_T3_vs_nic_T1/")
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_single_DE_file/nic_T3_vs_nic_T2/")

# # one at a time, read as input a DESeq2 results file (output of Pnic_DESeq2_analysis.R) and set corresponding wd

df = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_nic_T2_vs_nic_T1.csv", sep=",") #no sDEGs
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_single_DE_file/nic_T2_vs_nic_T1/")

# df = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_nic_T3_vs_nic_T1.csv", sep=",")
# setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_single_DE_file/nic_T3_vs_nic_T1/")

# df = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_nic_T3_vs_nic_T2.csv", sep=",")
# setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_single_DE_file/nic_T3_vs_nic_T2/")

## ----- Prep list of LFC and names of sDEGs from individual RNA-seq DESeq2 results  -----------------------------------------------------------
# # #https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/

#filter significantly DEGs by padj < thresh and |log2FC| > 1
significant <- subset(df, padj<0.1 & abs(log2FoldChange)>1)
# we want the log2 fold change; the values are saved as a vector
original_gene_list <- as.numeric(significant$log2FoldChange)
# name the vector (give gene ids back to the log2 fold change values)
names(original_gene_list) <- significant$X #the same as ALIAS from org.Pnicotinovorans.eg.db
# omit any NA values
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
# Extract significant results (padj < 0.1)
sig_genes = subset(significant, padj < 0.1)
# From significant results, we want to extract the foldchanges
genes <- sig_genes$log2FoldChange
# Name the vector (each fold change with the corresponding ID)
names(genes) <- sig_genes$X
# omit NA values
genes <- na.omit(genes)
# sort the list in decreasing order (required for clusterProfiler)
genes <- sort(genes, decreasing = TRUE)
# filter on min log2fold change (log2FoldChange > 1)
genes <- names(genes)[abs(genes) > 1]

## ----- Actual GSEA -----------------------------------------------------------

# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
##The clusterProfiler package supports all organisms that have KEGG annotation data available in the KEGG database. 
#to search KEGG db for Pnic: search_kegg_organism('pnv', by='kegg_code')
#we'll get one hit,  kegg_code  6953               scientific_name   pnv Paenarthrobacter nicotinovorans 

library(clusterProfiler)

gse <- gseGO(geneList=gene_list,
             ont ="ALL",
             keyType = "ALIAS",
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = "org.Pnicotinovorans.eg.db",
             pAdjustMethod = "BH")

# # GSE OUPUT - dotplot:
require(DOSE)
library(ggplot2)
library(Cairo)

png(file="dose_dotplot.png", res=100)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

# ------ PATHWAY OVER-REPRESENTATION ANALYSIS FOR MERGED DATA -----------------------------------------------------------
# #we'll generate an image/ pathway for nic vs control, respectively nic vs nic: 

## ----- Prep list of LFC and names of sDEGs & sDEPs from RNA-seq and MS/MS data: NIC VS NIC -----------------------------------------
library("dplyr")

###Process data for nic vs nic: first the three DESeq2 Results files, followed by the three lists for MS/MS
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_merged_data/")
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_nic_merged_data/")

# #start with RNA-seq data:
# #reading input from deseq2 - Wald test results for all relevant comparisons
df1 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_nic_T2_vs_nic_T1.csv", sep=",")
df1 <- df1[,c(1,3,7)]

df2 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_nic_T3_vs_nic_T1.csv", sep=",")
df2 <- df2[,c(1,3,7)]

df3 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_nic_T3_vs_nic_T2.csv", sep=",")
df3 <- df3[,c(1,3,7)]

rna_all_df <- merge(df1, df2, by="X", all=TRUE)
rna_all_df <- merge(rna_all_df, df3, by="X", all=TRUE)
make.unique(colnames(rna_all_df))
colnames(rna_all_df) <- c("X", "log2FoldChange.1", "padj.1", "log2FoldChange.2", "padj.2", "log2FoldChange.3", "padj.3")

significant_rna <- subset(rna_all_df, (padj.1 <0.1 & abs(log2FoldChange.1)>=1) | (padj.2 <0.1 & abs(log2FoldChange.2)>=1) )
# #REMOVE PADJ COLUMNS
significant_rna <- significant_rna[,c(1,2,4,6)]
# #sort the list in decreasing order (required for clusterProfiler)
significant_rna <- arrange(significant_rna, desc(log2FoldChange.1), desc(log2FoldChange.2), desc(log2FoldChange.3))
rescaled <- significant_rna

# #assign value of 0 to -1< LFC <1 (not significant LFC)
rescaled <- rescaled %>%
  mutate(across(.cols = c(2,3,4), .fns = function(y) ifelse(abs(y) <= 1, 0, y)))

# #now MS/MS data:
df.p1 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/nic_T2_vs_T1_rescaled.csv", sep=",")
df.p2 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/nic_T3_vs_T1_rescaled.csv", sep=",")
df.p3 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/nic_T3_vs_T2_rescaled.csv", sep=",")

prot_all_df <- merge(df.p1, df.p2, by="X", all=TRUE)
prot_all_df <- merge(prot_all_df, df.p3, by="X", all=TRUE)
significant_prot <- prot_all_df

# #now MERGE RNA-seq and MS/MS data for nic vs nic
significant_all <- merge(rescaled, significant_prot, by="X", all=TRUE)
# #remove duplicated rows
significant_all <- significant_all[!duplicated(significant_all), ]
original_gene_list_all <- significant_all[,2:7]
original_gene_list_all <- original_gene_list_all[,c("log2FoldChange.1", "log2FoldChange.x", "log2FoldChange.2", "log2FoldChange.y", "log2FoldChange.3", "log2FoldChange" )]
rownames(original_gene_list_all) <- significant_all$X #the same as ALIAS from org.Pnicotinovorans.eg.db
gene_list<-original_gene_list_all
genes <- rownames(original_gene_list_all)
genes = sort(genes, decreasing = TRUE) 

## ----- Now, for Pathway enrichment, executed the code under the header ## ----- ENRICHED KEGG PATHWAY ANALYSIS & VIZUALISATION -----------------------------------------


# #After finishing with the first comparison set, move on to:

## ----- Prep list of LFC and names of sDEGs & sDEPs from RNA-seq and MS/MS data: NIC VS CONTROL -----------------------------------------
## OVERWRITES FILES GENERATED FOR NIC VS NIC IF WD IS NOT CHANGED

###Process data for nic vs control: first the three DESeq2 Results files, followed by the three lists for MS/MS
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_no_nic_merged_data/")
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_funct_analysis/test/nic_vs_no_nic_merged_data/")

# #start with RNA-seq data:
df1 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_T1_nic_vs_no_nic.csv", sep=",")
df1 <- df1[,c(1,3,7)]

df2 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_T2_nic_vs_no_nic.csv", sep=",")
df2 <- df2[,c(1,3,7)]

df3 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_geneids/tables/res_geneids_T3_nic_vs_no_nic.csv", sep=",")
df3 <- df3[,c(1,3,7)]

rna_all_df <- merge(df1, df2, by="X", all=TRUE)
rna_all_df <- merge(rna_all_df, df3, by="X", all=TRUE)
make.unique(colnames(rna_all_df))
colnames(rna_all_df) <- c("X", "log2FoldChange.1", "padj.1", "log2FoldChange.2", "padj.2", "log2FoldChange.3", "padj.3")

significant_rna <- subset(rna_all_df, (padj.1 <0.1 & abs(log2FoldChange.1)>=1) | (padj.2 <0.1 & abs(log2FoldChange.2)>=1) )
# rna_all_df <- rna_all_df[,-c(4, 7, 10, 13, 16)]
# #REMOVE PADJ COLUMNS
significant_rna <- significant_rna[,c(1,2,4,6)]
# #sort the list in decreasing order (required for clusterProfiler)
significant_rna <- arrange(significant_rna, desc(log2FoldChange.1), desc(log2FoldChange.2), desc(log2FoldChange.3))
rescaled <- significant_rna

# #assign value of 0 to -1< LFC <1 (not significant LFC)
rescaled <- rescaled %>%
  mutate(across(.cols = c(2,3,4), .fns = function(y) ifelse(abs(y) <= 1, 0, y)))

# #now MS/MS data:
df.p1 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/T1_nic_vs_no_nic_rescaled.csv", sep=",")
df.p2 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/T2_nic_vs_no_nic_rescaled.csv", sep=",")
df.p3 = read.csv("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/T3_nic_vs_no_nic_rescaled.csv", sep=",")

prot_all_df <- merge(df.p1, df.p2, by="X", all=TRUE)
prot_all_df <- merge(prot_all_df, df.p3, by="X", all=TRUE)
significant_prot <- prot_all_df

# #now MERGE RNA-seq and MS/MS data for nic vs control
significant_all <- merge(rescaled, significant_prot, by="X", all=TRUE)
# #remove duplicated rows
significant_all <- significant_all[!duplicated(significant_all), ]
original_gene_list_all <- significant_all[,2:7]
original_gene_list_all <- original_gene_list_all[,c("log2FoldChange.1", "log2FoldChange.x", "log2FoldChange.2", "log2FoldChange.y", "log2FoldChange.3", "log2FoldChange" )]
rownames(original_gene_list_all) <- significant_all$X #the same as ALIAS from org.Pnicotinovorans.eg.db
gene_list<-original_gene_list_all
genes <- rownames(original_gene_list_all)
genes = sort(genes, decreasing = TRUE)

## Now, for Pathway enrichment, executed the code under the header ## ----- ENRICHED KEGG PATHWAY ANALYSIS & VIZUALISATION -----------------------------------------

# # EXECUTE ONCE for nic vs control, ONCE for nic vs nic

## ----- ENRICHED KEGG PATHWAY ANALYSIS & VIZUALISATION -----------------------------------------------------------
# Go here and select prefix: pnv, then search pathways https://www.genome.jp/kegg/pathway.html 
# available pathways for Pnic: https://www.kegg.jp/kegg-bin/search_pathway_text?map=pnv&keyword=&mode=1&viewImage=true
# CHECKED KEGG_enrich.csv to see which pathways are over-represented

## actual enrichment analysis:
library(clusterProfiler)

KEGG_enrich <- enrichKEGG(gene = genes, organism = "pnv", pvalueCutoff = 0.1, pAdjustMethod = "BH")
write.csv(KEGG_enrich, file="KEGG_enrich.csv")

## visualize KEGG Pathways:

library("pathview")

# #make a list of enriched pathways
sig_paths <- KEGG_enrich$ID

#legend positions depending on pathway
key.positions <- setNames(c("bottomright", "bottomright", "bottomright", "topright", "topleft", "bottomright", "topright", "topright", "topright", "topright", "bottomright", "bottomright", "topright", "topright", "bottomright", "bottomright", "bottomright", "bottomleft", "topright", "topright", "topright", "topright"), 
                          c("pnv00760", "pnv00650", "pnv00220", "pnv02020", "pnv01120", "pnv00790", "pnv00860", "pnv01240", "pnv00290", "pnv03018", "pnv03010", "pnv01230", "pnv00190", "pnv00020", "pnv01110", "pnv03070", "pnv00920", "pnv01210", "pnv00010", "pnv00250", "pnv00270", "pnv00230"))

# #loop through list of significant pathways and generate figure for each one:
for (i in sig_paths)
{
  skip_to_next <- FALSE

  tryCatch(
    pathview_res <- pathview(gene.data  = gene_list[,1:6],
                             pathway.id = i,
                             species    = "pnv",
                             gene.idtype = "KEGG",
                             kegg.native = T,
                             match.data = FALSE,
                             keys.align = "y",
                             key.pos = key.positions[i],
                             same.layer = T,
                             both.dirs = TRUE,
                             bins = 16,
                             limit = c(-8, 8),
                             high = "red", mid = "gray", low = "green",
                             # #na.col = "gray",
                             cex = 0.15,
                             out.suffix = "gene_ids_RNA_MS"),

    error = function(e)     { skip_to_next <<- TRUE}) 
    if(skip_to_next)  { next } # #some KEGG pathways are not mappable, this skips the null output
}

# graph for the P. nicotinovorans ATCC 49919 nicotine pathway - pnv00760
pathview_res <- pathview(gene.data  = gene_list[,1:6],
                         pathway.id = "pnv00760",
                         species    = "pnv",
                         gene.idtype = "KEGG",
                         kegg.native = T,
                         match.data = FALSE,
                         keys.align = "y",
                         key.pos = key.positions["pnv00760"],
                         same.layer = T,
                         both.dirs = TRUE,
                         bins = 16,
                         limit = c(-8, 8),
                         high = "red", mid = "gray", low = "green",
                         # #na.col = "gray",
                         cex = 0.15,
                         out.suffix = "gene_ids_RNA_MS")
