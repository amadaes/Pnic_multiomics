#if required, install and activate BiocManager
if (!require("stringr", quietly = TRUE)) install.packages("stringr")
if (!require("gsubfn", quietly = TRUE)) install.packages("gsubfn")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")

#activate general use libraries - for using strings, saving as pdf etc. 
library(stringr)
library(gsubfn)
library(data.table)
library(tidyverse)


#for csv files saved from Scaffold v5, these steps were first done in Microsoft Excel/ Libre Office: 
#removed first three columns and rows and all empty rows between Accession number and p-value column
#replace [ with nothing and ] with comma, then save as file.csv
#name Accession Number as Accession.Number
#rename p-value and FoldChange columns to padj and log2FoldChange

# #input path:
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/")
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/")

#create the list of input files
list <- list.files(path = "C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/proteomics/",
                   pattern = "*.csv",
                   full.names = T)

print(list) #check list

# #function to transform data:

for (x in list)
  {
  results = read.csv(x, sep="\t", header=TRUE)

  #save gene/ protein info from results column 1 in a vector
  info <- results[ ,1]

  #use this to get gene_ids (locus_tag)
  locus_tag <- strapplyc(info, "locus_tag=(JMY29_.*), protein=", simplify = TRUE)
  locus_tag <- gsub("\\,.*", "", locus_tag) #for protein clusters kept the first product gene in the respective list
  results$gene_id <- locus_tag

  #use this to get gene names, if available
  gene <- strapplyc(info, "gene=(.*), locus_tag=" , simplify = TRUE)
  results$gene_name<-gene

  #use this to get protein names
  protein <- strapplyc(info, "(protein=.*), protein_id", simplify = TRUE)
  head(protein)
  protein <- sub("protein=", "", protein)
  results$protein_name<-protein

  #use this to get protein ids
  protein_id <- strapply(info, "(protein_id=.*), location", simplify = TRUE)
  head(protein_id)
  protein_id <- sub("protein_id=", "", protein_id)
  results$protein_id <- protein_id

  #use this to get chromosome/ pAO1 location
  location <- strapply(results$Accession.Number, "(CP08929.*)_prot", simplify = TRUE)
  head(location)
  location <- sub("CP089293.1", "chrs", location)
  location <- sub("CP089294.1", "pAO1", location)
  # head(location)
  results$location <- location

  better_results <- results[,-1:-2]

  # this adds the protein name to the output data frame
  final_results <- better_results[,c(9,1:2)]
  final_results$padj <- gsub("< ", "", final_results$padj)

  final_results_ordered <- final_results[order(final_results$padj),]
  significant_results <- subset(final_results_ordered, padj<=0.1)
  significant_results$padj <- gsub("< ", "", significant_results$padj)
  rescaled <- significant_results

  # #rescale FC values to match scale for RNA-Seq LFC:
  # max values are [-8, 6] for RNA-seq, but the scale will go up to 8 in order to have 0 as the median value
  
  # #assign value of 100000000 to FC=inf
  rescaled$log2FoldChange<-replace(rescaled$log2FoldChange, rescaled$log2FoldChange==Inf, 100000000)
  
  # #assign 100 to insignificant FC (0.5, 2) => -1, 1
  rescaled <- rescaled %>%
    mutate(across(.cols = c(3), .fns = function(y) ifelse((y > 0.5 & y < 2), 100, y)))

  # #assign LFC to non-0 and significantly upregulated FC values and rescale to match RNA-SEQ scale
  # #max value for MS/MS is FC = 19; LOG2(19) = 4.24792751344;
  # #TO SCALE UP TO MIN LFC FROM RNA-SEQ, WHICH IS 6, I DID 6:4.24792751344 = 1.412453480201022

  rescaled <- rescaled %>%
    mutate(across(.cols = c(3), .fns = function(y) ifelse((y!=0 & y >= 2 & y<100), log2(y)*1.412453480201022, y)))

  # #assign max LFC value (8) to FC=INF
  rescaled <- rescaled %>%
    mutate(across(.cols = c(3), .fns = function(y) ifelse(y == 100000000, 8, y)))

  # #assign LFC to non-0 and significantly downregulated FC values and rescale to match RNA-SEQ scale
  # #min value for MS/MS is FC = 0.08; LOG2(0.08) = -3.64385618977;
  # #TO SCALE UP TO MIN LFC FROM RNA-SEQ, WHICH IS -8, used -8:-3.64385618977= 2.195476326003129
  rescaled <- rescaled %>%
    mutate(across(.cols = c(3), .fns = function(y) ifelse((y!=0 & y <= 0.5), log2(y)*2.195476326003129, y)))

   # #assign 0.05 to insignificant FC (0.5, 2)
  rescaled <- rescaled %>%
    mutate(across(.cols = c(3), .fns = function(y) ifelse(y==100, 0.05, y)))

  # #assign -8 to FC equal to 0
  rescaled <- rescaled %>%
    mutate(across(.cols = c(3), .fns = function(y) ifelse(y==0, -8, y)))

  # #change gene_id column to match input for functional analysis script
  colnames(rescaled)[which(names(rescaled) == "gene_id")] <- ""

  # #remove padj column for pathview input
  rescaled_no_padj <- rescaled[,-2]

  # #save files as csv in working directory
  fwrite(rescaled_no_padj, sep=',', gsub('.csv', '_rescaled.csv', x))

}
