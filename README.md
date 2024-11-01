**- replace_stringtie_ID_merge_counts.py** Python3 script for replacing the MSTRG IDs introduced by StringTie with the reference IDs/ names and combining them with an existing count matrix. 
  This script can be used, for example, after running nf-core/nanoseq or rnaseq etc. to obtain a single table containing the reference IDs and the count data. 
  The script can also be easily adjusted for processing any gtf/ gff file and to split and extract the information of interest from the features column.

**- Pnic_DESeq2_analysis.R** R script used for differential gene expression analysis and visualisation, mainly based on the DESeq2 Bioconductor package. 

**- Pnic_proteome_LFC_rescale.R** R script used to batch process differential protein expression analysis results files (output from Scaffold software). Reads multiple .csv files from a list, processes and outputs the modifications as .csv.

**- Pnic_functional_enrichment.R** R script used to generate a custom annotation database and perform functional enrichment analyses of merged RNA-Seq and MS/MS data, mainly based on the ClusterProfiler and Pathview packages available in Bioconductor. 
