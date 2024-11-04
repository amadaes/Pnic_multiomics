Python3 and R Scripts used for the multiomic analysis of Paenarthobacter nicotinovorans ATTC 49919, starting from direct-RNA sequencing and nanoLC-MS/MS data.
A raw count matrix for differential gene expression analysis was obtained using StringTie2 and featureCounts, integrated into the nf-core/nanoseq version 3 pipeline. 

**- replace_stringtie_ID_merge_counts.py**
    - Python3 script for replacing the MSTRG IDs introduced by StringTie with the reference IDs/ names and combining them with an existing count matrix
    - this script can be used, for example, after running nf-core/nanoseq or rnaseq etc. to obtain a single table containing the reference IDs and the count data
    - the script can also be easily adjusted for processing any gtf/ gff file, to split and extract the information of interest from the features column
       
**- Pnic_DESeq2_analysis.R**
    - read in count matrix from featureCounts (countData) and create colData file needed for DESeqDataSetFromMatrix
    - perform sample quality assurance using the regular log transformation and the limma removeBatchEffect function to control for variability between replicates (generate PCA plots, sample heatmaps etc.)
    - execute DGE analysis and generate results from RNA-seq data using DESeq2
    - generate plots for visualizing DE results (dispersion plots, volcano plots etc.)
      
**- Pnic_proteome_LFC_rescale.R**
    - R script used to batch process differential protein expression analysis results (here, output from Scaffold software)
    - reads multiple .csv files from a list
    - extracts gene/ protein IDs and names from a comma-separated field and adds them to a new data frame 
    - mutates/ re-scales the fold change values to be mappable on the same value scale as RNA-seq data (for generating integrated figures with Pathview using the script Pnic_functional_enrichment.R )
    - output files as .csv
         
**- Pnic_functional_enrichment.R**
    - create a custom AnnotationDBI database using NCBI Annotation package
    - read data frames containing DEGs (DESeq2 results - output of Pnic_DESeq2_analysis.R) and DEP (Scaffold results - output of Pnic_proteome_LFC_rescale.R)
    - generate the genes and gene_list objects necessary for functional enrichment
    - perform gene set and pathway enrichment analysis using ClusterProfiler (DOSE, KEGG etc.) from Bioconductor
    - use for loop to generate all the figures for significantly enriched pathways (for combined DEG and DEP data) using the Pathview package from Bioconductor
