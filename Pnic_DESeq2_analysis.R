#activate BiocManager, Biobase, BioGenerics, IRanges, S4Vectors
#Install DESeq2, apeglm, vsn packages from Bioconductor if required (not yet installed)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
if (!require("stringr", quietly = TRUE)) install.packages("stringr")
if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2") 
if (!require("vsn", quietly = TRUE)) BiocManager::install("vsn")
if (!require("limma", quietly = TRUE)) BiocManager::install("limma")
if (!require("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!require("Cairo", quietly = TRUE)) install.packages("Cairo")

library("tidyverse")
library("Cairo")

#Create and set output directory
library(stringr)
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing")

# Import & pre-process featureCounts raw counts (countData) and sample info (colData) -----------------------------------------------------
##import count table from featurecounts and assign it to countdata
countdata <- data.frame(read.table("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/used_data/RNA_Seq/nf-core_nanoseq_all_samples/ids_and_counts.txt", sep=",", header=TRUE, skip = 0))
##keep only the column with Genenames or Geneids and counts from countdata
##use Genenames for DE analysis and Geneids (KEGG IDs) for functional enrichment
countdata$Chr <- countdata$Start <- countdata$End <- countdata$Length <- countdata$Strand <- NULL
##remove .sorted.bam or .sorted.sam from filenames
colnames(countdata)[2:length(colnames(countdata))] <- unlist(lapply(strsplit(colnames(countdata)[2:length(colnames(countdata))], "\\."), "[[", 1))
countdata <- aggregate(countdata[, -1:-2],countdata["Genename"],sum)
##create a vector to keep the group and sample order  
sample.order<-c("Genename", 
                "no_nic_T1_R1", "no_nic_T1_R2", "no_nic_T1_R3", "no_nic_T1_R4", "no_nic_T1_R5","no_nic_T1_R6", "no_nic_T1_R7",
                "no_nic_T2_R1", "no_nic_T2_R2", "no_nic_T2_R3", "no_nic_T2_R4", "no_nic_T2_R5","no_nic_T2_R6", "no_nic_T2_R7",
                "no_nic_T3_R1", "no_nic_T3_R2", "no_nic_T3_R3", "no_nic_T3_R4", "no_nic_T3_R5","no_nic_T3_R6", "no_nic_T3_R7",
                "nic_T1_R8", "nic_T1_R9", "nic_T1_R10", "nic_T1_R11", "nic_T1_R12",
                "nic_T2_R8", "nic_T2_R9", "nic_T2_R10", "nic_T2_R11", "nic_T2_R12",
                "nic_T3_R8", "nic_T3_R9", "nic_T3_R10", "nic_T3_R11", "nic_T3_R12",
                "pA_R1", "pA_R2", "pA_R3",
                "TEX_R1", "TEX_R2", "TEX_R3")
##order sample columns as control (1,2,3,4,5 etc.) and experiment (1,2,3,4,5 etc.) using sample.order
countdata <- countdata[ ,sample.order]
##remove pA and TEX samples
countdata <- subset(countdata, select=-c(pA_R1, pA_R2, pA_R3, TEX_R1, TEX_R2, TEX_R3))
##remove Genename column from countdata table, as required for the DESeqDataSetFromMatrix input 
countdata_no_ids <- countdata[, -1]
rownames(countdata_no_ids) <-countdata[, 1]
##creating the coldata file: we need sample names as the first column and group classification as the second column
##use a vector to keep sample names from countdata table to add to coldata
sample <- colnames(countdata_no_ids)
##use a vector to get the group names from the sample names (removes _R1, _R2 etc. from the sample names)
group <- sub("(^[^-]+)_.*", "\\1", sample)
condition <- sub("(^[^-]+)_.*", "\\1", group)
batch <- str_sub(sample, start=-3)
batch <- sub("_", "", batch)
time <- str_sub(group[1:36], start= -2)
length(time) <- length (sample)
#create the coldata dataframe using the data from the sample and group vectors
coldata <- data.frame(group, condition, time, batch, row.names = sample)
# coldata <- data.frame(group, batch, row.names = sample)
##check that order of samples from coldata matches order of samples in count table; if not, columns will be arranged to match between files
if (!all(rownames(coldata) == colnames(countdata_no_ids))){coldata <- coldata[match(colnames(countdata_no_ids), rownames(coldata)), ]}

# Analysis with DESeq2------------------------------------------------------------------
library("DESeq2")

## ----experimental design factors------------------------------------------------------------------
coldata$condition <- factor(coldata$condition)
coldata$time <- factor(coldata$time)
coldata$group <- factor(coldata$group)
coldata$batch <- factor(coldata$batch)

## ---- create the DESeqDataSetFromMatrix------------------------------------------------------------------
# dds <- DESeqDataSetFromMatrix(countData = countdata_no_ids, colData = coldata, design = ~ condition + time + condition:time)
dds <- DESeqDataSetFromMatrix(countData = countdata_no_ids, colData = coldata, design = ~ group)

## ----factorlvl
# #by default, factor levels are alphabetically ordered; need to specify wanted order 
dds$condition <- factor(dds$condition, levels = c("no_nic","nic"))
dds$time <- factor(dds$time, levels = c("T1", "T2", "T3"))
dds$group <- factor(dds$group, levels = c("no_nic_T1", "no_nic_T2", "no_nic_T3", "nic_T1", "nic_T2", "nic_T3"))
print(dds$group)
dds$batch <- factor(dds$batch, levels = c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12"))
print(dds$batch)

# SAMPLE QUALITY CHECKS: all samples------------------------------------------------------------------
#create directory in which to save initial sample QC results
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/all_samples_QC/")
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/all_samples_QC/")

## ---- Regularized log transformation (rlog/rlogTransformation) for data visualization----------------------------------------------------
#if blind=FALSE rlog uses design formula to calculate the within-group variability (HERE = BIOLOGICAL REPLICATES)
#(https://rdrr.io/bioc/DESeq2/man/rlog.html)

rld_blind_F <- rlogTransformation(dds, blind=FALSE) # this will be later overwritten
rld_blind_F_no_limma <- rld_blind_F

#if blind=TRUE rlog calculates the across-all-samples variability
rld_blind_T <- rlogTransformation(dds, blind=TRUE)

#to visualize the transformed data with batch variation removed, use the removeBatchEffect function from limma
#this simply removes any shifts in the log2-scale expression data that can be explained by batch
#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#control-features-for-estimating-size-factors

library("limma")

mat <- assay(rld_blind_F)
mm <- model.matrix(~time+condition+condition:time, colData(rld_blind_F))
mat <- limma::removeBatchEffect(mat, batch=rld_blind_F$batch, design=mm)
assay(rld_blind_F) <- mat

## ----Effects of transformations on the variance
# The figure below plots the standard deviation of the transformed data, across samples, against the mean, using the rlog transformation
# Note that the vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions.
# While a flat curve of the square root of variance over the mean may seem like the goal of such transformations, 
# this may be unreasonable in the case of data sets with many true differences due to the experimental conditions.

library("vsn")

#this one shows effect of transformation using the experimental design without limma to remove batch variation
grDevices::cairo_pdf(file="stdev_rld_blind_F_no_limma_vs_mean.pdf")
meanSdPlot(assay(rld_blind_F_no_limma))
dev.off()

#this one shows effect of transformation using the experimental design with limma to remove batch variation
grDevices::cairo_pdf(file="stdev_rld_blind_F_vs_mean.pdf")
meanSdPlot(assay(rld_blind_F))
dev.off()

#this one shows effect of transformation across all samples without accounting for experimental design
grDevices::cairo_pdf(file="stdev_rld_blind_T_vs_mean.pdf")
meanSdPlot(assay(rld_blind_T))
dev.off()

## ---- Sample distance heatmaps ------------------------------------------------

# Colours for plots below
##use RColorBrewer or any from https://github.com/EmilHvitfeldt/r-color-palettes
#library("RColorBrewer") or make a vector for color values wanted (can also make a gradient using https://mycolor.space/gradient3
# #ff0000, #ff5400, #ff7b00, #ff9b19, #ffb73a, #ffab3f, #ff9f44, #ff944a, #ff5658, #ff0084, #d300c2, #0432ff
mycols <- c("#ff0000", "#FFD966","#0432FF", "#3b3b3b", "#777777", "#b9b9b9")

#this one for sample distance calculated using rlog transformation using the experimental design (blind=FALSE) without limma to remove batch variation
library(gplots)

sampleDists <- as.matrix(dist(t(assay(rld_blind_F_no_limma))))
ncol(sampleDists)

grDevices::cairo_pdf(file="qc-heatmap-samples_rld_blind_F_no_limma.pdf", width=20, height=20, pointsize=25)
heatmap.2(as.matrix(sampleDists), key=T, key.xlab="Value", key.ylab="count", symkey=FALSE, keysize= 0.2, trace="none", key.par = list(cex=0.8), lhei=c(1.5,3), lwid=c(3,6),
          col=colorpanel(50, "black", "white"),
          margin=c(10, 10))
dev.off()

#this one for sample distance calculated using rlog transformation using the experimental design to remove batch variation (blind=FALSE)

sampleDists <- as.matrix(dist(t(assay(rld_blind_F))))
ncol(sampleDists)

grDevices::cairo_pdf(file="qc-heatmap-samples_rld_blind_F.pdf", width=20, height=20, pointsize=25)
heatmap.2(as.matrix(sampleDists), key=T, key.xlab="Value", key.ylab="count", symkey=FALSE, keysize= 0.2, trace="none", key.par = list(cex=0.8), lhei=c(1.5,3), lwid=c(3,6),
          col=colorpanel(50, "black", "white"),
          margin=c(10, 10))
dev.off()

#this one for sample distance calculated using rlog transformation without accounting for design (blind=TRUE)

sampleDists <- as.matrix(dist(t(assay(rld_blind_T))))
ncol(sampleDists)

grDevices::cairo_pdf(file="qc-heatmap-samples_rld_blind_T.pdf", width=20, height=20, pointsize=25)
heatmap.2(as.matrix(sampleDists), key=T, key.xlab="Value", key.ylab="count", symkey=FALSE, keysize= 0.2, trace="none", key.par = list(cex=0.8), lhei=c(1.5,3), lwid=c(3,6),
          col=colorpanel(50, "black", "white"),
          margin=c(10, 10))
dev.off()

## ---- Principal components analysis (PCA) ----------------------------------------

## used a method from Stephen Turner, @genetics_blog

#using rlog transformation using the experimental design to remove batch variation (blind=FALSE) without removing batch effect with limma

png(file="qc-pca_rld_blind_F_no_limma.png", 1500, 1500, pointsize=20)
rld_pca <- function (rld_blind_F_no_limma, intgroup = c("condition","time"), ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, repel=TRUE, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  require(DESeq2)
  rv = rowVars(assay(rld_blind_F_no_limma))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld_blind_F_no_limma)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld_blind_F_no_limma)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "#ff0000")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
rld_pca(rld_blind_F_no_limma, colors=mycols, intgroup=c("condition","time"), xlim=c(-50, 50))
dev.off()

#using rlog transformation using the experimental design to remove batch variation (blind=FALSE)

png(file="qc-pca_rld_blind_F.png", 1500, 1500, pointsize=20)
rld_pca <- function (rld_blind_F, intgroup = c("condition","time"), ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  require(DESeq2)
  rv = rowVars(assay(rld_blind_F))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld_blind_F)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld_blind_F)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "#ff0000")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
rld_pca(rld_blind_F, colors=mycols, intgroup=c("condition","time"), xlim=c(-50, 50))
dev.off()

#using rlog transformation without accounting for design (blind=TRUE)
png(file="qc-pca_rld_blind_T.png", 1500, 1500, pointsize=20)
rld_pca <- function (rld_blind_T, intgroup = c("condition","time"), ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  require(DESeq2)
  rv = rowVars(assay(rld_blind_T))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld_blind_T)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld_blind_T)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "#ff0000")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
rld_pca(rld_blind_T, colors=mycols, intgroup=c("condition","time"), xlim=c(-50, 50))
dev.off()

# REMOVE OUTLIERS?------------------------------------------------------------------

#check the PCA and sample distance heatmaps and decide if there are sample outliers
#extract outliers from the countdata and coldata before continuing with analyzing DE
print(rownames(colData(dds)))
#removed from colData outlier samples (by specifying the rownumber from colData)
#outlier sample: no_nic_T2_R5
# dds <- dds[,-c(11,12,18)]
dds <- dds[,-12]

#check the samples kept in colData
print(rownames(colData(dds)))

# SAMPLE QUALITY CHECKS: kept samples------------------------------------------------------------------

#now we generate the PCA and sample distance heatmaps with the new dds object, which contains only samples to be kept for DE analysis etc.
# #create directory in which to save kept sample QC results
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_QC/")
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_QC/")

## ---- Regularized log transformation (rlog/rlogTransformation) for data visualization----------------------------------------------------
#if blind=FALSE rlog uses design formula to calculate the within-group variability (HERE = BIOLOGICAL REPLICATES)
#(https://rdrr.io/bioc/DESeq2/man/rlog.html)

rld_blind_F <- rlogTransformation(dds, blind=FALSE) # this will be later overwritten
rld_blind_F_no_limma <- rld_blind_F

#if blind=TRUE rlog calculates the across-all-samples variability
rld_blind_T <- rlogTransformation(dds, blind=TRUE)

#to visualize the transformed data with batch variation removed, use the removeBatchEffect function from limma
#this simply removes any shifts in the log2-scale expression data that can be explained by batch
#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#control-features-for-estimating-size-factors
#Why after VST are there still batches in the PCA Biplot?
library("limma")

mat <- assay(rld_blind_F)
mm <- model.matrix(~time+condition+condition:time, colData(rld_blind_F))
mat <- limma::removeBatchEffect(mat, batch=rld_blind_F$batch, design=mm)
assay(rld_blind_F) <- mat

## ----Effects of transformations on the variance
# The figure below plots the standard deviation of the transformed data, across samples, against the mean, using the rlog transformation
#Note that the vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions.
#While a flat curve of the square root of variance over the mean may seem like the goal of such transformations, this may be unreasonable in the case of datasets with many true differences due to the experimental conditions.

library("vsn")

#this one shows effect of transformation using the experimental design without limma to remove batch variation
grDevices::cairo_pdf(file="stdev_rld_blind_F_no_limma_vs_mean.pdf")
meanSdPlot(assay(rld_blind_F_no_limma))
dev.off()

#this one shows effect of transformation using the experimental design with limma to remove batch variation
grDevices::cairo_pdf(file="stdev_rld_blind_F_vs_mean.pdf")
meanSdPlot(assay(rld_blind_F))
dev.off()

#this one shows effect of transformation across all samples without accounting for experimental design
grDevices::cairo_pdf(file="stdev_rld_blind_T_vs_mean.pdf")
meanSdPlot(assay(rld_blind_T))
dev.off()

## ---- Sample distance heatmaps ------------------------------------------------

# Colours for plots below
##use RColorBrewer or any from https://github.com/EmilHvitfeldt/r-color-palettes
#library("RColorBrewer") or make a vector for color values wanted (can also make a gradient using https://mycolor.space/gradient3
# #ff0000, #ff5400, #ff7b00, #ff9b19, #ffb73a, #ffab3f, #ff9f44, #ff944a, #ff5658, #ff0084, #d300c2, #0432ff
mycols <- c("#ff0000", "#FFD966","#0432FF", "#3b3b3b", "#777777", "#b9b9b9")

#this one for sample distance calculated using rlog transformation using the experimental design (blind=FALSE) without limma to remove batch variation
library(gplots)

sampleDists <- as.matrix(dist(t(assay(rld_blind_F_no_limma))))
ncol(sampleDists)
grDevices::cairo_pdf(file="qc-heatmap-samples_rld_blind_F_no_limma.pdf", width=20, height=20, pointsize=25)
heatmap.2(as.matrix(sampleDists), key=T, key.xlab="Value", key.ylab="count", symkey=FALSE, keysize= 0.2, trace="none", key.par = list(cex=0.8), lhei=c(1.5,3), lwid=c(3,6),
          col=colorpanel(50, "black", "white"),
          margin=c(10, 10))
dev.off()

#this one for sample distance calculated using rlog transformation using the experimental design to remove batch variation (blind=FALSE)
sampleDists <- as.matrix(dist(t(assay(rld_blind_F))))
ncol(sampleDists)

grDevices::cairo_pdf(file="qc-heatmap-samples_rld_blind_F.pdf", width=20, height=20, pointsize=25)
heatmap.2(as.matrix(sampleDists), key=T, key.xlab="Value", key.ylab="count", symkey=FALSE, keysize= 0.2, trace="none", key.par = list(cex=0.8), lhei=c(1.5,3), lwid=c(3,6),
          col=colorpanel(50, "black", "white"),
          margin=c(10, 10))
dev.off()

#this one for sample distance calculated using rlog transformation without accounting for design (blind=TRUE)
sampleDists <- as.matrix(dist(t(assay(rld_blind_T))))
ncol(sampleDists)

grDevices::cairo_pdf(file="qc-heatmap-samples_rld_blind_T.pdf", width=20, height=20, pointsize=25)
heatmap.2(as.matrix(sampleDists), key=T, key.xlab="Value", key.ylab="count", symkey=FALSE, keysize= 0.2, trace="none", key.par = list(cex=0.8), lhei=c(1.5,3), lwid=c(3,6),
          col=colorpanel(50, "black", "white"),
          margin=c(10, 10))
dev.off()

## ---- Principal components analysis (PCA) ----------------------------------------

## used a method from Stephen Turner, @genetics_blog

#using rlog transformation using the experimental design to remove batch variation (blind=FALSE) without removing batch effect with limma

png(file="qc-pca_rld_blind_F_no_limma.png", 1500, 1500, pointsize=20)
rld_pca <- function (rld_blind_F_no_limma, intgroup = c("condition","time"), ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, repel=TRUE, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  require(DESeq2)
  rv = rowVars(assay(rld_blind_F_no_limma))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld_blind_F_no_limma)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld_blind_F_no_limma)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "#ff0000")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
rld_pca(rld_blind_F_no_limma, colors=mycols, intgroup=c("condition","time"), xlim=c(-50, 50))
dev.off()

#using rlog transformation using the experimental design to remove batch variation (blind=FALSE)

png(file="qc-pca_rld_blind_F.png", 1500, 1500, pointsize=20)
rld_pca <- function (rld_blind_F, intgroup = c("condition","time"), ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  require(DESeq2)
  rv = rowVars(assay(rld_blind_F))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld_blind_F)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld_blind_F)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "#ff0000")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
rld_pca(rld_blind_F, colors=mycols, intgroup=c("condition","time"), xlim=c(-50, 50))
dev.off()

#using rlog transformation without accounting for design (blind=TRUE)
png(file="qc-pca_rld_blind_T.png", 1500, 1500, pointsize=20)
rld_pca <- function (rld_blind_T, intgroup = c("condition","time"), ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  require(DESeq2)
  rv = rowVars(assay(rld_blind_T))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld_blind_T)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld_blind_T)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "#ff0000")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
rld_pca(rld_blind_T, colors=mycols, intgroup=c("condition","time"), xlim=c(-50, 50))
dev.off()

# RUN THE DESeq2 PIPELINE: GENERATE DE ANALYSIS RESULTS------------------------------------------------------------------

#recheck that colData has outlier samples removed
print(rownames(colData(dds)))

#perform DE analysis 
dds <- DESeq(dds)

## ---- Plot dispersions------------------------------------------------------------------
png(file="dds_qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, cex=1, main="Dispersion plot", xlab="mean of normalized counts", ylab="dispersion")
dev.off()

## ----Order DE results by increasing padj and save results tables------------------------------------------------------------------

dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/")
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/")
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/tables/")
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/tables/")

res_no_nic_T2_vs_no_nic_T1 <- results(dds, contrast=c("group","no_nic_T2","no_nic_T1"))
res_no_nic_T2_vs_no_nic_T1_Ordered <- res_no_nic_T2_vs_no_nic_T1[order(res_no_nic_T2_vs_no_nic_T1$padj),]
sig_res_no_nic_T2_vs_no_nic_T1_Ordered <- subset(res_no_nic_T2_vs_no_nic_T1_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_no_nic_T2_vs_no_nic_T1_Ordered), file="res_genenames_no_nic_T2_vs_no_nic_T1.csv")
write.csv(as.data.frame(sig_res_no_nic_T2_vs_no_nic_T1_Ordered), file="sig_res_no_nic_T2_vs_no_nic_T1_Ordered.csv")


res_no_nic_T3_vs_no_nic_T1 <- results(dds, contrast=c("group","no_nic_T3","no_nic_T1"))
res_no_nic_T3_vs_no_nic_T1_Ordered <- res_no_nic_T3_vs_no_nic_T1[order(res_no_nic_T3_vs_no_nic_T1$padj),]
sig_res_no_nic_T3_vs_no_nic_T1_Ordered <- subset(res_no_nic_T3_vs_no_nic_T1_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_no_nic_T3_vs_no_nic_T1_Ordered), file="res_genenames_no_nic_T3_vs_no_nic_T1.csv")
write.csv(as.data.frame(sig_res_no_nic_T3_vs_no_nic_T1_Ordered), file="sig_res_no_nic_T3_vs_no_nic_T1_Ordered.csv")


res_no_nic_T3_vs_no_nic_T2 <- results(dds, contrast=c("group","no_nic_T3","no_nic_T2"))
res_no_nic_T3_vs_no_nic_T2_Ordered <- res_no_nic_T3_vs_no_nic_T2[order(res_no_nic_T3_vs_no_nic_T2$padj),]
sig_res_no_nic_T3_vs_no_nic_T2_Ordered <- subset(res_no_nic_T3_vs_no_nic_T2_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_no_nic_T3_vs_no_nic_T2_Ordered), file="res_genenames_no_nic_T3_vs_no_nic_T2.csv")
write.csv(as.data.frame(sig_res_no_nic_T3_vs_no_nic_T2_Ordered), file="sig_res_no_nic_T3_vs_no_nic_T2_Ordered.csv")

res_T1_nic_vs_no_nic <- results(dds, contrast=c("group","nic_T1","no_nic_T1"))
res_T1_nic_vs_no_nic_Ordered <- res_T1_nic_vs_no_nic[order(res_T1_nic_vs_no_nic$padj),]
sig_res_T1_nic_vs_no_nic_Ordered <- subset(res_T1_nic_vs_no_nic_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_T1_nic_vs_no_nic_Ordered), file="res_genenames_T1_nic_vs_no_nic.csv")
write.csv(as.data.frame(sig_res_T1_nic_vs_no_nic_Ordered), file="sig_res_T1_nic_vs_no_nic_Ordered.csv")


res_T2_nic_vs_no_nic <- results(dds, contrast=c("group","nic_T2","no_nic_T2"))
res_T2_nic_vs_no_nic_Ordered <- res_T2_nic_vs_no_nic[order(res_T2_nic_vs_no_nic$padj),]
sig_res_T2_nic_vs_no_nic_Ordered <- subset(res_T2_nic_vs_no_nic_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_T2_nic_vs_no_nic_Ordered), file="res_genenames_T2_nic_vs_no_nic.csv")
write.csv(as.data.frame(sig_res_T2_nic_vs_no_nic_Ordered), file="sig_res_T2_nic_vs_no_nic_Ordered.csv")

res_T3_nic_vs_no_nic <- results(dds, contrast=c("group","nic_T3","no_nic_T3"))
res_T3_nic_vs_no_nic_Ordered <- res_T3_nic_vs_no_nic[order(res_T3_nic_vs_no_nic$padj),]
sig_res_T3_nic_vs_no_nic_Ordered <- subset(res_T3_nic_vs_no_nic_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_T3_nic_vs_no_nic_Ordered), file="res_genenames_T3_nic_vs_no_nic.csv")
write.csv(as.data.frame(sig_res_T3_nic_vs_no_nic_Ordered), file="sig_res_T3_nic_vs_no_nic_Ordered.csv")


res_nic_T2_vs_nic_T1 <- results(dds, contrast=c("group","nic_T2","nic_T1"))
res_nic_T2_vs_nic_T1_Ordered <- res_nic_T2_vs_nic_T1[order(res_nic_T2_vs_nic_T1$padj),]
sig_res_nic_T2_vs_nic_T1_Ordered <- subset(res_nic_T2_vs_nic_T1_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_nic_T2_vs_nic_T1_Ordered), file="res_genenames_nic_T2_vs_nic_T1.csv")
write.csv(as.data.frame(sig_res_nic_T2_vs_nic_T1_Ordered), file="sig_res_nic_T2_vs_nic_T1_Ordered.csv")


res_nic_T3_vs_nic_T1 <- results(dds, contrast=c("group","nic_T3","nic_T1"))
res_nic_T3_vs_nic_T1_Ordered <- res_nic_T3_vs_nic_T1[order(res_nic_T3_vs_nic_T1$padj),]
sig_res_nic_T3_vs_nic_T1_Ordered <- subset(res_nic_T3_vs_nic_T1_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_nic_T3_vs_nic_T1_Ordered), file="res_genenames_nic_T3_vs_nic_T1.csv")
write.csv(as.data.frame(sig_res_nic_T3_vs_nic_T1_Ordered), file="sig_res_nic_T3_vs_nic_T1_Ordered.csv")


res_nic_T3_vs_nic_T2 <- results(dds, contrast=c("group","nic_T3","nic_T2"))
res_nic_T3_vs_nic_T2_Ordered <- res_nic_T3_vs_nic_T2[order(res_nic_T3_vs_nic_T2$padj),]
sig_res_nic_T3_vs_nic_T2_Ordered <- subset(res_nic_T3_vs_nic_T2_Ordered, padj<0.1 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_nic_T3_vs_nic_T2_Ordered), file="res_genenames_nic_T3_vs_nic_T2.csv")
write.csv(as.data.frame(sig_res_nic_T3_vs_nic_T2_Ordered), file="sig_res_nic_T3_vs_nic_T2_Ordered.csv")

## ----Volcano plots without and with labels ---------------------------------------------------

setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/")
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/figures/")

#SAVE Volcano plots with unlabeled genes
dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/figures/unlabeled")
setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/figures/unlabeled")


#SAVE Volcano plots with "significant" genes labeled
# dir.create("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/figures/labeled")
# setwd("C:/Users/Amada/Desktop/Pnic_multiomic_study_nicotine/R_data_processing/kept_samples_DE_results/with_genenames/figures/labeled")

#generate volcano plot for each comparison of interest
volcanoplot <- function (res_nic_T3_vs_nic_T2_Ordered, lfcthresh=1, sigthresh=0.1, main=NULL, legendpos="top", labelsig=TRUE, textcx=1, ...) {
  with(res_nic_T3_vs_nic_T2_Ordered, plot(log2FoldChange, -log10(pvalue), main=main, pch=20, col="gray", cex=1.2, ...))
    # Not Significant because padj > sigthresh
  with(subset(res_nic_T3_vs_nic_T2_Ordered, padj>sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="#4D4D4D", bg="#4D4D4D", cex=1.2,...))
    # Not Significant because LFC < sigthresh 
  with(subset(res_nic_T3_vs_nic_T2_Ordered, padj<sigthresh & abs(log2FoldChange)<lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="black", bg="black", cex=1.2, ...))
    # significantly up-regulated: padj < 0.1, LFC > 1
  with(subset(res_nic_T3_vs_nic_T2_Ordered, (padj<sigthresh & log2FoldChange>lfcthresh)), points(log2FoldChange, -log10(pvalue), pch=24, bg="red", col="red", cex=1.5, ...))
    # significantly down-regulated: padj < 0.1, LFC < -1
  with(subset(res_nic_T3_vs_nic_T2_Ordered, (padj<sigthresh & log2FoldChange<(-lfcthresh))), points(log2FoldChange, -log10(pvalue), pch=25, bg="green", col="green", cex=1.5, ...))
  
    # add gene labels
  if (labelsig) 
  {
    require(calibrate)
    # make a vector of significant rownames and italicise them
    gene_names_italics <- lapply(rownames(subset(res_nic_T3_vs_nic_T2_Ordered, padj<sigthresh & abs(log2FoldChange)>lfcthresh)), function(x) bquote(italic(.(x))))
    
    # #add gene labels to significant genes
    # with(subset(res_nic_T3_vs_nic_T2_Ordered, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=as.expression(gene_names_italics), cex=1, ...))
  }
  
  # add horizontal line - p-adj thresh
  abline(h = -log10(sigthresh), col = "black", lty=2, lwd=4)
  text(x=-7, y=3, 'p-value < 0.1', cex=1.2)
  # text(x=-7, y=1.08, 'p-value < 0.1', cex=1.2)
  # text(x=-7, y=3, 'p-value < 0.1', cex=1.2)
  
  # add vertical lines - LFC thresh
  abline(v = lfcthresh, col = "black", lty=3, lwd=4)
  text(x=-1.2, y=30, srt=90, 'FC < -2', cex =1.2)
  # text(x=-1.2, y=1.98, srt=90, 'FC < -2', cex =1.2)
  # text(x=-1.2, y=8, srt=90, 'FC < -2', cex =1.2)
  
  abline(v = -lfcthresh, col = "black", lty=3, lwd=4)
  text(x=1.2, y=30, srt=90, 'FC > 2', cex=1.2)
  # text(x=1.2, y=1.98, srt=90, 'FC > 2', cex=1.2)
  # text(x=1.2, y=8, srt=90, 'FC > 2', cex=1.2)
  
  # add legend
    legend(legendpos, cex=1.4, text.width=c(strwidth("Up-regulated          "), strwidth("Down-regulated           "), strwidth("p-adj < 0.1        "), 
                                          strwidth("-2 < FC > 2        "), strwidth("Not significant")), 
         bty='n', yjust=0.5, xpd=TRUE, horiz=TRUE, inset=-0.15, 
         legend=c(paste("Up-regulated"), paste("Down-regulated"), paste("p-adj < 0.1"), paste("-2 < FC > 2"), paste("Not significant")), 
         pch=c(24,25,20,20,20), col=c("red","green","black","#4D4D4D","gray"), pt.bg=c("red","green","black","#4D4D4D","gray"))
}
png(file="volcanoplot_nic_T3_vs_nic_T2_no_labs.png", 5000, 3000, pointsize=40, res=140)
volcanoplot(res_nic_T3_vs_nic_T2_Ordered, mgp=c(2.5,0.8,0), cex.axis=1.5, cex.lab=1.5, lfcthresh=1, sigthresh=0.1, textcx=1, xlim=c(-8, 8), xlab="log2FoldChange", ylab="-log10(p-value)")
dev.off()
