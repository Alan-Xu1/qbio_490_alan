# Author: Alan Xu
# Assignment: Midterm Project
# Class: QBIO 490 Directed Research

# Set Working Directory
setwd("/Users/alanx/OneDrive/Desktop/QBIO490/qbio_490_alan/midterm_project_Xu")
dir.create("outputs")
setwd("outputs")

# Library Loading
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library(BiocManager)
if (!require("TCGAbiolinks", quietly = TRUE)){
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)
if (!require("maftools", quietly = TRUE)){
  BiocManager::install("maftools")
}
library(maftools)
if (!require("maftools", quietly = TRUE)){
  BiocManager::install("DESeq2")
}
library(DESeq2)
if (!require(dplyr)){
  install.packges("dplyr")
}
library(dplyr)
if (!require(survival)){
  install.packges("survival")
}
library(survival)
if (!require(survminer)){
  install.packges("survminer")
}
library(survminer)
if (!require(ggplot2)){
  install.packges("ggplot2")
}
library(ggplot2)
if (!require(EnhancedVolcano)){
  BiocManager::install('EnhancedVolcano')
}
library(EnhancedVolcano)

# Data Loading
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")

colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <-
  "Tumor_Sample_Barcode"
maf_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf,
                       clinicalData = clinic,
                       isTCGA = TRUE)

rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

### Start of Data Analysis
# The objective of this analysis is to see if there is any mutations associated
# with HER2 Immunohistocehmistry level
# First we want to identify if there is a correlation between HER2 immunohistochemistry
# level and survivability of cancer
# Remove all NA HER2
clinic_her2 <- subset(clinic, clinic$her2_immunohistochemistry_level_result != ""
                      & is.na(clinic$her2_immunohistochemistry_level_result) == FALSE)
# making a survival time column
clinic_her2$survival_time <- ifelse(is.na(clinic_her2$days_to_death) ,
                                    clinic_her2$survival_time <- clinic_her2$days_to_last_followup,
                                    clinic_her2$survival_time <- clinic_her2$days_to_death)
# remove all -Inf in survival time
inf_mask <- ifelse(clinic_her2$survival_time == -Inf, F,T)
clinic_her2 <- clinic_her2[inf_mask,]

# adding death event
clinic_her2$death_event <- ifelse(clinic_her2$vital_status == "Alive",
                                  clinic_her2$death_event <- FALSE,
                                  clinic_her2$death_event <- TRUE)
# Plotting a bar plot of patients based on their her2 expression level
ggplot(clinic_her2, aes(x = her2_immunohistochemistry_level_result)) +
  geom_bar() + ggtitle("Number of patients by HER2 Immunhistochemistry Level") +
  xlab("Semi-Quantitative Immunohistochemistry Level") + ylab("Number of Patients")
ggsave("barplot.jpeg")
# Making KM Plot
survival_object <- Surv(time = clinic_her2$survival_time,
                        event = clinic_her2$death_event)
fit_object <- survfit(survival_object ~ clinic_her2$her2_immunohistochemistry_level_result,
                      data = clinic_her2)
survplot <- ggsurvplot(fit_object, pval= TRUE, ggtheme = 
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))
KM_plot
ggsave("KM_Plot.jpeg")
### MAF Object Processing
# Observe if there is different gene mutations in these HER2_low (0,1+) and
# HER2_high (2+,3+) groups

# Extract the MAFs that contain HER2 Information
has_her2_mask <- ifelse(maf_object@clinical.data$her2_immunohistochemistry_level_result != "" 
                        & is.na(maf_object@clinical.data$her2_immunohistochemistry_level_result) == FALSE,
                        TRUE,
                        FALSE)
has_her2_barcode <- maf_object@clinical.data[has_her2_mask,Tumor_Sample_Barcode]
maf_her2 <- subsetMaf(maf = maf_object, tsb = has_her2_barcode)
# Split the MAFs that have HER2 Information Into the Low and High HER2 Group
# Low HER2 Group
her2_low_mask <- ifelse(maf_her2@clinical.data$her2_immunohistochemistry_level_result == "0"
                        | maf_her2@clinical.data$her2_immunohistochemistry_level_result == "1+",
                        TRUE,
                        FALSE)
her2_low_barcode <- maf_her2@clinical.data[her2_low_mask, Tumor_Sample_Barcode]
maf_low_her2 <- subsetMaf(maf = maf_her2, tsb = her2_low_barcode)

# High HER2 Group
her2_high_mask <- ifelse(maf_her2@clinical.data$her2_immunohistochemistry_level_result == "2+"
                         | maf_her2@clinical.data$her2_immunohistochemistry_level_result == "3+",
                         TRUE,
                         FALSE)
her2_high_barcode <- maf_her2@clinical.data[her2_high_mask, Tumor_Sample_Barcode]
maf_high_her2 <- subsetMaf(maf = maf_her2, tsb = her2_high_barcode)
# Making the co-oncoplot
coOncoplot(m1 = maf_low_her2,
           m2 = maf_high_her2,
           m1Name = "Patients that express low level HER2 (0,1+)",
           m2Name = "Patients that express high level HER2 (2+,3+)",
           borderCol = NA)
ggsave("co-oncoplot.jpeg")
# PIK3CA and TP53 appears mutation status appears to be different for the HER2 low
# and high group
# See if the genomic mutational difference carries to transcriptional data
### Summarized Experiment (RNA-Seq) Processing

# Data Processing
# Getting the clinical information, gene information and count matrix
rna_clinical <- as.data.frame(rna_se@colData)
rna_genes <- as.data.frame(rna_se@rowRanges@elementMetadata)
rna_counts <- rna_se@assays@data$unstranded

# Getting only the data from cancerous tissues
cancer_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", FALSE, TRUE)
rna_clinical <- rna_clinical[cancer_mask,]
rna_counts <- rna_counts[, cancer_mask]

# We will control the tissue_or_organ_of_origin, ajcc_pathologic_stage and age_category for
# differential gene analysis
rna_clinical$age_category <- ifelse(rna_clinical$age_at_index > 58, "Old", "Young")
rna_clinical$age_category <- factor(rna_clinical$age_category)
rna_clinical$ajcc_pathologic_stage <- as.factor(rna_clinical$ajcc_pathologic_stage)
rna_clinical$tissue_or_organ_of_origin <- as.factor(rna_clinical$tissue_or_organ_of_origin)

# There appears to be some NAs in the info that require masking
sum(is.na(rna_clinical$age_category))
sum(is.na(rna_clinical$ajcc_pathologic_stage))
sum(is.na(rna_clinical$tissue_or_organ_of_origin))

# Mask out the NAs in both rna_clinical and rna_counts
na_mask <-  ifelse(is.na(rna_clinical$age_category), FALSE, ifelse(is.na(rna_clinical$ajcc_pathologic_stage), FALSE, ifelse(is.na(rna_clinical$tissue_or_organ_of_origin), FALSE, TRUE)))
rna_clinical <-  rna_clinical[na_mask,]
rna_counts <- rna_counts[,na_mask]
# We will also filter genes that have < 10 counts
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums < 10, FALSE, TRUE)
rna_counts <- rna_counts[low_counts_mask,]
rna_genes <- rna_genes[low_counts_mask,]
# Adding the patient barcode and Ensembl gene ID to the count matrix
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id

# Running DESeq2
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~tissue_or_organ_of_origin + ajcc_pathologic_stage + age_category)
dds_obj <- DESeq(dds)
results <- results(dds_obj, format = "DataFrame", contrast = c("age_category", "Young", "Old")) 

# Adding a gene name column for easy recognition
results$gene_name <- rna_genes$gene_name
# Adjusted p value cut-off at 0.05
EnhancedVolcano(results, x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                lab = results$gene_name)
ggsave("volcano.jpeg")