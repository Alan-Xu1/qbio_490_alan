# Author: Alan Xu
# Intro to Clinical Data: Part II
# Set Working Directory
setwd("/Users/alanx/OneDrive/Desktop/QBIO490/qbio_490_alan/analysis_data")

# Library Loading
if (!require(BiocManager)){
  install.packges("BiocManager")
}
library(BiocManager)
if (!require(TCGAbiolinks)){
  install.packges("TCGAbiolinks")
}
library(TCGAbiolinks)
if (!require(maftools)){
  install.packges("maftools")
}
library(maftools)
if (!require(ggplot2)){
  install.packges("ggplot2")
}
library(ggplot2)

# Reading in Clinical Data and MAF File
clinical <- read.csv("/Users/alanx/OneDrive/Desktop/QBIO490/qbio_490_alan/analysis_data/brca_clinical_data.csv")
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
maf <- GDCprepare(maf_query, directory = "/Users/alanx/OneDrive/Documents/GDCdata")
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

# Start of Data Analysis
# Only take the MAF info of patients who have radiation_therapy information
has_rad_mask <- ifelse(maf_object@clinical.data$radiation_therapy != "",
                         TRUE,
                         FALSE)
has_rad_barcode <- maf_object@clinical.data[has_rad_mask,Tumor_Sample_Barcode]
maf_has_rad <- subsetMaf(maf = maf_object,
                           tsb = has_rad_barcode)
maf_has_rad@clinical.data$radiation_therapy <-
  as.factor(maf_has_rad@clinical.data$radiation_therapy)
# Split the has radiation_therapy info group into 2 categories
# Those who underwent radiation therapy and those who didn't
rad_yes_mask <- ifelse(maf_has_rad@clinical.data$radiation_therapy == "YES",
                      TRUE,
                      FALSE)
rad_yes_barcode <- maf_has_rad@clinical.data[rad_yes_mask, Tumor_Sample_Barcode]
maf_yes_rad <- subsetMaf(maf = maf_has_rad, tsb = rad_yes_barcode)
rad_no_mask <- ifelse(maf_has_rad@clinical.data$radiation_therapy == "NO",
                       TRUE,
                       FALSE)
rad_no_barcode <- maf_has_rad@clinical.data[rad_no_mask, Tumor_Sample_Barcode]
maf_no_rad <- subsetMaf(maf = maf_has_rad, tsb = rad_no_barcode)
# Determine top 10 most mutated genes in the group
oncoplot(maf = maf_yes_rad, top = 10)
# Creating Co-Oncoplot
coOncoplot(m1 = maf_yes_rad, 
           m2 = maf_no_rad,
           genes = c("PIK3CA", "TP53", "CDH1", "TTN", "RYR2", "KMT2C", "MAP3K1",
                     "TBX3", "ANKFN1", "SYNE1"),
           m1Name = "Patients that underwent radiation therapy", 
           m2Name = "Patients that did not undergo radiation therapy", 
           borderCol = NA)
# PIK3CA is significantly more frequently mutated in patients that underwent 
# radiation therapy. 
# PIK3CA overactivates PI3K enzyme which allows for cancer proliferation.
# In breast cancer, PIK3CA mutational status is correlated with response to 
# therapy where mutations can confer resistance to some antibody-based
# therapy.
# See: Samuels and Walkman, "Oncogenic Mutations of PIK3CA in Human Cancers"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3164550/
# See: https://www.cancer.gov/publications/dictionaries/cancer-terms/def/pik3ca-gene
# Since mutations in PIK3CA confer cancer proliferation, more aggresive
# therapy such as radiation would be warranted in patients with mutations in
# PIK3CA.

# Fisher's Exact Test
# getting the number of patient that underwent radiation therapy
barcode_rad_yes <- maf_yes_rad@clinical.data$Tumor_Sample_Barcode
len_rad_yes <- length(barcode_rad_yes)
# getting the number of patient that has mutation in PIK3CA
genePIK_maf <- subsetMaf(maf = maf_has_rad,
                       genes = "PIK3CA")
barcode_PIK_yes <- genePIK_maf@clinical.data$Tumor_Sample_Barcode
len_PIK_yes <- length(barcode_PIK_yes)
# intersect two population to get patients that have PIK3CA mutation and
# underwent radiation therapy
both_PIK_and_Rad <- intersect(barcode_rad_yes,barcode_PIK_yes)
len_PIK_and_Rad <- length(both_PIK_and_Rad)
# getting number of patient that only underwent radiation therapy or only have
# mutation in PIK3CA
len_rad_only <- len_rad_yes - len_PIK_and_Rad
len_PIK_only <- len_PIK_yes - len_PIK_and_Rad
# getting number of patient that do not have mutation in PIK3CA nor went through
# radiation therapy
len_no_PIK_rad <- length(maf_has_rad@clinical.data$Tumor_Sample_Barcode) - len_PIK_and_Rad - len_PIK_only - len_rad_only
# constructing the contingency table
contig <- matrix(c(len_PIK_and_Rad, 
                   len_rad_only,
                   len_PIK_only,
                   len_no_PIK_rad), 
                   nrow=2)
contig
# making the mosaic plot
mosaicplot(contig)
# conducting fisher exact test
fisher_test <- fisher.test(contig)
fisher_test

# Using the convention of having 0.05 as the p-value cutoff, there is no 
# statistically significant correlation between PIK3CA mutation status
# and whether the patient goes through radiation therapy.
# The odds ration of 1.6 indicates that PIK3CA mutation and whether patient goes
# through radiation therapy is co-occurent.

# co-lollipop
lollipopPlot2(m1 = maf_yes_rad, 
              m2 = maf_no_rad, 
              m1_name = "Patients that underwent radiation treatment",
              m2_name = "Patients that did not receieve radiation treatment",
              gene = "PIK3CA")
# In PIK3CA gene, patients who did not undergo radiation treatment tend to have
# missense mutation in PI3K_classI_alpha region.

# KM plot
# making the survival status a factor
maf_object@clinical.data$Overall_Survival_Status <-
  ifelse(maf_object@clinical.data$vital_status == "Alive",
         maf_object@clinical.data$Overall_Survival_Status <- TRUE,
         maf_object@clinical.data$Overall_Survival_Status <- FALSE)
# constructing the KM plot
mafSurvival(maf = maf_object,
            genes = "PIK3CA", ## pick a gene of your choosing
            time = "days_to_last_followup", ## name of the column in maf_object@clinical.data containing survival time
            Status = "Overall_Survival_Status", ## name of the column that contains a boolean value for death events, you may need to recreate this... 
            isTCGA = TRUE)
# Number of patients that are available for analysis
length(maf_object@clinical.data$Tumor_Sample_Barcode)
# Using the conventional p-value cutoff of 0.05, the mutation status in PIK3CA is not
# significant against the survivability of cancer. From the analysis and literature
# research done above,
# it does seem like the mutation status in PIK3CA is correlated with more
# severe form of cancer. One possibility is that the population size is too
# small to detect a difference in survivability based on the mutational status
# of PIK3CA, so it could be necessary to expand the size of the study.