# library loading
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

# setting working directory
setwd("/Users/alanx/OneDrive/Desktop/QBIO490/qbio_490_alan/analysis_data")
# data loading
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
clinic <- read.csv("/Users/alanx/OneDrive/Desktop/QBIO490/qbio_490_alan/analysis_data/brca_clinical_data.csv")
clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug",
                                   directory = "/Users/alanx/OneDrive/Documents/GDCdata")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation",
                                  directory = "/Users/alanx/OneDrive/Documents/GDCdata")

sum(is.na(clinic$age_at_initial_pathologic_diagnosis))
# The variable chosen for clinic data is age_at_initial_pathologic_diagnosis
# This variable is a discrete variable
# The variable describes how old the patient is when first received diagnosis
clinical_rad$radiation_dosage
# The variable chosen for clinical_rad data is radiation_dosage
# The variable is continuous
# The variable describes how much radiation the patient has received through treatment.

# The younger the patient is at initial diagnosis, the higher the radiation dosage.
# The younger the patient is at initial diagnosis, the higher the survivability.
# The higher the radiation dosage, the higher the survivability.

# getting only the first occurrence of bcr_patient_barcode in clinical_rad
unique_rad_bcr <- clinical_rad[match(unique(clinical_rad$bcr_patient_barcode), 
                                     clinical_rad$bcr_patient_barcode), ]
# first subset the rows that have radiation information
clinic_has_rad_info <- subset(clinic, clinic$has_radiations_information == "YES")
# only extract the first occurrence of a particular bcr barcode
unique_clinic_has_rad <- clinic_has_rad_info[match(unique(clinic_has_rad_info$bcr_patient_barcode), 
                                                    clinic_has_rad_info$bcr_patient_barcode),]

plot(unique_clinic_has_rad$age_at_initial_pathologic_diagnosis, 
     unique_rad_bcr$radiation_dosage, xlab = "Age at Which Patient Was First Diagnosed", 
     ylab = "Radiation Dosage Patient Received", 
     main = "Radiation Dosage Patient Received vs. Age at Which Patient Was First Diagnosed", 
     col = 4, pch = 4)
# A scatter plot is chosen since both variables are of numeric nature and are not
# categorical. There is not a clear relationship between radiation dosage and age
# of initial diagnosis
#####################################################################################
# KM Plot
# remove NA
age_na <- ifelse(is.na(clinic$age_at_initial_pathologic_diagnosis), F, T)
clinic_no_na <- clinic[age_na, ]
# determine what it means to be young, middle and old
quantile(clinic_no_na$age_at_initial_pathologic_diagnosis)
# < 49 young, < 68 middle, > 68 old
young_mask <- ifelse(clinic_no_na$age_at_initial_pathologic_diagnosis <= 49, T, F)
middle_mask <- ifelse(clinic_no_na$age_at_initial_pathologic_diagnosis > 49 & 
                        clinic_no_na$age_at_initial_pathologic_diagnosis <= 68,
                        T, F)
old_mask <- ifelse(clinic_no_na$age_at_initial_pathologic_diagnosis > 68, T, F)
clinic_no_na$age_status <- ifelse(young_mask, "Young", ifelse(middle_mask, "Middle", "Old"))

# making a survival time column
clinic_no_na$survival_time <- ifelse(is.na(clinic_no_na$days_to_death),
                                     clinic_no_na$survival_time <- clinic_no_na$days_to_last_followup,
                                     clinic_no_na$survival_time <- clinic_no_na$days_to_death)
# remove all -Inf in survival time
inf_mask <- ifelse(clinic_no_na$survival_time == -Inf, F,T)
clinic_no_na <- clinic_no_na[inf_mask,]

# adding death event
clinic_no_na$death_event <- ifelse(clinic_no_na$vital_status == "Alive",
                                   clinic_no_na$death_event <- FALSE,
                                   clinic_no_na$death_event <- TRUE)
# Making the Plot
survival_object <- Surv(time = clinic_no_na$survival_time,
                        event = clinic_no_na$death_event)
fit_object <- survfit(survival_object ~ clinic_no_na$age_status,
                      data = clinic_no_na)
survplot <- ggsurvplot(fit_object, pval= TRUE, ggtheme = 
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))
KM_plot
ggsave("KM_plot_Age.jpeg", KM_plot, device = "jpeg", 
       path = "/Users/alanx/OneDrive/Desktop/QBIO490/qbio_490_alan/analysis_data")

#####################################################################################
# Second KM Plot
# subset out all the blank data
rad_na_mask <- ifelse(unique_rad_bcr$radiation_dosage == "", F, T)
rad_no_na <- unique_rad_bcr[rad_na_mask, ]

# Create Columns for high, middle and low dosage
quantile(as.numeric(rad_no_na$radiation_dosage))
# low dosage <= 78
# medium dosage <= 102
# high dosage > 102
low_mask <- ifelse(as.numeric(rad_no_na$radiation_dosage) <= 78, T, F)
medium_mask <- ifelse(as.numeric(rad_no_na$radiation_dosage) > 78 & 
                        as.numeric(rad_no_na$radiation_dosage) <= 102,
                      T, F)
high_mask <- ifelse(as.numeric(rad_no_na$radiation_dosage) > 102, T, F)
rad_no_na$radiation_status <- ifelse(low_mask, "Low", ifelse(medium_mask, "Medium", "High"))
# merge info from clinic to rad_no_na by barcode
rad_no_na <- merge(rad_no_na, clinic, by.x = "bcr_patient_barcode",
                   by.y = "bcr_patient_barcode")
# making survival time column
rad_no_na$survival_time <- ifelse(is.na(rad_no_na$days_to_death),
                                     rad_no_na$survival_time <- rad_no_na$days_to_last_followup,
                                     rad_no_na$survival_time <- rad_no_na$days_to_death)

# remove all -Inf in survival time
inf_mask <- ifelse(rad_no_na$survival_time == -Inf, F,T)
rad_no_na <- rad_no_na[inf_mask,]

# adding death event
rad_no_na$death_event <- ifelse(rad_no_na$vital_status == "Alive",
                                rad_no_na$death_event <- FALSE,
                                rad_no_na$death_event <- TRUE)
# Making the Plot
survival_object <- Surv(time = rad_no_na$survival_time,
                        event = rad_no_na$death_event)
fit_object <- survfit(survival_object ~ rad_no_na$radiation_status,
                      data = rad_no_na)
survplot <- ggsurvplot(fit_object, pval= TRUE, ggtheme = 
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")
KM_plot_dosage <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))
KM_plot_dosage
ggsave("KM_plot_dosage.jpeg", KM_plot_dosage, device = "jpeg", 
       path = "/Users/alanx/OneDrive/Desktop/QBIO490/qbio_490_alan/analysis_data")

# KM plot for diagnosis age suggests that lower the initial diagnosis age, the more
# survivable the cancer is.
# KM plot for radiation dosage does not paint a clear picture as whether there is
# a correlation between radiation dosage and survivability of cancer.
# The p-value for the KM plot for diagnosis age is 0.0013.
# The p-value for the KM plot for radiation dosage is 0.18.
# With a conventional p-value cut off of 0.05, it shows initial diagnosis age is
# significant for the survivability of cancer and the radiation dosage is not
# significant for the survivability of cancer.