# This script runs PAM50 package that predicts PAM50 subtypes.
# Project page: https://ccchang0111.github.io/PAM50/
# It can predict 5 subtypes for breast cancer patients: Luminal A, Luminal B, HER2, Basal, and Normal. 
# Patients with Luminal A subtype usually have the best prognosis.

# load packages
library(readr)
library(tidyverse)
library(dplyr) # for Error: > batch1.pred <- PAM50(batch1_cpm_df, cutoff = 0)
# Error: unable to find an inherited method for function 'select' for signature 'x = "tbl_df"'
##### re-install corrected PAM50
#devtools::install_github("ccchang0111/PAM50")
library(PAM50)
library(Biobase) ## for parsing ExpressionSet data
library(biomaRt) ## fo converting Ensembl ID to Entrez ID, PAM50 needs Entrez IDs for prediction
library(edgeR) ## for generating CPM for raw readcounts data 
# CPM normalizes for sequencing depth, Comparing gene expression across samples

#' PAM50 main function
#' @description
#' This is what you need for getting PAM50 prediction for each sample. It takes
#' RNA-seq values (raw or CPM) and returns PAM50 subtype for each sample.

############################################################
############################################################

# Step 1: Input raw read-counts from batch 1 and 2
batch1_raw <- read.delim(here::here("data","rawRNA","batch1_combined_counts.txt"), comment.char = "#", check.names = FALSE)
batch2_raw <- read.delim(here::here("data","rawRNA","batch2_combined_counts.txt"), comment.char = "#", check.names = FALSE)

# clean column names for both batchs
# Batch 1:
## Keep only gene name and count columns in the raw read counts table
batch1_raw <- batch1_raw[, c(1, 7:ncol(batch1_raw))]
rownames(batch1_raw) <- batch1_raw$Geneid
## Clean column names/sample names (remove file paths and extensions)
# e.g. turn `../WT_1Aligned.sortedByCoord.out.bam` into `WT_1`
colnames(batch1_raw) <- gsub(".*/|Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(batch1_raw))

## modify sample names to "bng_bx_1" format
colnames(batch1_raw)
# Extract "bng_bx_*" or "mal_bx_*" from each sample column
modified_sample_names <- sub(".*_(bng_bx_\\w+|mal_bx_\\w+)_.*", "\\1", colnames(batch1_raw)[-1])
colnames(batch1_raw) <- c(colnames(batch1_raw)[1], modified_sample_names)
print("Batch 1 colnames:")
colnames(batch1_raw)

# Batch 2:
## Keep only gene name and count columns in the raw read counts table
batch2_raw <- batch2_raw[, c(1, 7:ncol(batch2_raw))]
rownames(batch2_raw) <- batch2_raw$Geneid
## Clean column names/sample names (remove file paths and extensions)
# e.g. turn `../WT_1Aligned.sortedByCoord.out.bam` into `WT_1`
colnames(batch2_raw) <- gsub(".*/|.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(batch2_raw))
## modify sample names to "SAM_1" format
colnames(batch2_raw)
## map sample names
# Load the mapping table
map_df <- read.delim(here::here("data","rawRNA","batch2_filenames_SampleIDs.txt"), sep = "\t", header = TRUE)
# Extract the relevant mapping from 'Index Name' + 'Library ID' → 'Description'
map_df$FullID <- paste0(sub("SM_*","",map_df$Index.Name), "_", map_df$Library.ID)
map_df$SAM <- sub(".*: ", "", map_df$Description)
# Remove rows where mapping isn't present
map_df <- map_df[!is.na(map_df$SAM), ]
# Get current column names of the count matrix (excluding "Geneid")
orig_cols <- colnames(batch2_raw)[-1]
rename_cols <- map_df$SAM[match(orig_cols, map_df$FullID)] 
print("rename_cols:")
rename_cols
# Apply the new names
colnames(batch2_raw)[-1] <- rename_cols
print("Batch 2 colnames:")
colnames(batch2_raw)


# Step 2: CPM conversion using cpm() function in the `edgeR` package
# log = TRUE: applies log2 transfromation to the CPM values
# proir.count = 1: adds 1 before log2 to avoid log1(0)
print("Batch 1 raw readcounts Number of genes:")
length(batch1_raw[[1]]) 
#63086
batch1_mat <- as.data.frame(batch1_raw[,-1]) # convert raw read counts to count matrix
rownames(batch1_mat) <- batch1_raw[[1]] # set Ensembl IDs to be rownames
batch1_cpm <- cpm(batch1_mat, log = TRUE, prior.count = 1) 

print("Batch 2 raw readcounts Number of genes:")
length(batch2_raw[[1]])
batch2_mat <- as.data.frame(batch2_raw[,-1])
rownames(batch2_mat) <- batch2_raw[[1]]
batch2_cpm <- cpm(batch2_mat, log = TRUE, prior.count = 1) 

length(rownames(batch1_cpm))
length(rownames(batch2_cpm))

# Step 5: Run PAM50 prediction
# cutoff: only the subtype with > cutoff probability and is the max among all subtypes will be assigned the final PAM50 subtype.
# if none of the subtypes produce peobability > cutoff, the final subtyoe = NA

# check for duplicated Entrez IDs
duplicated_ids <- rownames(batch1_cpm)[duplicated(rownames(batch1_cpm))]
table(duplicated_ids)

# Remove duplicated Entrez IDs
batch1_cpm_df <- as.data.frame(batch1_cpm[!duplicated(rownames(batch1_cpm)), ])

# Error in PAM50(): Error: unable to find an inherited method for function 'select' for signature 'x = "tbl_df"'
# Solved by issue: https://github.com/ccchang0111/PAM50/issues/1
# Before running PAM50, run: fix(PAM50) and change select(-SampleID) into dplyr::select(-SmapleID)

# PAM50 predict
batch1.pred <- PAM50(batch1_cpm_df, cutoff = 0)
table(batch1.pred[[1]])
# LumA Normal
# 15     77
batch1.pred$Sample <- colnames(batch1_cpm_df)
write_csv(batch1.pred, here::here("results","PAM50", "batch1_PAM50_predictions.csv"))
print("Batch1 prediction done!")

### Batch 2
# check for duplicated Entrez IDs
duplicated_ids <- rownames(batch2_cpm)[duplicated(rownames(batch2_cpm))]
table(duplicated_ids)

print("batch2_cpm:")
head(batch2_cpm)

# Remove duplicated Entrez IDs and predict
batch2_cpm_df <- as.data.frame(batch2_cpm[!duplicated(rownames(batch2_cpm)), ])
batch2.pred <- PAM50(batch2_cpm_df, cutoff = 0)
table(batch2.pred[[1]])
batch2.pred$Sample <- colnames(batch2_cpm_df)
write_csv(batch2.pred, here::here("results","PAM50", "batch2_PAM50_predictions.csv"))

