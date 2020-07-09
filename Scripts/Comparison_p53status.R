# Script to compare TP53 hits in pancancer analysis after filtering out
# cell-lines based on TP53 status using BayesDel score.

library(slidr)
library(dplyr)
setwd("~/Downloads/Slidr_Results_new/")

# Load the pancancer data
load("./PanCan8pc/pancan.Rdata")
rm(list=setdiff(ls(), "pc_data"))

# Hits without stratifying
hits_pancan <- read.delim("./PanCan8pc/Hit_List/SL_hits_pan_cancer.txt", stringsAsFactors = FALSE)

# Reading in the list based on TP53 status
tp53_annot <- read.delim("./Tp53classes/tp3_status.csv", sep = ",", stringsAsFactors = FALSE)

# Chose BayesDel to stratify since the cell-lines with tolerated mutations are a superset of
# the ones from Polyphen2 = B (except "kyse410") and REVEL < 0.5.
# More details on https://p53.iarc.fr/Manual.aspx#ProtDescription
tol_celllines <- tp53_annot %>%
                    dplyr::filter(BayesDel < 0.16) %>%
                    dplyr::pull(cellline)

# Remove these cell-lines from the mutations
filt_data <- pc_data
celllines <- colnames(filt_data$mutations)
filt_celllines <-  celllines[-which(grepl(paste(tol_celllines, collapse="|"), celllines))]

# Data object after filtering out the cell-lines with tolerated mutations in  TP53
filt_data$mutations     <- filt_data$mutations[, filt_celllines, drop = FALSE]
filt_data$viabilities   <- filt_data$viabilities[, filt_celllines]
filt_data$CNalterations <- filt_data$CNalterations[, filt_celllines]

# Run Slidr on the filtered object
path_results <- "./Tp53classes/"

hits_filt <- slidr::identifySLHits(canc_data = filt_data, path_results = path_results,
                                   WT_pval_thresh = 0, fp_thresh = 1)

# Compare the TP53 hits before and after filtering based on status
tp53_unfilt <- hits_pancan %>% dplyr::filter(driver_gene == "TP53")
tp53_filt   <- hits_filt %>% dplyr::filter(driver_gene == "TP53")

all(paste(tp53_unfilt[,1],tp53_unfilt[,2], sep = "_") %in%
    paste(tp53_filt[,1],tp53_filt[,2], sep = "_"))
