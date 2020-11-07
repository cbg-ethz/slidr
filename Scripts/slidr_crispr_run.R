# Crispr analysis

library(slidr)
library(dplyr)
setwd("~/Downloads/Slidr_Results_new/")

# Load pan cancer object
load("./PanCan8pc/pancan.Rdata")
rm(list=setdiff(ls(), c("all_data", "hits", "pc_data")))

# Load Crispr data from depmap
crispr_data <- read.csv("./DepMap/Achilles_gene_effect.csv",
                        header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Sample information
sample_info <- read.csv("./DepMap/sample_info.csv",
                        header = TRUE, stringsAsFactors = FALSE)

# Processed dataframe with cell line names for crispr data
proc_data <- crispr_data %>%
              left_join(select(sample_info, DepMap_ID, CCLE_Name), by = "DepMap_ID")

rownames(proc_data)   <- tolower(proc_data$CCLE_Name)
colnames(proc_data)   <- sapply(colnames(proc_data), function(x){strsplit(x, " ")[[1]][1]})

# Median centering the data
proc_data <- t(proc_data[-c(1,18121)])
proc_data <- t(apply(proc_data, 1, function(x){x - median(x, na.rm = TRUE)}))

# imputing missing values
proc_data <- as.data.frame(t(apply(proc_data,
                      1,
                      function(x){
                        x[which(is.na(x))] <- mean(x, na.rm = TRUE)
                        x})))

# Function for preparing slidr object for crispr data
getCRISPRdata <- function(slidr_obj){

  crispr_so <- slidr_obj

  common_KO        <- intersect(rownames(slidr_obj$viabilities), rownames(proc_data))
  common_celllines <- intersect(colnames(slidr_obj$mutations), colnames(proc_data))

  # Subset of common cell lines and genes for viability, mutations and copy number
  crispr_so$viabilities   <- proc_data[common_KO, common_celllines]
  crispr_so$mutations     <- crispr_so$mutations[ ,common_celllines]
  crispr_so$CNalterations <- crispr_so$CNalterations[ ,common_celllines]

  return(crispr_so)
}

# Getting slidr object for cancer type and cancer-type specific
crispr_pc_data <- getCRISPRdata(pc_data)

crispr_cs_data <- lapply(all_data, getCRISPRdata)

# running slidr on crispr data for pancancer data
path_results <- "./DepMap/Crispr/PanCan/"

hits_crispr_pc <- slidr::identifySLHits(canc_data = crispr_pc_data,
                                        path_results = path_results,
                                        WT_pval_thresh = 0,
                                        fp_thresh = 1)


# Running slidr on crispr data for cancer type specific
path_results <- "./DepMap/Crispr/CancerSpecific/"

hits_crispr_cs <- lapply(crispr_cs_data[1],
                         function(x){slidr::identifySLHits(canc_data = x,
                                                           path_results = path_results,
                                                           WT_pval_thresh = 0.1)})

names(hits_crispr_cs) <- names(crispr_cs_data)
