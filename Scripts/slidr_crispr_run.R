# Crispr analysis

library(slidr)
library(dplyr)
setwd("~/Downloads/Slidr_Results_new/")

# Load pan cancer object
load("./PanCan8pc/pancan.Rdata")
rm(list=setdiff(ls(), c("all_data", "hits", "pc_data")))

# Load Crispr data from depmap
crispr_dat <- read.csv("./DepMap/Achilles_gene_effect.csv",
                        header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Sample information
sample_info <- read.csv("./DepMap/sample_info.csv",
                        header = TRUE, stringsAsFactors = FALSE)

# Dataframe with cell line names for crispr data
temp <- crispr_dat %>%
            left_join(select(sample_info, DepMap_ID, CCLE_Name), by = "DepMap_ID")

rownames(temp)   <- tolower(temp$CCLE_Name)
colnames(temp)   <- sapply(colnames(temp), function(x){strsplit(x, " ")[[1]][1]})

# imputing missing values
temp <- t(temp[-c(1,18121)])
temp <- as.data.frame(t(apply(temp,
                      1,
                      function(x){
                        x[which(is.na(x))] <- mean(x, na.rm = TRUE)
                        x})))

# Preparing slidr object for crispr data
crispr_pc_data   <- pc_data
common_KO        <- intersect(rownames(pc_data$viabilities), rownames(temp))
common_celllines <- intersect(colnames(pc_data$viabilities), colnames(temp))

# Subset of common cell lines and genes for viability, mutatins and copy number
crispr_pc_data$viabilities   <- temp[common_KO, common_celllines]
crispr_pc_data$mutations     <- crispr_pc_data$mutations[ ,common_celllines]
crispr_pc_data$CNalterations <- crispr_pc_data$CNalterations[ ,common_celllines]

# running slidr on crispr data
path_results <- "./DepMap/"

hits_crispr_pc <- slidr::identifySLHits(canc_data = crispr_pc_data,
                                        path_results = path_results,
                                        WT_pval_thresh = 1e-6,
                                        fp_thresh = 1)
