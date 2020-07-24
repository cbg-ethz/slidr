# Comparing with ISLE

library(slidr)
library(dplyr)
setwd("~/Downloads/Slidr_Results_new/")

# Dataframe for gold standard. Retaining known true SL
isle_df <- read.delim("./GoldStandard/ISLE_goldstandard.txt", stringsAsFactors = FALSE)
isle_df <- isle_df %>% dplyr::filter(SL == 1)
#  changing the order to ensure gene 1 is always the mutated gene
temp_df <- isle_df %>%
            dplyr::filter(gene2.perturbation == "mut") %>%
            dplyr::relocate(gene2, gene1, gene2.perturbation, gene1.perturbation)
colnames(temp_df) <- colnames(isle_df)
isle_df_new <- rbind(isle_df %>%
                       dplyr::filter(gene2.perturbation != "mut"),
                     temp_df)


########################
# Cancer type-specific #
########################
# NOTE: Skip this part and directly move to pan-cancer while running this code.
# This is only for checking the cancer type specific hits and is exploratory only.

# Load the cancer-specific data
load("./ContCN/ProcessedData.Rdata")
rm(list=setdiff(ls(), c("all_data", "hits")))

# Running Slidr for different cancer types
canc_types <- cbind.data.frame(unique(isle_df_new$cancer.type.tested),
                    c("large_intestine", "ovary", "endometrium", "lung",
                      "breast", "skin", "kidney", "haematopoietic_and_lymphoid_tissue"))

names(canc_types) <- c("tcga_site", "drive_site")

# Comparing for COAD
x <- canc_types$tcga_site[1]
y <- canc_types$drive_site[1]
# Subset of COAD and genes overlapping with our list of mutations and perturbed genes
isle_sub <- isle_df_new %>%
              dplyr::filter(cancer.type.tested %in% x) %>%
              dplyr::filter(gene1 %in% rownames(all_data[[y]]$mutations)) %>%
              dplyr::filter(gene2 %in% rownames(all_data[[y]]$viabilities))


# Checking the p-value for all the gold standard pairs in our dataset
pval_df <- NULL

for(i in 1:nrow(isle_sub)){
  pval_df <- rbind.data.frame(pval_df,
                              slidr::getPval(all_data[[y]],isle_sub$gene1[i], isle_sub$gene2[i]))
}

pval_df$sc_pvalue <- pval_df$mut_pvalue * nrow(all_data[[y]]$viabilities) * nrow(all_data[[y]]$mutations)



###############
# Pan cancer #
###############

# Load the cancer-specific data
load("./PanCan8pc/pancan.Rdata")
rm(list=setdiff(ls(), c("pc_data", "isle_df", "isle_df_new")))

# Load pan-cancer hits
hits_pancan <- read.delim("./PanCan8pc/Hit_List/SL_hits_pan_cancer.txt",stringsAsFactors = F)

# Subset of gold standard hits overlapping with our list of mutations and perturbed genes
isle_sub <- isle_df_new %>%
              dplyr::filter(gene1 %in% rownames(pc_data$mutations)) %>%
              dplyr::filter(gene2 %in% rownames(pc_data$viabilities))


# Checking the p-value for all the gold standard pairs in our dataset
pval_df <- NULL

for(i in 1:nrow(isle_sub)){
  pval_df <- rbind.data.frame(pval_df,
                              slidr::getPval(pc_data,isle_sub$gene1[i], isle_sub$gene2[i]))
}

pval_df$sc_pvalue <- pval_df$mut_pvalue * nrow(pc_data$viabilities) * nrow(pc_data$mutations)

# Use slidr to get the top hits with mut p-val less than 5%
pval_thresh <- 0.05
avg_fp <- round(pval_thresh *  nrow(pc_data$viabilities) * nrow(pc_data$mutations))

# Removing the mut and copy number annotation to skip the plotting step in SLIDR
new_data <- pc_data
new_data$CNalterations <- NULL
new_data$mutation_annot <- NULL

# CAUTION: Computationally intensive step. Takes 24 Hrs
hits_gs <- slidr::identifySLHits(canc_data = new_data,
                              path_results = "./GoldStandard/",
                              WT_pval_thresh = 0,
                              fp_thresh = avg_fp)

save.image("./GoldStandard/ComparisonGS.Rdata")

# Hyper-geometric test to see the degree of overlap
common_pairs <- unique(intersect(paste(hits_gs$driver_gene, hits_gs$sl_partner_gene, sep = "_"),
                                 paste(pval_df$driver_gene, pval_df$sl_partner_gene, sep = "_")))
Overlap   <- length(common_pairs)
group2    <- nrow(hits_gs) # All pairs with p-val < pval_thresh
Total     <- nrow(pc_data$viabilities) * nrow(pc_data$mutations)
group1    <- nrow(isle_sub) # All the true positives

phyper(Overlap-1, group2, Total-group2, group1, lower.tail= FALSE)
