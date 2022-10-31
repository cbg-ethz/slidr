# Script to compare to CRISPR literature gold-standard for the 2021 revisions.
# Note: Need to run Comparison_goldstd.R before running this.

library(tidyverse)

load("~/Downloads/Slidr_Results_new/GoldStandard/ComparisonGS.Rdata")
# Compiled CRISPR literature gold-standard
lit_goldstd <- read.delim("~/Downloads/Slidr_Results_new/GoldStandard/lit_goldstd.txt",
                          sep = "\t", header = TRUE )
# Creating a subset of dataset by collapsing references and cell lines
temp_df <- lit_goldstd %>%
            dplyr::rowwise() %>%
            dplyr::mutate(Pairs = paste0(sort(c(Gene_A, Gene_B)), collapse = '_')) %>%
            dplyr::distinct_at(vars(Pairs, Hit_Cell_Line, Authors))

# Combine the cell lines for common pairs within a reference
temp_df <- temp_df %>%
            dplyr::group_by_at(vars(Pairs,Authors))  %>%
            dplyr::summarise(Cell_line = paste0(Hit_Cell_Line, collapse = ", "))

# Combine references for common pairs
temp_df <- temp_df %>%
            dplyr::group_by(Pairs) %>%
            dplyr::summarise(across(Cell_line:Authors, ~ paste0(., collapse = "; "))) %>%
            separate(Pairs, c("Gene_A", "Gene_B"), sep = "_", remove = FALSE)

# Subset of the gold-standard dataframe overlapping with our list of possible interactions
lit_goldstd_sub1 <- temp_df %>%
                      dplyr::filter(Gene_A %in% rownames(pc_data$mutations)) %>%
                      dplyr::filter(Gene_B %in% rownames(pc_data$viabilities))

lit_goldstd_sub2 <- temp_df %>%
                      dplyr::filter(Gene_B %in% rownames(pc_data$mutations)) %>%
                      dplyr::filter(Gene_A %in% rownames(pc_data$viabilities))

lit_goldstd_sub <- rbind.data.frame(lit_goldstd_sub1, lit_goldstd_sub2)

# Remove duplicates
lit_goldstd_sub <- distinct(lit_goldstd_sub)

# All the common pairs
common_pairs <-   unique(intersect(paste(hits_gs$driver_gene, hits_gs$sl_partner_gene, sep = "_"),
                                   c(paste(lit_goldstd_sub$Gene_A, lit_goldstd_sub$Gene_B, sep = "_"),
                                     paste(lit_goldstd_sub$Gene_B, lit_goldstd_sub$Gene_A, sep = "_"))))

# C <- data.frame(t(apply(lit_goldstd_sub[1:2], 1, sort)))
# D <- C %>% dplyr::distinct()
# common_pairs <-   unique(intersect(paste(hits_gs$driver_gene, hits_gs$sl_partner_gene, sep = "_"),
#                                     paste(D[,1], D[,2], sep = "_")))
#
# common_pairs <- unique(c(common_pairs,
#                    unique(intersect(paste(hits_gs$driver_gene, hits_gs$sl_partner_gene, sep = "_"),
#                                   paste(D[,2], D[,1], sep = "_")))))

# lit_goldstd_final <- lit_goldstd_sub %>%
#                         dplyr::rowwise() %>%
#                         dplyr::mutate(Pairs = paste0(sort(c(Gene_A, Gene_B)), collapse = '_')) %>%
#                         dplyr::select(Pairs, Hit_Cell_Line, Authors) %>%
#                         dplyr::distinct()
#
# Test <- lit_goldstd_final %>%
#                         dplyr::group_by(Authors) %>%
#                         dplyr::distinct(across(contains("Pairs")), .keep_all = TRUE) %>%
#                         dplyr::group_by(Pairs) %>%
#                         dplyr::summarise(across(Hit_Cell_Line:Authors, ~ paste0(., collapse = "; ")))


# common_pairs_alt <- unique(intersect(lit_goldstd_sub$Pairs,
#                        c(paste(hits_gs$driver_gene, hits_gs$sl_partner_gene, sep = "_"),
#                        paste(hits_gs$sl_partner_gene, hits_gs$driver_gene, sep = "_"))))


# If we want to remove the duplicates in our possible pairs space
# namespace <- expand.grid(x = rownames(pc_data$mutations), y = rownames(pc_data$viabilities), stringsAsFactors = FALSE)
# colnames(namespace) <- c("Gene_A", "Gene_B")
# sorted_namespace <- namespace %>%
#                       dplyr::rowwise() %>%
#                       dplyr::mutate(Pairs = paste0(sort(c(Gene_A, Gene_B)), collapse = '_'))
# poss_pairs <- unique(sorted_namespace$Pairs)

Overlap   <- length(common_pairs)
group2    <- nrow(hits_gs) # All pairs with p-val < pval_thresh
Total     <- nrow(pc_data$viabilities) * nrow(pc_data$mutations)
# All the unique true positives since isle_sub has duplicate pairs in diff cancers
group1    <- nrow(lit_goldstd_sub)

phyper(Overlap-1, group2, Total-group2, group1, lower.tail= FALSE)

