# Script to FDR correct and use controls for PRISM validation

library(dplyr)
library(tidyr)
library(parallel)

load("~/Downloads/Slidr_Results_new/DepMap/Prism/29Nov2021/prism_validation.Rdata")

valid_ls <- cs_valid
valid_ls$pancancer <- pc_valid

# FDR correction
for(tissue in names(valid_ls)){

  res_df     <- valid_ls[[tissue]]$valid_hits
  res_df$fdr <- p.adjust(res_df$p_val, method = "fdr")

  if(nrow(res_df) > 0){
    res_df <- res_df %>%
                dplyr::filter(p_val < 0.1) %>%
                dplyr::arrange(p_val) %>%
                dplyr::select(driver_gene, sl_partner_gene,
                    WT_pvalue, mut_pvalue,
                    broad_id, name, target, p_val, fdr)

    colnames(res_df) <- c("Driver gene","SL partner gene",
                          "WT pvalue",	"mut pvalue", "Broad ID",
                          "Drug name", "Drug target(s)", "Validation pvalue", "FDR")
    # write.table(x = res_df,
    #             file = paste0("~/Downloads/Slidr_Results_new/DepMap/Prism/30Nov2021_FDR/prism_validation_", tissue, ".txt"),
    #             quote = FALSE,
    #             row.names = FALSE,
    #             sep = "\t")
  }
}


##### Compare with random controls #####
set.seed(83926) #85976

# Function to compare the fraction of significant hits from random controls
compare_control <- function(canc_data, cancer_type, seed_val , nruns = 10000){

  set.seed(seed_val)

  # Initialize for each cancer type
  frac_hits_ctrl <- integer(nruns)
  viabilities    <- canc_data$viabilities #all_data$pancreas$viabilities
  mutations      <- canc_data$mutations  #all_data$pancreas$mutations
  mut_genes      <- rownames(mutations)
  partner_genes  <- setdiff(rownames(viabilities), mut_genes)

  # retrieving the drugs for all possible SL partner genes
  subset_drugs <- do.call(rbind.data.frame,
                          sapply(partner_genes,
                                 function(x){drugs_red[which(grepl(paste0('\\b',x,'\\b'), drugs_red$target)),]},
                                 simplify = FALSE))

  subset_drugs$sl_target_gene <- sapply(row.names(subset_drugs), function(x){strsplit(x,"[.]")[[1]][1]})

  # All possible SL partner genes in DRIVE with drugs in PRISM
  poss_target_genes <- unique(subset_drugs$sl_target_gene)

  # All possible pairs
  poss_pairs <- tidyr::expand_grid(mut_genes, poss_target_genes) %>% dplyr::distinct()

  # Number of tests
  pred_df <- valid_ls[[cancer_type]]$valid_hits
  n_tests <- pred_df %>% nrow()
  n_pairs <- pred_df %>% dplyr::distinct(driver_gene, sl_partner_gene) %>% nrow()
  n_sig_hits <- pred_df %>% dplyr::filter(p_val < 0.1) %>% nrow()

  # # Testing in controls. Running nruns times
  for(i in 1:nruns){
    # Initialize variables
    results    <- data.frame()
    idx        <- sample(1:nrow(poss_pairs), n_pairs, replace = FALSE)
    cand_pairs <- poss_pairs[idx, ]

    # Testing for each candidate pair
    for(j in 1:nrow(cand_pairs)){
      target_gene <- cand_pairs$poss_target_genes[j]
      driver_gene <- cand_pairs$mut_genes[j]

      # Subset of drugs targeting sl_partner_gene
      temp <- dplyr::filter(subset_drugs, target_gene == sl_target_gene)

      # Stratifying based on driver gene mutation status
      WT_celllines  <- colnames(mutations)[which(mutations[driver_gene,] == 0)]
      mut_celllines <- colnames(mutations)[which(mutations[driver_gene,] == 1)]

      # Testing all the drugs for each target
      for(k in 1:nrow(temp)){
        compound <- temp$column_name[k]

        # viabilities in mutated and WT cell-lines from the drug
        mut_viability <- drug_screen[mut_celllines,compound]
        WT_viability  <- drug_screen[WT_celllines,compound]

        # Ensuring that the missing data is less than half
        if(sum(is.na(mut_viability)) < length(mut_celllines)/2 && sum(is.na(WT_viability)) < length(WT_celllines)/2){
          # Testing if the viability in mutated cell lines is significantly less than in WT from the drug
          p_val <- t.test(mut_viability, WT_viability, alternative =  "less")$p.value

          results <- rbind.data.frame(results,
                                      cbind(cand_pairs[j,], temp[k,], p_val))
        }
      }
    }
    # vector of number of significant hits
    # frac_hits_ctrl[i] <- results %>%
    #                     dplyr::filter(p_val < 0.1) %>%
    #                     #dplyr::distinct(mut_genes, poss_target_genes) %>%
    #                     nrow()

    # For each run store the fraction of significant hits
    temp_n_hits  <- results %>%
                    dplyr::filter(p_val < 0.1) %>% nrow()
    temp_n_tests <- results %>% nrow()
    frac_hits_ctrl[i] <- temp_n_hits/temp_n_tests

    # Empirical probability of finding more hits than observed at random
    emp_pval <- mean(frac_hits_ctrl >= (n_sig_hits/n_tests))

    return(list(frac_hits_ctrl = frac_hits_ctrl,
                n_tests = n_tests,
                n_sig_hits = n_sig_hits,
                emp_pval = emp_pval))
  }
}

# Cancer types with more than 0 significant hits in prism
sig_cancer_types <- valid_ls[sapply(valid_ls, function(x){x$valid_hits %>%
                                                           nrow() > 0})]

sig_cancer_types <- sapply(sig_cancer_types, function(x){x$valid_hits %>%
                                                          dplyr::filter(p_val < 0.1) %>%
                                                          nrow()})
sig_cancer_types <- sig_cancer_types[which(sig_cancer_types > 0)] %>% names()

# seeds vector
seed_vec <- sample.int(1e6, length(sig_cancer_types))

# List of data for cancer types with significant hits
data_ls           <- all_data
data_ls$pancancer <- pc_data
data_ls           <- data_ls[sig_cancer_types]

# Running for different cancers in parallel
# cluster <- makeCluster(4)
# result_ls <- clusterMap(cluster, compare_control,
#                         data_ls[1:4], sig_cancer_types[1:4], seed_vec[1:4])
result_ls <- mcmapply(compare_control,
                      data_ls, sig_cancer_types, seed_vec,
                      mc.set.seed = FALSE,
                      mc.cores = 6,
                      SIMPLIFY = FALSE)

sapply(result_ls, function(x){x$emp_pval})
# mean(frac_hits_ctrl >= (n_sig_hits/n_tests))
