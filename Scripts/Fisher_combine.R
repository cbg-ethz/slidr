# Script to get the robust hits from both crispr and DRIVE data using Fisher combined p-value

library(slidr)
library(dplyr)
setwd("~/Downloads/Slidr_Results_new/")

load("./DepMap/Crispr/crisprProcessedData.RData")

pc_hits  <- read.delim2("./PanCan8pc/Hit_List/SL_hits_pan_cancer.txt",
                        stringsAsFactors = FALSE, header = TRUE)

# Function to get fishers combined p-value
getCombPval <- function(crispr_data, rnai_data, crispr_hits, rnai_hits){

  # Splitting fusion transcripts in DRIVE to include only their main gene
  idx        <- grep(",.*-", rnai_hits$sl_partner_gene)
  temp_genes <- rnai_hits$sl_partner_gene
  temp_genes[idx] <- unlist(sapply(temp_genes[idx], function(x){strsplit(x, ",")[[1]][1]}))
  rnai_hits$alias_partner_gene <- temp_genes

  # Common pairs between both
  common_hits <- intersect(paste(crispr_hits$driver_gene, crispr_hits$sl_partner_gene, sep = "_"),
                           paste(rnai_hits$driver_gene, rnai_hits$alias_partner_gene, sep = "_"))

  # Taking the union of all the pairs in CRISPR and shRNA
  union_df <- full_join(rnai_hits, crispr_hits,
                        by = c("driver_gene" = "driver_gene",
                               "alias_partner_gene" = "sl_partner_gene")) %>%
              dplyr::select(driver_gene, sl_partner_gene, alias_partner_gene )

  # Setting the sl_partner_gene for crispr hits not in drive hits
  na_idx <- which(is.na(union_df$sl_partner_gene))
  union_df$sl_partner_gene[na_idx] <- union_df$alias_partner_gene[na_idx]

  # union_hits <- union(paste(crispr_hits$driver_gene, crispr_hits$sl_partner_gene, sep = "_"),
  #                     paste(rnai_hits$driver_gene, temp_genes, sep = "_"))
  #
  # union_df   <- do.call(rbind.data.frame, strsplit(union_hits, "_"))

  common_res <- data.frame()
  unique_res <- data.frame()

  # Get the Fishers Pvalue for all the SL pairs
  for(i in 1:nrow(union_df)){

    driver  <- union_df$driver_gene[i]
    partner <- union_df$sl_partner_gene[i]
    alias   <- union_df$alias_partner_gene[i] # Since Crispr doesn't have fusion transcripts

    # Getting the crispr KO and mutations
    rnai_KO    <- rownames(rnai_data$viabilities)
    crispr_KO  <- rownames(crispr_data$viabilities)
    crispr_mut <- rownames(crispr_data$mutations)

    # If the partner has only fusion transcripts in the DRIVE data then choose the first best match
    if(!partner %in% rnai_KO){
      partner <- rnai_KO[grep(paste0("^",partner), rnai_KO)][1]
    }

    # To check if the corresponding driver and partner are in crispr data
    if(alias %in% crispr_KO & driver %in% crispr_mut){
      # Get the p-values for each pair in each dataset
      # This part is redundant but used for sanity check
      temp <- cbind.data.frame(getPval(rnai_data, driver, partner),
                               getPval(crispr_data, driver, alias))
      # Get the pvals in mutated cell-lines
      mut_pvals <- temp[,c(4,8)]
      WT_pvals  <- temp[,c(3,7)]

      # Combine the pvalues using Fishers method
      temp$comb_WT_pval  <- pchisq(-2 * sum(log(WT_pvals)), df = length(WT_pvals) * 2, lower.tail = F)
      temp$comb_mut_pval <- pchisq(-2 * sum(log(mut_pvals)), df = length(mut_pvals) * 2, lower.tail = F)

      common_res <- rbind.data.frame(common_res, temp[,-5])
    }else{
      unique_res <- dplyr::bind_rows(unique_res,
                                     cbind.data.frame(getPval(rnai_data, driver, partner)))
    }
  }

  colnames(common_res) <- c("driver_gene", "sl_partner_gene", "drive_WT_pval", "drive_mut_pval",
                            "alias_partner_gene", "crispr_WT_pval", "crispr_mut_pval",
                            "comb_WT_pval","comb_mut_pval")

  # Scaling combined pvalue in mutated cases by total number of tests in drive
  total_tests <- nrow(rnai_data$viabilities) * nrow(rnai_data$mutations)
  common_res$comb_sc_pval <- common_res$comb_mut_pval * total_tests

  # Sometimes there are no unique hits specific to only DRIVE
  if(nrow(unique_res) > 0){
    colnames(unique_res) <- c("driver_gene", "sl_partner_gene", "drive_WT_pval", "drive_mut_pval")
    robust_res           <- dplyr::bind_rows(common_res, unique_res)
  }else
    robust_res           <- common_res


  return(list(robust_hits = robust_res,
              common_hits = common_hits))
}


pc_robust <- getCombPval(crispr_data = crispr_pc_data,
                         rnai_data = pc_data,
                         crispr_hits = hits_crispr_pc,
                         rnai_hits = pc_hits)

temp <- pc_robust$robust_hits %>%
          dplyr::filter(comb_sc_pval < 1)

write.table(temp,
            file = "./DepMap/Crispr/RobustHits/robust_hits_pancancer.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

cs_robust <- lapply(names(crispr_cs_data), function(x){
                                            getCombPval(crispr_data = crispr_cs_data[[x]],
                                                        rnai_data = all_data[[x]],
                                                        crispr_hits = hits_crispr_cs[[x]],
                                                        rnai_hits = hits[[x]])})

names(cs_robust) <- names(crispr_cs_data)

# writing the robust hits after filtering based on thresholds
for(i in names(cs_robust)){

  temp <- cs_robust[[i]]$robust_hits %>%
            dplyr::filter(comb_sc_pval < 1) %>%
            dplyr::filter(drive_WT_pval >= 0.1 | crispr_WT_pval >= 0.1) %>%
            dplyr::arrange(comb_sc_pval)

  write.table(temp,
              file = paste0("./DepMap/Crispr/RobustHits/robust_hits_",i,".txt"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
}

save.image("./DepMap/Crispr/RobustHits/FishersRobustHits.Rdata")
