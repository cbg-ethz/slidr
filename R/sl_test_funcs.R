#' CDF for Irwin-Hall distribution
#'
#' The function computes the CDF for Irwin-Hall distribution which is the p-value for the
#' rank of viablities in the mutated samples.
#'
#' @param x normalised rank sum of viabilities in the mutated samples
#' @param n number of mutated samples
#' @return the CDF value
#' @export

IH_CDF <- function(x, n) {
  if(n <= 20){
    X <-  floor(x)
    k <- seq(from = 0, to = X)
    # compute the cdf or p-value for n <= 100
    s <-  (-1)^k * choose(n, k)*( (x-k)^n)
    return(sum(s)/factorial(n))
  }else{
    # approximation for large n
    return(pnorm(x, n/2,sqrt(n/12), lower.tail = TRUE))
  }
}

#' Retrieve list of mutation-specific synthetic lethal partners for each type of cancer
#'
#' Rank test based on Irwin-Hall distribution to identify synthetic lethal (SL) partners
#' in the mutated samples. The function also performs a two sided wilcoxon or t-test
#' on the wild type samples to give a p-value required for filtering.
#'
#' @import tidyr
#' @import rDGIdb
#' @param canc_data Processed data object for a given cancer type
#' @param fp_thresh The average number of false positives allowed. Default = 1
#' @param path_results The path to where the results should be stored. Default working directory.
#' @param WT_pval_thresh Discard SL pairs with WT p-values less than this threshold. Default = 0.1
#' Required by `plotSLBoxplot` function
#' @return a dataframe of driver gene with its corresponding SL partner, the p-value in WT samples,
#' p-value in mutated samples, the corresponding value after scaling, and available drugs if any.
#' @export

identifySLHits <- function(canc_data, fp_thresh = 1, path_results = NULL, WT_pval_thresh = 0.1){

  if(is.null(path_results)){
    output_folder = paste(getwd(),"/Hit_List/", sep = "")
  }else
    output_folder = paste(path_results,"Hit_List/", sep = "")

  if(!dir.exists(output_folder)){
    dir.create(output_folder, recursive = TRUE, showWarnings = TRUE)
  }

  viabilities <- canc_data$viabilities
  mutations   <- canc_data$mutations

  # Creating an output dataframe
  results     <- data.frame(driver=character(),
                            sl_partner=character(),
                            WT_pvalue=double(),
                            mut_pvalue=double(),
                            sc_pvalue=double(),
                            stringsAsFactors=FALSE)

  for(driver_gene in rownames(mutations)) {
    # Choosing WT cell lines and mutated cell lines
    WT_celllines  <- colnames(mutations)[which(mutations[driver_gene,] == 0)]
    mut_celllines <- colnames(mutations)[which(mutations[driver_gene,] == 1)]

    # Ranking the viabilities of all genes in the mutated cell lines
    mut_viabilities_ranks <- apply(viabilities[, mut_celllines], 2, rank)
    # Normalised rank sum for mutated cell lines
    genes_rank_sum_mut    <- sort(apply(mut_viabilities_ranks, 1, sum))/nrow(viabilities)

    # Ranking the viabilities of all genes in the WT cell lines
    WT_viabilities_ranks <- apply(viabilities[, WT_celllines], 2, rank)
    # Normalised rank sum for WT cell lines
    genes_rank_sum_WT    <- sort(apply(WT_viabilities_ranks, 1, sum))/nrow(viabilities)

    # Number of mutated cases and the maximum rank or total number of samples
    n        <- length(mut_celllines)
    n_tests  <- nrow(viabilities) * nrow(mutations)
    max_rank <- max(mut_viabilities_ranks)
    # Number of WT cases
    m        <- length(WT_celllines)

    # Initialize
    k = 1
    mut_pvalue <- IH_CDF(genes_rank_sum_mut[k], n)
    # Perform tests until a given threshold on number of acceptable FP
    while((mut_pvalue * n_tests) < fp_thresh && k <= max_rank){
      sl_partner_gene <- names(genes_rank_sum_mut[k])
      sc_pvalue       <- mut_pvalue * n_tests # Scaling p-value to check if it's less than average FP
      # IH test for SL partner gene in WT cell lines
      WT_cdf          <- IH_CDF(genes_rank_sum_WT[sl_partner_gene], m)
      WT_pvalue       <- 2 * min(WT_cdf, 1-WT_cdf)

      results         <- rbind.data.frame(results,
                                          cbind(driver_gene,
                                                sl_partner_gene,
                                                WT_pvalue,
                                                mut_pvalue,
                                                sc_pvalue))
      k               <- k + 1
      if(k <= max_rank){
        mut_pvalue <- IH_CDF(genes_rank_sum_mut[k], n) # Performing Irwin Hall test for next gene
      }
    }

  }

  cat("Finished performing all tests successfully!\n")
  results$WT_pvalue       <- as.numeric(as.character(results$WT_pvalue))
  results$mut_pvalue      <- as.numeric(as.character(results$mut_pvalue))
  results$sc_pvalue      <- as.numeric(as.character(results$sc_pvalue))
  results$driver_gene     <- as.character(results$driver_gene)
  results$sl_partner_gene <- as.character(results$sl_partner_gene)
  results                 <- results[order(results$mut_pvalue),]
  drugs_df                <- detailedResults(queryDGIdb(unique(results$sl_partner_gene)))
  if(nrow(drugs_df) > 0){
    drugs_df                <- drugs_df %>%
                                dplyr::group_by(Gene) %>%
                                dplyr::summarise(Drugs  = paste(Drug, collapse = ","))
    drugs_df$Gene           <- as.character(drugs_df$Gene)
    results                 <- dplyr::left_join(results, drugs_df, by = c("sl_partner_gene" = "Gene"))
  }

  write.table(results,
              file = paste(output_folder, "SL_hits_", canc_data$primary_site, ".txt", sep = ""),
              quote = F, row.names = F, sep = "\t")

  if(!is.null(canc_data$CNalterations) && !is.null(canc_data$mutation_annot)){
    slidr::plotSLBoxplot(canc_data = canc_data,
                         hits_df = results,
                         path_results = path_results,
                         WT_pval_thresh = WT_pval_thresh)
  }

  return(results)
}


#' P-value for a given pair of genes
#'
#' Given the processed data object, a driver gene, and it's prospective partner, the
#' function returns the p-values for the wild type and mutated samples.
#'
#' @import tidyr
#' @param canc_data Processed data object for a given cancer type.
#' @param driver_gene The target mutated gene.
#' @param sl_partner_gene  The corresponding SL partner
#' @return a dataframe of driver gene with its corresponding SL partner, the p-value in WT samples,
#' p-value in mutated samples.
#' @export

getPval <- function(canc_data, driver_gene, sl_partner_gene){
  viabilities <- canc_data$viabilities
  mutations   <- canc_data$mutations

  WT_celllines  <- colnames(mutations)[which(mutations[driver_gene,] == 0)]
  mut_celllines <- colnames(mutations)[which(mutations[driver_gene,] == 1)]

  # Ranking the viabilities of all genes in the mutated cell lines
  mut_viabilities_ranks <- apply(viabilities[, mut_celllines], 2, rank)
  # Normalised rank sum for mutated cell lines
  genes_rank_sum_mut    <- sort(apply(mut_viabilities_ranks, 1, sum))/nrow(viabilities)

  # Ranking the viabilities of all genes in the WT cell lines
  WT_viabilities_ranks <- apply(viabilities[, WT_celllines], 2, rank)
  # Normalised rank sum for WT cell lines
  genes_rank_sum_WT    <- sort(apply(WT_viabilities_ranks, 1, sum))/nrow(viabilities)

  # Number of mutated cases and index of sl parnter in ranked list
  n        <- length(mut_celllines)
  # Number of WT cases and index of sl parnter in ranked list
  m        <- length(WT_celllines)

  # IH test for SL partner gene in mutated cell lines
  mut_pvalue      <- IH_CDF(genes_rank_sum_mut[sl_partner_gene], n)
  # IH test for SL partner gene in WT cell lines
  WT_cdf          <- IH_CDF(genes_rank_sum_WT[sl_partner_gene], m)
  WT_pvalue       <- 2 * min(WT_cdf, 1-WT_cdf)

  cbind.data.frame(driver_gene,
                   sl_partner_gene,
                   WT_pvalue,
                   mut_pvalue)
}
