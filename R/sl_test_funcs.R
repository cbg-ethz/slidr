#' CDF for Irwin-Hall distribution
#' The function computes the CDF for Irwin-Hall distribution which is the p-value for the
#' rank of viablities in the mutated samples.
#' @param x normalised rank sum of viabilities in the mutated samples
#' @param n number of mutated samples
#' @return the CDF value
#' @export

IH_CDF <- function(x, n) {
  X <-  floor(x)
  k <- seq(from = 0, to = X)
  # compute the cdf which would be the p-value
  s <-  (-1)^k * choose(n, k)*( (x-k)^n)
  return(sum(s)/factorial(n))
}

#' Retrieve list of mutation-specific synthetic lethal partners for each type of cancer
#' Rank test based on Irwin-Hall distribution to identify synthetic lethal (SL) partners
#' in the mutated samples. The function also performs a two sided wilcoxon or t-test
#' on the wild type samples to give a p-value required for filtering.
#' @import tidyr
#' @param canc_data Processed data object for a given cancer type
#' @param n_cand first 'n' number of SL partners to be tested after ranking. Default = 100
#' This is just to save on computation time since after multiple testing correction we will be
#' filtering a lot of them out.
#' @param qval_thresh The number of false positives allowed after correction. Default = 1
#' @param path_results The path to where the results should be stored
#' @param WT_pval_thresh Discard SL pairs with WT p-values less than this threshold. Default = 0.2
#' Required by `plotSLBoxplot` function
#' @return a dataframe of driver gene with its corresponding SL partner, the p-value in WT samples,
#' p-value in mutated samples and the corresponding value after correction.
#' @export

identifySLHits <- function(canc_data, n_cand = 100, qval_thresh = 1, path_results, WT_pval_thresh = 0.2){

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
                            mut_qvalue=double(),
                            stringsAsFactors=FALSE)

  for(driver_gene in rownames(mutations)) {
    # Choosing WT cell lines and mutated cell lines
    WT_celllines  <- colnames(mutations)[which(mutations[driver_gene,] == 0)]
    mut_celllines <- colnames(mutations)[which(mutations[driver_gene,] == 1)]

    # Ranking the viabilities of all genes in the mutated cell lines
    all_viabilities_ranks <- apply(viabilities[, mut_celllines], 2, rank)
    # Normalised rank sum
    genes_rank_sum        <- sort(apply(all_viabilities_ranks, 1, sum))/nrow(viabilities)

    # Number of mutated cases and the maximum rank or total number of samples
    n        <- length(mut_celllines)
    max_rank <- max(all_viabilities_ranks)

    results <- rbind.data.frame(results,
                                do.call(rbind,
                                lapply(1:n_cand,function(k){
                                  sl_partner_gene <- names(genes_rank_sum[k])
                                  mut_pvalue      <- IH_CDF(genes_rank_sum[k], n) # Performing Irwin Hall test
                                  mut_qvalue      <- mut_pvalue * nrow(viabilities) * ncol(mutations) # Correcting for multiple testing
                                  x               <- as.numeric(as.character(viabilities[sl_partner_gene,WT_celllines]))
                                  WT_pvalue       <- min(wilcox.test(x, mu = 0, alternative = "two.sided")$p.value,
                                                       t.test(x, mu = 0, alternative = "two.sided")$p.value) # Testing WT cell lines for sl partner gene
                                  cbind(driver_gene,
                                      sl_partner_gene,
                                      WT_pvalue,
                                      mut_pvalue,
                                      mut_qvalue)})))
  }

  results$WT_pvalue       <- as.numeric(as.character(results$WT_pvalue))
  results$mut_pvalue      <- as.numeric(as.character(results$mut_pvalue))
  results$mut_qvalue      <- as.numeric(as.character(results$mut_qvalue))
  results$driver_gene     <- as.character(results$driver_gene)
  results$sl_partner_gene <- as.character(results$sl_partner_gene)
  results                 <- results %>% dplyr::filter(mut_qvalue < qval_thresh)
  results                 <- results[order(results$mut_pvalue),]

  slidR::plotSLBoxplot(canc_data = canc_data,
                       hits_df = results,
                       path_results = path_results,
                       WT_pval_thresh = 0.2)
  write.table(results,
              file = paste(output_folder, "SL_hits_", canc_data$primary_site, ".txt", sep = ""),
              quote = F, row.names = F, sep = "\t")
  return(results)
}
