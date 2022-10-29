# Script to simulate data and assess performance of slidr

library(slidr)
library(dplyr)


# Load processed data
load("~/Downloads/Slidr_Results_new/ContCN/ProcessedData.Rdata")
rm(list=setdiff(ls(), "all_data"))

#####################
#  DATA SIMULATION #
####################

# Function for simulating data
simulateData <- function(site, n_sl){
  # Using the mutation matrix for liver and bone cancer
  mut_data   <- all_data[[site]]$mutations
  n_mut      <- dim(mut_data)[1]
  mut_genes  <- rownames(mut_data)
  cells      <- colnames(mut_data)
  n_cells    <- length(cells)

  # perturbed genes
  n_pert     <- dim(all_data[[site]]$viabilities)[1]
  pert_genes <- paste0("PG",1:n_pert)

  # Viability summary statistics
  v_summ  <- summary(unlist(all_data[[site]]$viabilities))
  v_sd    <- sd(unlist(all_data[[site]]$viabilities))

  # Sampling the wild type distribution
  pert_data <- matrix(rnorm(n_pert * n_cells , v_summ[4], v_sd),
                      n_pert,
                      n_cells)
  dimnames(pert_data) <- list(pert_genes, cells)


  # True SL pairs for each cancer
  true_sl <- cbind(sample(mut_genes, n_sl, replace = TRUE),
                   sample(pert_genes, n_sl, replace = FALSE))

  # Off-setting the viabilities of true sl pairs
  for(j in 1:n_sl){
    mut <- true_sl[j,1]
    tar <- true_sl[j,2]
    mut_cells <- cells[which(mut_data[mut,] == 1)]
    # mean for the new distribution as the mean between min and mean of WT
    mut_mean <- (v_summ[1]+v_summ[4])/2
    pert_data[tar, mut_cells] <- rnorm(length(mut_cells), mut_mean, 1.2*v_sd)
  }

  temp <- list(mutations = mut_data,
               viabilities = pert_data,
               CNalterations = NULL,
               mutation_annot = NULL,
               primary_site = site)

  return(list(data = temp, true_sl = true_sl))
}

# Simulate 30 sl for liver and bone
RNGkind(sample.kind = "Rounding")
set.seed(56429) # use: RNGkind(sample.kind = "Rounding") if running on R >= 3.6
n_sl  <- 30
sites <-  c("liver", "bone", "ovary")
sim_data <- mapply(simulateData, sites, rep(n_sl, 3), SIMPLIFY = FALSE)
names(sim_data) <- sites

#####################
#  COMPARISONS     #
####################

fp_thresh = 1

# Run slidr
path_results <- "~/Downloads/"

slidr_hits  <- lapply(sim_data,
               function(x){slidr::identifySLHits(canc_data = x$data,
                                                 fp_thresh = Inf,
                                                 path_results = path_results,
                                                 WT_pval_thresh = 0.1)})
names(slidr_hits) <- sites

filt_slidr_hits   <- lapply(slidr_hits, function(x){ x %>%
                                        dplyr::filter(sc_pvalue < fp_thresh & WT_pvalue >= 0.1)})

# Run Wilcoxon
Wilcox_hits  <- lapply(sim_data,
                       function(x){slidr::WilcoxSLHits(canc_data = x$data)})

names(Wilcox_hits) <- sites

# save.image("~/Downloads/Slidr_Results_new/SimulationStudy/Wilcox_Slidr_comp.Rdata")


# Filtering those that are more than FP_threshold
fpr_wilcox   <- lapply(Wilcox_hits, function(x){ x %>%
                                  dplyr::filter(sc_pvalue < fp_thresh*50 & WT_pvalue >= 0.1)})

# Filtering based on FDR
fdr_wilcox   <- Wilcox_hits

# Compute FDR
for(site in sites){
  fdr_wilcox[[site]]$fdr <- p.adjust(fdr_wilcox[[site]]$mut_pvalue, method = "fdr")
}

fdr_wilcox <- lapply(fdr_wilcox, function(x){ x %>%
                             dplyr::filter(fdr < 0.2 & WT_pvalue >= 0.1)})

# Run t-test
ttest_hits  <- lapply(sim_data,
                       function(x){slidr::ttestSLHits(canc_data = x$data)})

names(ttest_hits) <- sites

# Filtering those that are more than FP_threshold
fpr_ttest   <- lapply(ttest_hits, function(x){ x %>%
                      dplyr::filter(sc_pvalue < fp_thresh*50 & WT_pvalue >= 0.1)})

save.image("~/Downloads/Slidr_Results_new/Simulation_comparisons.Rdata")

#####################
#   EVALUATION      #
####################
# Function for F1 score
getF1 <- function(hits, site){

  true_pairs <- paste(sim_data[[site]]$true_sl[,1],
                      sim_data[[site]]$true_sl[,2],
                      sep = "_")
  pred_pairs <- paste(hits$driver_gene,
                      hits$sl_partner_gene,
                      sep = '_')
  # if(length(pred_pairs) < length(true_pairs)){
  #   pred_pairs <- c(pred_pairs, rep("none", length(true_pairs) - length(pred_pairs)))
  # }else{
  #   true_pairs <- c(true_pairs, rep("none", length(pred_pairs) - length(true_pairs)))
  # }

  prec <- sum(pred_pairs %in% true_pairs)/length(pred_pairs)
  rec  <- sum(pred_pairs %in% true_pairs)/n_sl
  f1   <- 2 * (prec * rec)/(prec + rec)
  #if(is.nan(f1)) browser()
  return(f1)
}

f1scores_df <- rbind.data.frame(mapply(getF1, filt_slidr_hits, sites),
                                mapply(getF1, fpr_wilcox, sites),
                                mapply(getF1, fpr_ttest, sites))

colnames(f1scores_df) <- sites


##########################################################
############## Plotting PR curve #########################
##########################################################

library(pROC)

plabs        <- letters[1:3]
names(plabs) <- sites

pdf("./SimulationStudy/Performance_plots.pdf", width = 10, height = 6)

layout(matrix(c(1, 3, 5, 2, 4, 6, 7, 7, 7), ncol=3, byrow = TRUE), heights = c(4, 4, 0.8))

par(mar = c(4, 4, 2, 2))

for(site in sites){
  # Ground truth SL pairs
  true_pairs <- paste(sim_data[[site]]$true_sl[,1],
                      sim_data[[site]]$true_sl[,2],
                      sep = "_")
  # All possible pairs
  drivers    <- rownames(sim_data[[site]]$data$mutations)
  pert_genes <- rownames(sim_data[[site]]$data$viabilities)
  all_pairs  <- sort(apply(expand.grid(drivers, pert_genes), 1, paste, collapse = "_"))

  # The true outcomes vector
  n        <- nrow(sim_data[[site]]$data$mutations) * nrow(sim_data[[site]]$data$viabilities)
  outcomes <- rep(0, n)
  names(outcomes) <- all_pairs
  outcomes[which(names(outcomes) %in% true_pairs)] <- 1

  # Constructing the dataframes for each method
  pred_slidr <- cbind.data.frame(paste(slidr_hits[[site]]$driver_gene,
                                       slidr_hits[[site]]$sl_partner_gene, sep = '_'),
                                 slidr_hits[[site]]$sc_pvalue)
  colnames(pred_slidr) <- c("pairs", "scores")

  pred_wilcox <- cbind.data.frame(paste(Wilcox_hits[[site]]$driver_gene,
                                        Wilcox_hits[[site]]$sl_partner_gene, sep = '_'),
                                  Wilcox_hits[[site]]$sc_pvalue)
  colnames(pred_wilcox) <- c("pairs", "scores")

  pred_ttest <- cbind.data.frame(paste(ttest_hits[[site]]$driver_gene,
                                       ttest_hits[[site]]$sl_partner_gene, sep = '_'),
                                 ttest_hits[[site]]$sc_pvalue)
  colnames(pred_ttest) <- c("pairs", "scores")

  # Sorting so that the pair order in prediction matches ground truth
  pred_slidr  <- pred_slidr %>% dplyr::arrange(pairs)
  pred_wilcox <- pred_wilcox %>% dplyr::arrange(pairs)
  pred_ttest  <- pred_ttest %>% dplyr::arrange(pairs)

  # Generate the roc objects
  roc_obj_slidr  <- roc(outcomes, pred_slidr$scores, percent = FALSE)
  roc_obj_wilcox <- roc(outcomes, pred_wilcox$scores, percent = FALSE)
  roc_obj_ttest <- roc(outcomes, pred_ttest$scores, percent = FALSE)

  # Plotting the roc curves and pr curves
  plot((1-roc_obj_slidr$specificities), roc_obj_slidr$sensitivities,
       type = "l", col = "#437BAA", xlab = "1 - Specificity", ylab = "Sensitivity",
       lwd = 2, cex.lab = 1.2, cex.axis = 1.2, main = toupper(site))
  lines((1-roc_obj_wilcox$specificities), roc_obj_wilcox$sensitivities,
        type = "l", col = "#BC3B3B", lwd = 2)
  lines((1-roc_obj_ttest$specificities), roc_obj_ttest$sensitivities,
        type = "l", col = "#6FAA2F", lwd = 2)
  mtext(plabs[site], side = 3, line = 0.6, cex = 0.8, adj = -0.2, col = "black", font = 2)
  #title(outer = TRUE, adj = 0, main = "A", cex = 1.1, col = "black", font = 2, line=-1)

  plot(precision ~ recall,
       coords(roc_obj_slidr, "all", ret = c("recall", "precision"), transpose = FALSE),
       type = "l", col = "#437BAA", lwd = 2, xlab = "Recall", ylab = "Precision",
       cex.lab = 1.2, cex.axis = 1.2)#, main = toupper(site))
  lines(precision ~ recall,
        coords(roc_obj_wilcox, "all", ret = c("recall", "precision"), transpose = FALSE),
        col = "#BC3B3B", lwd = 2)
  lines(precision ~ recall,
        coords(roc_obj_ttest, "all", ret = c("recall", "precision"), transpose = FALSE),
        col = "#6FAA2F", lwd = 2)
  abline(h = 0.5 , col = "black", lwd = 2, lty = 2)

}

par(mai=c(0,0,0,0))
plot.new()
legend(x = "center", legend=c("SLIdR","t-test", "Wilcoxon test"),
       col=c("#437BAA", "#6FAA2F", "#BC3B3B"), lwd =  2, cex = 1.2,
       xpd = TRUE, horiz = TRUE, bty = "n", title="Methods")

dev.off()
