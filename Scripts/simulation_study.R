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
set.seed(56429)
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
                                                 fp_thresh = fp_thresh,
                                                 path_results = path_results,
                                                 WT_pval_thresh = 0.1)})
names(slidr_hits) <- sites

slidr_hits        <- lapply(slidr_hits, function(x){ x %>%
                                                    dplyr::filter(WT_pvalue >= 0.1)})

# Run Wilcoxon
Wilcox_hits  <- lapply(sim_data,
                       function(x){slidr::WilcoxSLHits(canc_data = x$data)})

names(Wilcox_hits) <- sites

save.image("~/Downloads/Slidr_Results_new/SimulationStudy/Wilcox_Slidr_comp.Rdata")


# Filtering those that are more than FP_threshold
fpr_hits   <- lapply(Wilcox_hits, function(x){ x %>%
                                  dplyr::filter(sc_pvalue < fp_thresh*50 & WT_pvalue >= 0.1)})

# Filtering based on FDR
fdr_hits   <- Wilcox_hits

# Compute FDR
for(site in sites){
  fdr_hits[[site]]$fdr <- p.adjust(fdr_hits[[site]]$mut_pvalue, method = "fdr")
}

fdr_hits <- lapply(fdr_hits, function(x){ x %>%
                             dplyr::filter(fdr < 0.2 & WT_pvalue >= 0.1)})


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

f1scores_df <- rbind.data.frame(mapply(getF1, slidr_hits, sites),
                                mapply(getF1, fpr_hits, sites))

colnames(f1scores_df) <- sites
