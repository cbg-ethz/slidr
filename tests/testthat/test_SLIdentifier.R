context("SL pairs idenitification")
library(slidr)
library(rDGIdb)

#  Creating toy dataset
test_data <- NULL
set.seed(496)

n_celllines <- 10
n_perturbs  <- 4000
n_mutations <- 4

celllines   <- paste0(rep("CL", n_celllines), 1:n_celllines)
perturbs    <- paste0(rep("PG", n_perturbs), 1:n_perturbs)
mutations   <- paste0(rep("MG", n_mutations), 1:n_mutations)

# Creating toy mutation matrix
mut_mat                      <- matrix(0,n_mutations, n_celllines)
dimnames(mut_mat)            <- list(mutations,celllines)
mut_mat["MG1",c("CL2", "CL3", "CL4", "CL5")] <- 1
mut_mat["MG2",c("CL7", "CL9", "CL10")]       <- 1
mut_mat["MG3",c("CL1", "CL6", "CL5", "CL8")] <- 1
mut_mat["MG4",c("CL1", "CL3", "CL10")]   <- 1
test_data$mutations          <- mut_mat

# Ceating toy viabilities matrix
per_mat                    <- matrix(rnorm(n_celllines * n_perturbs, 0, 0.4), n_perturbs, n_celllines)
dimnames(per_mat)          <- list(perturbs,celllines)
mut1_celllines             <- celllines[which(mut_mat[1,] == 1)]
per_mat[8, mut1_celllines] <- rnorm(length(mut1_celllines), -2., 0.25)
test_data$viabilities      <- as.data.frame(per_mat)


test_that("identifySLHits function works correctly", {
  expect_equal(identifySLHits(test_data)$mut_pvalue,
               IH_CDF(0.001, 4)) # 0.001 comes from sum of rank of PG8 for all 4 cases
  expect_equal(identifySLHits(test_data)$WT_pvalue,
               2 * min(IH_CDF(3.8155, 6), 1 - IH_CDF(3.8155, 6))) # 3.8155 comes from sum of rank of PG8 for all 6 WT cases
})

test_that("getPval works correctly", {
  expect_equal(getPval(test_data, "MG1", "PG8")$mut_pvalue,
               identifySLHits(test_data)$mut_pvalue)
  expect_equal(getPval(test_data, "MG1", "PG8")$WT_pvalue,
               identifySLHits(test_data)$WT_pvalue)
  expect_equal(getPval(test_data, "MG4", "PG3658")$mut_pvalue, IH_CDF(0.1295, 3)) # 0.1295 is sum of ranks for PG3658 in MG4 mutated cell lines
  expect_equal(getPval(test_data, "MG4", "PG3658")$WT_pvalue,
               2 * min(IH_CDF(4.689, 7), 1 - IH_CDF(4.689, 7)))# 4.689 is sum of ranks for PG3658 in MG4 WT cell lines
})
