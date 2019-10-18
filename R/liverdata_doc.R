#' SLIdR object for liver cancer.
#'
#' A dataset containing the viabilities, mutation profiles,
#' copy number alterations, and  mutation annotation dataframes
#' for liver cancer.
#'
#' @format A list with 5 data structures:
#' \describe{
#'   \item{viabilities}{ATARiS viability scores of 6393 perturbed genes (rows) across 15 cell lines (cols).}
#'   \item{mutations}{mutation matrix of 7 driver genes (rows) across 13 cell lines (cols).}
#'   \item{CNalterations}{copy number matrix of 54 genes (rows) across 15 cell lines (cols).}
#'   \item{mutation_annot}{data frame of mutation annotations of all the cell lines for the box plots.}
#'   \item{primary_site}{Name of the primary site.}
#' }
"LiverData"
