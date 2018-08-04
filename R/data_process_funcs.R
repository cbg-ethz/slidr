#' Retrieve cell lines for corresponding cancer type
#'
#' The function takes in the annotation file, meta data and the cancer type and returns
#' a set of all cell lines for that cancer type.
#'
#' @import tidyr
#' @param x primary site
#' @param cellline_annot cell line annotation data frame
#' @param meta_data information on different sub types for each primary site
#' @return A vector of cell lines for each cancer type
#' @export

getCelllines <- function(x, cellline_annot, meta_data){
  # Choosing only cellines for cancer of interest x
  celllines <- cellline_annot[cellline_annot$PRIMARY_SITE == x & cellline_annot$in.CCLE == "Y",]

  # Merging celllines of different subtypes
  if(subset(meta_data, Primary_site == x, select = Additional_filters, drop = TRUE) == "None"){
      celllines <- celllines$CELLLINE
  }else{
      celllines <- unlist(celllines %>%
                            dplyr::filter(PATHOLOGIST_ANNOTATION %in%
                                            unlist(strsplit(subset(meta_data,Primary_site==x,select=Additional_filters,drop=TRUE),"[;]"))) %>%
                            dplyr::select(CELLLINE))
  }
  return(celllines)
}


#' Prepare integer valued copy number matrix from Gistic based copy number
#'
#' The function takes in Gistic values of copy number as obtained from CCLE/TCGA,
#' cell lines of interest and top drivers and returns an integer valued copy
#' number matrix (0 = Homozygous deletion; 1 = Heterozygous deletion; 2 = Normal; etc.).
#'
#' @param CN_df_gistic Copy number in Gistic format
#' @param celllines vector of cell lines of interest
#' @param top_drivers vector of driver or mutated genes of interest
#' @param x primary site
#' @return Integer valued copy number matrix for selected cell lines and drivers
#' @export

prepareCNMatGistic <- function(CN_df_gistic, celllines, top_drivers, x){

  CCLE_cellline_name     <- toupper(paste(celllines,x, sep = "_"))
  rownames(CN_df_gistic) <-  CN_df_gistic$Hugo_Symbol
  CN_mat_subset <- CN_df_gistic[, -c(1:2)]

  # Converting the numbers to same format as other copy number function 0,1,2...
  CN_mat_subset <- CN_df_gistic %>% dplyr::select_(.dots = intersect(CCLE_cellline_name, colnames(CN_df_gistic))) + 2 #CN_df_gistic[,which(colnames(CN_df_gistic) %in% CCLE_cellline_name)] + 2
  colnames(CN_mat_subset) <- tolower(colnames(CN_mat_subset))

  # Integer values for copy numbers 0,1,2 ...
  CN_subset_int <-  matrix(NA,nrow=length(top_drivers),ncol=length(celllines))
  colnames(CN_subset_int) <- tolower(CCLE_cellline_name)
  rownames(CN_subset_int) <- top_drivers

  # Choosing only top driver genes
  indexes <-  which(top_drivers %in% rownames(CN_mat_subset))
  CN_subset_int[indexes, colnames(CN_mat_subset)] <- as.matrix(CN_mat_subset[top_drivers[indexes], ])

  return(CN_subset_int)
}


#' Get a list of essential genes
#'
#' The function takes in RSA data and returns a vector of essential genes
#'
#' @param data data frame of RSA normalised values
#' @param x primary site
#' @param celllines subset of interested celllines
#' @param thresh threshold for defining essentiality. Default = -3
#' @return a vector of essential genes
#' @export

getEssentialGenes <- function(x, celllines, data, thresh = -3){
  data_subset <- data[,paste(celllines,x,sep="_")]

  # Filling missing data with mean values
  data_subset <- as.data.frame(t(apply(data_subset,
                                    1,
                                    function(x){
                                      x[which(is.na(x))] <- mean(x, na.rm = TRUE)
                                      x})))

  # Use RSA data to get the essential genes with viability <= -3 in more than 50%
  essential_genes <- rownames(data_subset)[which(apply(data_subset,
                                                    1,
                                                    function(x){
                                                      if(sum(x<=thresh)>=round(length(x)/2)){
                                                        return(TRUE)
                                                      }else return(FALSE)}))]

  return(essential_genes)
}


#' Prepare integer valued copy number matrix
#'
#' The function takes in log2 ratios of copy number as obtained from CCLE,
#' cell lines of interest and top drivers and returns an integer valued copy
#' number matrix
#'
#' @param CN_df log2 transformed copy number matrix
#' @param celllines vector of cell lines of interest
#' @param top_drivers vector of driver or mutated genes of interest
#' @param x primary site
#' @return Integer valued copy number matrix for selected cell lines and drivers
#' @export

prepareCNMat <- function(CN_df, celllines, top_drivers, x){

  CCLE_cellline_name <- toupper(paste(celllines,x, sep = "_"))

  # Copy numbers for cell lines of interest
  CN_mat_subset <- CN_df %>% dplyr::select_(.dots = c("SYMBOL",intersect(CCLE_cellline_name, colnames(CN_df)))) #CN_df[,c(2, which(colnames(CN_df) %in% CCLE_cellline_name))]
  colnames(CN_mat_subset) <- tolower(colnames(CN_mat_subset))
  CN_mat_subset <- data.frame(CN_mat_subset[,-1], row.names=CN_mat_subset[,1])

  # Converting continuous values to integers
  CN_mat_subset  <- round(2*2^CN_mat_subset)

  # Integer values for copy numbers 0,1,2 ...
  CN_subset_int <-  matrix(NA,nrow=length(top_drivers),ncol=length(celllines))
  colnames(CN_subset_int) <- tolower(CCLE_cellline_name)
  rownames(CN_subset_int) <- top_drivers

  # Choosing only top driver genes
  indices <- which(top_drivers %in% rownames(CN_mat_subset))
  CN_subset_int[indices, colnames(CN_mat_subset)] <- as.matrix(CN_mat_subset[top_drivers[indices], ])

  return(CN_subset_int)
}

#' Prepare mutation matrix
#'
#' Function takes in a list of driver genes, samples and mutation/maf file
#' and returns a binary muration matrix
#'
#' @param x primary site
#' @param driver_genes vector of interested driver genes
#' @param samples vector of cell lines or patient samples
#' @param all_cancers_mut_df dataframe of mutations with Hugo_symbol, Tumor_Sample_Barcode and type of mutation
#' @return a mutation matrix for the given samples and genes
#' @export
#'
prepareMutMat <- function(x, driver_genes, samples, all_cancers_mut_df){

  # Generate mutation matrix of cell line and mutated genes file from each cell lines mutation file
  mut_mat <-  matrix(0,nrow=length(driver_genes),ncol=length(samples))

  for(i in 1:length(samples)){
    temp_mut <- subset(all_cancers_mut_df,
                       Tumor_Sample_Barcode == toupper(paste(samples[i],x,sep="_")),
                       Hugo_Symbol)
    mut_idx  <- which(driver_genes %in% unique(temp_mut$Hugo_Symbol))
    mut_mat[mut_idx,i] <- 1
  }

  colnames(mut_mat) <- paste(samples,x,sep="_")
  rownames(mut_mat) <- driver_genes

  return(mut_mat)
}

#' Prepare the data object required for downstream analysis
#'
#' The function processes the pan-cancer data and returns an object with viabilities matrix, mutation matrix,
#' mutation annotations and primary site for different types of cancers.
#'
#' @import tidyr
#' @import reshape2
#' @import biclust
#' @param data input data frame of cell line viabilities for different gene knockdowns
#' @param x primary site
#' @param fdr fdr cut-off for choosing the top drivers from mutSig2 list of drivers. Default = 0.05
#' @param min_Nmut lower bound of number of cell lines with mutations. Default = 2
#' @param all_cancers_mut_df MAF file from CCLE
#' @param CN_df copy number dataframe from CCLE
#' @param gistic Logical variable checking if copy number is based on Gistic. Default = FALSE
#' @param CN_Thr threshold for using CN data. Values: 0 = Homozygous and heterozygous deletions ; 1 = Homozygous deletions only; 2 = No copy number used (default)
#' @param minNrcelllines lower bound of number of cell lines. Default = 5
#' @param meta_data information on different sub types for each primary site
#' @param celllines vector of interested celllines
#' @param essential_genes vector of essential genes
#' @return An object for each cancer type
#' \describe{
#'    \item{viabilities}{dataframe of viabilities for each cancer type}
#'    \item{mutations}{matrix of mutations in drivers for each cancer type}
#'    \item{CNalterations}{matrix of non-negative copy number alterations of drivers for each cancer type}
#'    \item{mutation_annot}{annotations of the mutations}
#'    \item{primary_site}{cancer type}
#' }
#' @export
#'
#'
prepareDataObjects <- function(data, x, fdr = 0.05, min_Nmut = 2, all_cancers_mut_df, CN_df, gistic = FALSE,
                               CN_Thr = 2, minNrcelllines = 5, celllines, meta_data, essential_genes = NULL) {

  # Choosing only cellines for cancer of interest x
  #celllines <- cellline_annot[cellline_annot$PRIMARY_SITE == x & cellline_annot$in.CCLE == "Y",]

  if(length(celllines) < minNrcelllines) {
    return(NULL)
  }else{

    # Getting the top driver genes
    top_drivers <- NULL

    for(i in unlist(strsplit(subset(meta_data, Primary_site == x, Driver_gene_file, drop = TRUE),"[;]"))){
      # Loading driver mutation list from FireBrowse
      driver_genes <- read.delim(paste(base_folder, "DriverGenes/",i,sep=""),stringsAsFactors = F)
      # Choosing top driver genes with threshold of q <= fdr
      top_drivers <- c(top_drivers, driver_genes $gene[driver_genes $q <= fdr])
    }
    top_drivers <- unique(top_drivers)

    # Get mutation matrix for selected cell lines and top driver genes
    mut_mat <- prepareMutMat(x, top_drivers, celllines, all_cancers_mut_df)

    # Get the binary CN matrix of cell lines and  top driver genes.
    if(gistic == TRUE){
      CN_subset_int <- prepareCNMatGistic(CN_df, celllines, top_drivers, x)
    }else{
      CN_subset_int <- prepareCNMat(CN_df, celllines, top_drivers, x)
    }

    # threshold= 2 no copy number;threshold= 1 - only deep deletions; threshold= 0 - any deletion
    CN_subset_bin <- binarize(2-CN_subset_int, threshold = CN_Thr)
    CN_subset_bin[is.na(CN_subset_bin)] <- 0

    # Updating copy numbers in mutation matrix
    mut_mat <- binarize((CN_subset_bin + mut_mat), threshold = 0)

    # Filtering out genes not seen in any of the cell lines or cell lines without any mutations
    mut_mat <- mut_mat[,colSums(mut_mat) > 0]
    mut_mat <- mut_mat[rowSums(mut_mat) >= min_Nmut,]
    mut_mat <- mut_mat[(ncol(mut_mat) - rowSums(mut_mat)) >= min_Nmut, ,drop=FALSE]

    mut_mat_annot <- melt(mut_mat)
    mut_mat_annot <- mut_mat_annot %>% dplyr::filter(value != 0)
    mut_mat_annot <- cbind(mut_mat_annot,
                           t(mapply(function(a,b){
                             c(all_cancers_mut_df %>%
                                 dplyr::filter(Tumor_Sample_Barcode == toupper(b) & Hugo_Symbol == a) %>%
                                 dplyr::summarise(Variant_Classification = paste(Variant_Classification,collapse=";")),CN_subset_int[a,b])
                           }, as.character(mut_mat_annot$Var1), as.character(mut_mat_annot$Var2))))
    colnames(mut_mat_annot) <- c("Hugo_Symbol","Cell_Line","Mut_Status","Variant_Classification","CN_Value")

    data_cancer <- data[,paste(celllines,x,sep="_")]
    data_cancer <- as.data.frame(t(apply(data_cancer,
                                         1,
                                         function(x){
                                           x[which(is.na(x))] <- mean(x, na.rm = TRUE)
                                           x})))
    data_cancer <- data_cancer[!rownames(data_cancer) %in% essential_genes,]

    return(list(viabilities = data_cancer,
                mutations = mut_mat,
                CNalterations = CN_subset_int,
                mutation_annot = mut_mat_annot,
                primary_site = x))
  }
}

