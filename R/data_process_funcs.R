#' Get a list of essential genes
#'
#' The function takes in RSA data and returns a vector of essential genes
#'
#' @param data data frame of RSA normalised values
#' @param x primary site
#' @param thresh threshold for defining essentiality. Default = -3
#' @return a vector of essential genes
#' @export

getEssentialGenes <- function(x, data, thresh = -3){
  data_RSA <- data[,paste(celllines,x,sep="_")]
  
  # Filling missing data with mean values
  data_RSA <- as.data.frame(t(apply(data_RSA,
                                    1,
                                    function(x){
                                      x[which(is.na(x))] <- mean(x, na.rm = TRUE)
                                      x})))
  
  # Use RSA data to get the essential genes with viability <= -3 in more than 50%
  essential_genes <- rownames(data_RSA)[which(apply(data_RSA,
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
#' @param CN_mat log2 transformed copy number matrix 
#' @param celllines vector of cell lines of interest
#' @param top_drivers vector of driver or mutated genes of interest
#' @param x primary site
#' @return Integer valued copy number matrix for selected cell lines and drivers
#' @export

prepareCNMat <- function(CN_mat, celllines, top_drivers, x){
  
  CCLE_cellline_name <- toupper(paste(celllines,x, sep = "_"))
  
  # Copy numbers for cell lines of interest
  CN_mat_subset <- CN_mat[,c(2, which(colnames(CN_mat) %in% CCLE_cellline_name))]
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

#' Prepare the data object required for downstream analysis
#' 
#' The function processes the pan-cancer data and returns an object with viabilities matrix, mutation matrix,
#' mutation annotations and primary site for different types of cancers.
#' 
#' @param data input data frame of cell line viabilities for different gene knockdowns
#' @param x primary site
#' @param fdr fdr cut-off for choosing the top drivers from mutSig2 list of drivers. Default = 0.05
#' @param min_Nmut lower bound of number of cell lines with mutations. Default = 2
#' @param all_cancers_mut_df MAF file from CCLE
#' @param CN_mat copy number matrix from CCLE
#' @param minNrcelllines lower bound of number of cell lines. Default = 5
#' @param meta_data information on different sub types for each primary site
#' @param cellline_annot cell line annotation data frame
#' @param essential_genes vector of essential genes
#' @return An object for each cancer type
#' \describe{
#'    \item{viabilities}{matrix of viabilities for each cancer type}
#'    \item{mutations}{matrix of mutations in drivers for each cancer type}
#'    \item{mutation_annot}{annotations of the mutations}
#'    \item{primary_site}{cancer type}
#' }
#' @export
#' 
#' 
prepareDataObjects <- function(data, x, fdr = 0.05, min_Nmut = 2, all_cancers_mut_df, CN_mat, minNrcelllines = 5, 
                               cellline_annot, meta_data, essential_genes = NULL) { 
  
  # Choosing only cellines for cancer of interest x
  celllines <- cellline_annot[cellline_annot$PRIMARY_SITE == x & cellline_annot$in.CCLE == "Y",]
  
  if(nrow(celllines) < minNrcelllines) {
    return(NULL);
  }else{
    # Merging celllines of different subtypes
    if(subset(meta_data, Primary_site == x, select = Additional_filters, drop = TRUE) == "None"){
      celllines <- celllines$CELLLINE
    }else{
      celllines <- unlist(celllines %>%
                            dplyr::filter(PATHOLOGIST_ANNOTATION %in% 
                                            unlist(strsplit(subset(meta_data,Primary_site==x,select=Additional_filters,drop=TRUE),"[;]"))) %>%
                            dplyr::select(CELLLINE))
    }
    
    # Getting the top driver genes
    top_drivers <- NULL
    
    for(i in unlist(strsplit(subset(meta_data, Primary_site == x, Driver_gene_file, drop = TRUE),"[;]"))){
      # Loading driver mutation list from FireBrowse
      driver_genes <- read.delim(paste(base_folder, "DriverGenes/",i,sep=""),stringsAsFactors = F)
      # Choosing top driver genes with threshold of q <= fdr
      top_drivers <- c(top_drivers, driver_genes $gene[driver_genes $q <= fdr])
    }
    top_drivers <- unique(top_drivers)
    
    # Get the binary CN matrix of cell lines and  top driver genes. 
    CN_subset_int <- prepareCNMat(CN_mat, celllines, top_drivers, x)
    # threshold= 2 no copy number;threshold= 1 - only deep deletions; threshold= 0 - any deletion
    CN_subset_bin <- binarize(2-CN_subset_int, threshold = 2)
    CN_subset_bin[is.na(CN_subset_bin)] <- 0
    
    # Generate mutation matrix of cell line and mutated genes file from each cell lines mutation file
    mut_mat <-  matrix(0,nrow=length(top_drivers),ncol=length(celllines))
    
    for(i in 1:length(celllines)){
      temp_mut <- subset(all_cancers_mut_df,Tumor_Sample_Barcode == toupper(paste(celllines[i],x,sep="_")), Hugo_Symbol)
      mut_idx  <- which(top_drivers %in% unique(temp_mut$Hugo_Symbol))
      mut_mat[mut_idx,i] <- 1 
    }
    
    colnames(mut_mat) <- paste(celllines,x,sep="_")
    rownames(mut_mat) <- top_drivers
    
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
                mutation_annot = mut_mat_annot,
                primary_site = x))
  }
}

#############################
# Loading Cell-lines annotation file from the project DRIVE paper
cellline_annot <-  read.csv(paste0(base_folder, "ProjectDRIVE/TableS2.csv"),stringsAsFactors = F)
# Loading meta datat file for the drivers
meta_data <- read.csv(paste0(base_folder, "MutationFiles/File_metadata.csv"), header = TRUE, stringsAsFactors = FALSE)

data <- readRDS(paste0(base_folder, "ProjectDRIVE/DRIVE_RSA_data.RDS"))

# Remove the essential data from ATARiS normalised data
data <- readRDS(paste0(base_folder, "ProjectDRIVE/DRIVE_ATARiS_data.RDS"))


#add description dependencies
# rerun for ataris drive and validate
