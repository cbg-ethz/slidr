library(msigdbr)
library(GOstats)
library(hgu95av2.db)
library(org.Hs.eg.db)
library(dplyr)

load("~/Downloads/Slidr_Results/ContCN/ProcessedData.Rdata")

# Reading the hit list for 17 cancers
hit_files <- list.files("~/Downloads/Slidr_Results/ContCN/Hit_List/", full.names = TRUE)
hits        <- lapply(hit_files, function(x){read.delim(x, stringsAsFactors = FALSE, sep = "\t")})
names(hits) <- lapply(hit_files, function(x){sub("SL_hits_", "",strsplit(x, "\\/|\\.")[[1]][9])})
hits        <- lapply(hits, function(x){ x %>% dplyr::filter(WT_pvalue >= 0.2)})
hits        <- hits[Primary_sites]

# curated KEGG pathways from msigdb
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
m_df <- m_df[grep("SIGNAL",m_df$gs_name),]

# Defining the universe
universe <- as.character(unique(m_df$entrez_gene))

# Function for running hypergeometric test and correct for multiple testing
HG_test <- function(hit){

  # selected genes
  selectedGenes <- unlist(lapply(unique(hit$sl_partner_gene), function(x){unlist(strsplit(x, ",")[[1]])}))
  geneMap       <- select(org.Hs.eg.db, selectedGenes, "ENTREZID", "SYMBOL")
  # When symbol maps to many entrez remove multiple maps
  selectedGenes <- geneMap[!duplicated(geneMap$SYMBOL), "ENTREZID"]

  # if intersection between selected genes and universe is not null
  if(any(selectedGenes %in% universe)){
    # Defining the params for hypergeometric test
    params <- new("GOHyperGParams",
                  geneIds = selectedGenes,
                  universeGeneIds = universe,
                  annotation="hgu95av2.db",
                  ontology = "BP",
                  pvalueCutoff = 1,
                  testDirection = "over")

    hgOver <- hyperGTest(params)

    # correcting for multiple testing
    res <- summary(hgOver) %>% tbl_df() %>%
      dplyr::mutate(Pvalue = p.adjust(Pvalue, method = "BH")) %>%
      dplyr::filter(Pvalue < 0.1) %>%
      dplyr::arrange(Pvalue)

    return(res)
  }else{
    return(tibble(a = character(), b = character()))
  }
}

# Running GO on all the hits
GO_res  <- mclapply(hits, HG_test, mc.cores = 2)
# Removing cancers with no enrichment
GO_res  <- lapply(GO_res, function(x){if(nrow(x) > 0){x}})
GO_res  <- GO_res[-which(sapply(GO_res, is.null))]

# Choosing only signalling pathways
sig_pathways <- unlist(lapply(GO_res, function(x){x$Term})) %>% unique()
sig_pathways <- sig_pathways[grep("signal", sig_pathways)]

# create a data frame with all signalling pathways
cs_pathways <- do.call(rbind.data.frame, lapply(1:length(GO_res), function(x){
                  temp <- GO_res[[x]] %>%
                            dplyr::filter(Term %in% sig_pathways)
                  rows <- nrow(temp)
                  if(rows > 0){
                    cbind.data.frame(rep(names(GO_res)[x], rows),
                                      temp)
                  }
}))
