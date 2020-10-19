library(dplyr)

# Load pan cancer object
load("~/Downloads/Slidr_Results_new/PanCan8pc/pancan.Rdata")
rm(list=setdiff(ls(), c("all_data", "hits", "pc_data")))

# Load hits
pc_hits  <- read.delim2("~/Downloads/Slidr_Results_new/PanCan8pc/Hit_List/SL_hits_pan_cancer.txt",
                       stringsAsFactors = FALSE, header = TRUE)
# Filter out self-dependent hits
pc_hits <- dplyr::filter(pc_hits, driver_gene != sl_partner_gene)
# Load the primary drug screens data
drugs_df <- read.csv("~/Downloads/Slidr_Results_new/DepMap/primary-screen-replicate-collapsed-treatment-info.csv",
                     stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
drug_screen <- read.csv("~/Downloads/Slidr_Results_new/DepMap/primary-screen-replicate-collapsed-logfold-change.csv",
                        header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(drug_screen)[1] <- "DepMap_ID"
# Removing the 10 FAILED_STR entries
drug_screen <- drug_screen[-grep("FAILED_STR", drug_screen$DepMap_ID),]

sample_info <- read.csv("~/Downloads/Slidr_Results_new/DepMap/sample_info.csv",
                        header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# Converting DepMAp_ID to cell names for rownames
drug_screen <- drug_screen %>%
                  left_join(select(sample_info, DepMap_ID, CCLE_Name), by = "DepMap_ID")

rownames(drug_screen) <- tolower(drug_screen$CCLE_Name)

# fewer columns selected for downstream
drugs_red <- drugs_df %>%dplyr::select(column_name, broad_id, name, target)

# retrieving the subset of drugs for the SL partner genes
subset_drugs <- do.call(rbind.data.frame,
                        sapply(unique(pc_hits$sl_partner_gene),
                                     function(x){drugs_red[which(grepl(x, drugs_red$target)),]}, simplify = FALSE))

subset_drugs$sl_partner_gene <- sapply(row.names(subset_drugs), function(x){strsplit(x,"[.]")[[1]][1]})

# initializing variables for the results
results   <- data.frame()
mutations <- pc_data$mutations

for(i in 1:nrow(pc_hits)){
  target_gene <- pc_hits$sl_partner_gene[i]
  driver_gene <- pc_hits$driver_gene[i]

  if(target_gene %in% subset_drugs$sl_partner_gene){
    temp <- dplyr::filter(subset_drugs, target_gene == sl_partner_gene)

    # Stratifying based on driver gene mutation status
    WT_celllines  <- colnames(mutations)[which(mutations[driver_gene,] == 0)]
    mut_celllines <- colnames(mutations)[which(mutations[driver_gene,] == 1)]

    # For multiple compounds targeting the same partner gene
    for(j in 1:nrow(temp)){
      compound <- temp$column_name[j]

      # viabilities in mutated and WT cell-lines from the drug
      mut_viability <- drug_screen[mut_celllines,compound]
      WT_viability  <- drug_screen[WT_celllines,compound]

      # Testing if the viability in mutated cell lines is significantly less than in WT from the drug
      p_val <- t.test(mut_viability, WT_viability, alternative =  "less")$p.value

      results <- rbind.data.frame(results,
                                  cbind(pc_hits[i,], temp[j,], p_val))
    }
  }
}
