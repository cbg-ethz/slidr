library(dplyr)

########### Load slidr pred data ###############
# Load pan cancer object
load("~/Downloads/Slidr_Results_new/PanCan8pc/pancan.Rdata")
rm(list=setdiff(ls(), c("all_data", "hits", "pc_data")))

# Load hits
pc_hits  <- read.delim2("~/Downloads/Slidr_Results_new/PanCan8pc/Hit_List/SL_hits_pan_cancer.txt",
                        stringsAsFactors = FALSE, header = TRUE)

########### Prism screen data ###############
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

########### Validation on prism ###############

# Function to test and validate drug targets for the predicted SL pairs on prism data
getSigDrugs <- function(canc_data, pred_hits){

  # Filter out self-dependent hits
  pred_hits <- dplyr::filter(pred_hits, driver_gene != sl_partner_gene)

  # retrieving the subset of drugs for the SL partner genes
  subset_drugs <- do.call(rbind.data.frame,
                          sapply(unique(pred_hits$sl_partner_gene),
                                 function(x){drugs_red[which(grepl(paste0('\\b',x,'\\b'), drugs_red$target)),]}, simplify = FALSE))

  subset_drugs$sl_target_gene <- sapply(row.names(subset_drugs), function(x){strsplit(x,"[.]")[[1]][1]})

  # initializing variables for the results
  results   <- data.frame()
  mutations <- canc_data$mutations

  for(i in 1:nrow(pred_hits)){
    target_gene <- pred_hits$sl_partner_gene[i]
    driver_gene <- pred_hits$driver_gene[i]

    if(target_gene %in% subset_drugs$sl_target_gene){
      temp <- dplyr::filter(subset_drugs, target_gene == sl_target_gene)

      # Stratifying based on driver gene mutation status
      WT_celllines  <- colnames(mutations)[which(mutations[driver_gene,] == 0)]
      mut_celllines <- colnames(mutations)[which(mutations[driver_gene,] == 1)]

      # For multiple compounds targeting the same partner gene
      for(j in 1:nrow(temp)){
        compound <- temp$column_name[j]

        # viabilities in mutated and WT cell-lines from the drug
        mut_viability <- drug_screen[mut_celllines,compound]
        WT_viability  <- drug_screen[WT_celllines,compound]

        # Ensuring that the missing data is less than half
        if(sum(is.na(mut_viability)) < length(mut_celllines)/2 && sum(is.na(WT_viability)) < length(WT_celllines)/2){
          # Testing if the viability in mutated cell lines is significantly less than in WT from the drug
          p_val <- t.test(mut_viability, WT_viability, alternative =  "less")$p.value

          results <- rbind.data.frame(results,
                                      cbind(pred_hits[i,], temp[j,], p_val))
        }
      }
    }
  }
  return(list(avail_drugs = subset_drugs,
              valid_hits = results))
}

########### run pan-cancer ###############
# Validate hits for pancancer data
pc_valid <- getSigDrugs(pc_data, pc_hits)

res_df <- pc_valid$valid_hits %>%
            dplyr::filter(p_val < 0.1) %>%
            dplyr::arrange(p_val) %>%
            dplyr::select(driver_gene, sl_partner_gene,
                          WT_pvalue, mut_pvalue,
                          broad_id, name, target, p_val)

colnames(res_df) <- c("Driver gene","SL partner gene",
                      "WT pvalue",	"mut pvalue", "Broad ID",
                      "Drug name", "Drug target(s)", "Validation pvalue")

write.table(x = res_df,
            file = paste0("~/Downloads/Slidr_Results_new/DepMap/Prism/prism_validation_pancancer.txt"),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

########### run cancer type-specific ###############
# Validate hits for cancer-type specific
cs_hits <- lapply(hits, function(x){ x %>% dplyr::filter(WT_pvalue >= 0.1)})
cs_valid <- apply(mapply(getSigDrugs, all_data, cs_hits), 2, as.list)

for(tissue in names(cs_valid)){

  res_df <- cs_valid[[tissue]]$valid_hits

  if(nrow(res_df) > 0){
    res_df <- res_df %>%
                dplyr::filter(p_val < 0.1) %>%
                dplyr::arrange(p_val) %>%
                dplyr::select(driver_gene, sl_partner_gene,
                              WT_pvalue, mut_pvalue,
                              broad_id, name, target, p_val)

    colnames(res_df) <- c("Driver gene","SL partner gene",
                          "WT pvalue",	"mut pvalue", "Broad ID",
                          "Drug name", "Drug target(s)", "Validation pvalue")
    write.table(x = res_df,
                file = paste0("~/Downloads/Slidr_Results_new/DepMap/Prism/prism_validation_", tissue, ".txt"),
                quote = FALSE,
                row.names = FALSE,
                sep = "\t")
  }

}

save.image("~/Downloads/Slidr_Results_new/DepMap/Prism/prism_validation.Rdata")
