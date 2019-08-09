# Script for identifying true SL pairs from several confounding mutations

library(tableone)
library(ipw)
library(Matching)
library(dplyr)

set.seed(478)
# Loading all the pancancer data
hits_pancan <- read.delim("~/Downloads/Slidr_Results_new/PanCan8pc//Hit_List/SL_hits_pan_cancer.txt",stringsAsFactors = F)
load("~/Downloads/Slidr_Results_new/PanCan8pc/pancan.Rdata")

# creating a data frame for results
causal_hits <- data.frame("driver_gene" = character(),
                          "sl_partner_gene" = character(),
                          "pval" = numeric())

# grouping all the drivers with the same sl partner gene
grouped_hits <- hits_pancan %>%
                  dplyr::group_by(sl_partner_gene) %>%
                  summarise(drivers = paste(driver_gene, collapse = ","))

# set the caliper for matching
set_cal = 0.1

for(i in 1:nrow(grouped_hits)){
  temp_drivers <- sapply(grouped_hits$drivers[i], function(x){strsplit(x,",")[[1]]})
  if(length(temp_drivers) >= 2){ # running for two or more genes
    # defining the causal df
    sl_gene    <- grouped_hits$sl_partner_gene[i]
    causal_df  <- t(pc_data$mutations[temp_drivers,])
    causal_df  <- cbind.data.frame(t(pc_data$viabilities[sl_gene,rownames(causal_df)]), causal_df)
    # To avoid confusion between driver and same gene as target
    colnames(causal_df)[1]  <- "Viability"

    #causal_df[,"Viability"] <- 1 - binarize(causal_df[,"Viability"], threshold = median(causal_df[,"Viability"]))
    colnames(causal_df)     <- gsub('-', "_", colnames(causal_df))
    temp_drivers            <- gsub('-', "_", temp_drivers)

    # For each driver
    for(j in 1: length(temp_drivers)){
      # t_res <- NULL
      t_res <- data.frame("smd_sum" = numeric(),
                          "p_val" = numeric())
      # Since matching depends on the order repeat it 50 times and choose the min pvalue
      for(k in 1:50){
        # Shuffling the rows for matching
        temp_df   <- causal_df[sample(1:nrow(causal_df)),]
        # values of viabilities and treatment var
        y         <- temp_df$Viability
        tr_gene   <- temp_drivers[j]
        # confounding mutations as covariates
        xvars     <- temp_drivers[-j]
        table1    <- CreateTableOne(vars = xvars,
                                    strata = tr_gene,
                                    data = temp_df,
                                    test = FALSE)
        #print(table1,smd=TRUE)

        # Propensity matching
        psmodel <- glm(as.formula(paste0(tr_gene," ~" ,paste(xvars, collapse = " + "))),
                       family = binomial(),
                       data = temp_df)

        #create propensity score
        pscore <- psmodel$fitted.values

        #do greedy matching on logit(PS) using Match with a caliper
        logit   <- function(p) {log(p)-log(1-p)}
        err_mes <- try(Match(Tr = temp_df[,tr_gene],
                             M = 1,
                             X = logit(pscore),
                             replace = FALSE,
                             caliper = set_cal,
                             estimand = "ATT"), silent = FALSE)
        if(is.na(err_mes[[1]])){
          # t_res <- c(t_res,NaN)
          t_res <- rbind.data.frame(t_res,
                                    cbind.data.frame(smd_sum = NaN, p_val = NaN))
        }else{
          psmatch <- Match(Tr = temp_df[,tr_gene],
                           M = 1,
                           X = logit(pscore),
                           replace = FALSE,
                           caliper = set_cal,
                           estimand = "ATT")
          matched <- temp_df[unlist(psmatch[c("index.treated","index.control")]), ]

          #get standardized differences
          matchedtab1 <- CreateTableOne(vars=xvars, strata =tr_gene,
                                        data=matched, test = FALSE)
          #print(matchedtab1, smd = TRUE)
          temp_smd_sum <- sum(ExtractSmd(matchedtab1))

          #outcome analysis
          y_trt <- matched$Viability[matched[,tr_gene] == 1]
          y_con <- matched$Viability[matched[,tr_gene] == 0]

          #pairwise difference
          diff_y <- y_trt - y_con

          # paired t-test
          if(class(try(t.test(diff_y), silent = TRUE)) == "try-error"){
            t_res <- rbind.data.frame(t_res,
                                      cbind.data.frame(smd_sum = temp_smd_sum, p_val = NaN))
          }else
            t_res <- rbind.data.frame(t_res,
                                      cbind.data.frame(smd_sum = temp_smd_sum, p_val = t.test(diff_y)$p.value))
        }
      }

      t_res           <- t_res[!is.nan(t_res$p_val),]

      # if all viabilities are equal then p-value will be NaN and hence we set the value to 1
      if(nrow(t_res) != 0){
        # setting p-val to mean p-val when you have multiple samples with same min smd
        if(sum(t_res$smd_sum == min(t_res$smd_sum)) > 1){
          fin_pval <- t_res %>%
                        dplyr:: filter(smd_sum == min(t_res$smd_sum)) %>%
                        dplyr::summarize(x = mean(p_val))
        }else
          fin_pval <- t_res$p_val[which.min(t_res$smd_sum)]
      }else{
        fin_pval <- Inf # Invalid number for removing entries that could not computed
      }
      causal_hits <- rbind.data.frame(causal_hits,
                                      cbind.data.frame(driver_gene = tr_gene,
                                                       sl_partner_gene = sl_gene,
                                                       pval = as.numeric(fin_pval)))
    }
  }

}

# retaining pairs for which true causal hit exists
best_causal_hits <- causal_hits %>%
                      dplyr::group_by(sl_partner_gene) %>%
                      dplyr::slice(which.min(pval))

sig_causal_hits <- causal_hits %>%
                        dplyr::filter(pval <= 0.05)

# filtering out the pairs for which true causal hit exists
insig_causal_hits <- causal_hits %>%
                        dplyr::filter(pval > 0.05)

# Remove the insignificant causal pairs from the total pancancer hits
filt_hits_pancan <- hits_pancan %>%
                      dplyr::anti_join(insig_causal_hits, by = c("driver_gene", "sl_partner_gene"))

save.image("~/Downloads/Slidr_Results_new/PanCan8pc/causal_res.Rdata")
