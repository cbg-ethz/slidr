#' Function to plot annotated box-plots of SL pairs for a given cancer type
#'
#' @import RColorBrewer
#' @import ggplot2
#' @param canc_data Processed data object for a given cancer type
#' @param hits_df dataframe of SL pairs with their p-values returned by `identifySLHits` function
#' @param path_results The path to where the results should be stored
#' @param WT_pval_thresh Discard SL pairs with WT p-values less than this threshold. Default = 0.2
#' @export
#'
plotSLBoxplot <- function(canc_data, hits_df, path_results, WT_pval_thresh = 0.2){
  # Creating the output folder if it does not exist
  output_folder = paste(path_results,"Plots/", canc_data$primary_site, "/", sep = "")
  if(!dir.exists(output_folder)){
    dir.create(output_folder, recursive = TRUE, showWarnings = TRUE)
  }

  set.seed(1234)

  # Defining the mutation types
  mut_levels <- c("Missense_Mutation", "Splice_Site", "Nonsense_Mutation", "De_novo_Start_Out_Of_Frame",
                  "Start_Codon_SNP", "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                  "In_Frame_Ins", "Stop_Codon_Del", "Stop_Codon_Ins", "Start_Codon_Del", "Start_Codon_Ins", "Normal")

  # Colour code for the mutations
  mut_colors <- c(colorRampPalette(brewer.pal(7,"Set2"))(length(mut_levels)-1),"white")
  names(mut_colors) <- mut_levels

  # Copy number colours
  CNA_colors <- c("#CA1A1A","#807C7C", "#055C95")
  names(CNA_colors) <- c("Amplification","Neutral","Deep Deletion")

  # Filtering out the pairs with significant viabilities in WT
  hits_df <- hits_df %>% dplyr::filter(WT_pvalue > WT_pval_thresh)

  if(nrow(hits_df) > 0){
    viabilities   <- canc_data$viabilities
    mutations     <- canc_data$mutations
    celllines     <- colnames(mutations)
    mut_annot     <- data.frame(canc_data$mutation_annot, stringsAsFactors = F)
    CNalterations <- canc_data$CNalterations

    # Replacing all NAs and NULLs with 1 for neutral
    CNalterations[is.na(CNalterations)]   <- 1
    CNalterations[is.null(CNalterations)] <- 1

    for(j in 1:nrow(hits_df)){
      # Boxplot dataframe for a given SL pair
      temp_data <- rbind.data.frame(mutations[hits_df$driver_gene[j],celllines],
                                    viabilities[hits_df$sl_partner_gene[j],celllines])
      temp_data <- as.data.frame(t(temp_data))
      colnames(temp_data)  <- c("mut_status","viabilities")
      # Mutation status of the cell lines
      temp_data$mut_status <- ifelse(temp_data$mut_status == 0,"WT","Mut")

      # Using a subset of mutation annotation
      temp_annot <- mut_annot[as.character(mut_annot$Hugo_Symbol) == hits_df$driver[j],]
      temp_annot$CN_Value[is.na(temp_annot$CN_Value)]   <- 1 # replacing NA with neutral
      temp_annot$CN_Value[is.null(temp_annot$CN_Value)] <- 1 # replacing NULL with neutral

      # Annotating the different CN alterations for the box plot labels
      temp_data <- cbind(temp_data,
                         do.call(rbind, lapply(1:nrow(temp_data),function(x){
                           if(temp_data$mut_status[x] == "Mut"){# Annotation labels for mut samples
                             df <- temp_annot %>%
                                    dplyr::filter(Cell_Line == rownames(temp_data)[x])

                             y <- strsplit(df$Variant_Classification[[1]],";")[[1]][1]
                             if(df$CN_Value < 1){
                               cbind(y,"Deep Deletion")
                             }else if(df$CN_Value %in% c(1,2) ){
                               cbind(y,"Neutral")
                             }else if(df$CN_Value > 2){
                               cbind(y,"Amplification")
                             }
                           }else{ # Annotation labels for WT samples
                             tempCN <-  CNalterations[hits_df$driver_gene[j],rownames(temp_data)[x]]
                             if(tempCN < 1){
                               cbind("Normal","Deep Deletion")
                             }else if(tempCN %in% c(1,2)){
                               cbind("Normal","Neutral")
                             }else if(tempCN > 2){
                               cbind("Normal","Amplification")
                             }
                           }
                         })))
      colnames(temp_data) <- c("mut_status","viabilities","mut_type","CN_type")

      # Using subset of the colours for the boxplot based on the annotation labels
      subset_CNA_colors <- CNA_colors[levels(temp_data$CN_type)]
      subset_colors     <- mut_colors[levels(temp_data$mut_type)]

      p <- ggplot(temp_data,aes(factor(mut_status),viabilities)) +
           geom_boxplot(fill = "white", color = "grey30",
                        alpha=0.4,
                        size = 0.5,
                        width = 0.3,
                        outlier.shape=NA) +
           geom_jitter(aes(colour = CN_type, fill = mut_type),
                    position = position_jitter(width = 0.2, height = 0.05),
                    shape = 21,
                    size = 2,
                    stroke = 1.5) +
           ylab(paste(hits_df$sl_partner_gene[j],"Viability score")) +
           xlab(hits_df$driver[j])+
           theme_bw()+
           ylim(-4,1.5) +
           theme(panel.grid.major= element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.x = element_line(size = 0.5,color="#525252"),
              axis.text.y=element_text(angle=0, vjust=0.3, hjust=0.5,size=12,colour="#525252"),
              axis.text.x=element_text(angle=0, vjust=0.1, hjust=0,size=12,colour="#525252"),
              axis.title.x=element_text(angle=0, size=14,vjust=0,colour="#525252"),
              axis.title.y=element_text(angle=90, size=14,vjust=0,colour="#525252"),
              panel.border = element_rect(colour = "gray48"),
              legend.position="right",
              strip.text.y=element_blank(),
              strip.text.x=element_text(angle=0,size = 12, colour = "#525252"),
              strip.background = element_blank())+
           scale_color_manual(name = "CNA", values=subset_CNA_colors, labels = names(subset_CNA_colors))+
           scale_fill_manual(name="Mutation type",values=subset_colors, labels = gsub("_"," ", names(subset_colors)))

           # cat(hits_df$driver[j], "\t", hits_df$sl_partner_gene[j], "\n")
           # print(p)
           ggsave(paste(output_folder, hits_df$driver_gene[j], "_", hits_df$sl_partner_gene[j], ".pdf",sep = ""), p, width = 6, height = 5, limitsize = FALSE)
    }
  }
  cat("Finished plotting for", canc_data$primary_site, "data\n")
}

