# Script for plotting top PRISM hits
library(dplyr)
library(stringi)
library(ggplot2)

# Read the files
hit_ls <- list.files("~/Downloads/Slidr_Results_new/DepMap/Prism/30Nov2021_FDR/", full.names = TRUE, recursive = TRUE)

# Extracting the cancer names for the plot
cancer_names <- list.files("~/Downloads/Slidr_Results_new/DepMap/Prism/30Nov2021_FDR/", full.names = FALSE, recursive = TRUE)
cancer_names <- gsub("prism_validation_", "", cancer_names)
cancer_names <- gsub(".txt", "", cancer_names)
cancer_names <- gsub("_", " ", cancer_names) %>% stri_trans_totitle()

# Dataframes for plotting
df_bar <- data.frame()

for(i in 1:length(cancer_names)){
  temp_df <- read.delim(hit_ls[i], header = TRUE, stringsAsFactors = FALSE)
  if(nrow(temp_df) > 0){
    # Updating dataframe for the bar chart
    temp <- bind_rows(bind_cols(cancer_names[i], "FDR", sum(temp_df$FDR <= 0.2), nrow(temp_df)),
                      bind_cols(cancer_names[i], "Hits", nrow(temp_df) - sum(temp_df$FDR <= 0.2), nrow(temp_df)))
    df_bar <- bind_rows(df_bar, temp)
  }
}

colnames(df_bar) <- c("Primary_Site", "Significance", "value", "tot_hits")

# Plotting stacked bar plot
p <- ggplot(df_bar, aes(fill = Significance, x = reorder(Primary_Site, -tot_hits) , y = value)) +
      geom_bar(stat = "identity", width=0.7) +
      ylim(0,50) +
      theme_bw() +
      ylab("Number of significant candidate drugs") +
      xlab("Primary Site") +
      #labs(fill = "Primary site")+
      theme(panel.grid= element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size=0.5,color="#525252"),
        axis.line = element_line(colour="#525252"),
        axis.text.y=element_text(angle=0, vjust=0.5, hjust=1,size=7.5,colour="#525252"),
        axis.text.x=element_text(angle=45, hjust=1,size=7.5,colour="#525252"),
        axis.title.x=element_text(angle=0,size=9,face="bold",vjust=5,colour="#525252"),
        axis.title.y=element_text(angle=90, size=9,face="bold",vjust=1.5,colour="#525252"),
        legend.position = "none") +
      scale_fill_manual(values = c("#5a92b0","#9ecae1"))

ggsave("~/Downloads/Slidr_Results_new/DepMap/Prism/FDR_Barplot.pdf", p, width = 6, height = 4, units = "in")
