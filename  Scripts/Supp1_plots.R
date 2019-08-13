# Plotting supplementary figure for pan-cancer hits

library(ggplot2)
library(dplyr)
library(cowplot)

load("~/Downloads/Slidr_Results_new/PanCan8pc/causal_res.Rdata")

# Confounders plot for the top causally inferred hits
confounders   <- setdiff(sig_causal_hits$sl_partner_gene, c("PRMT5","MAT2A","RPL22L1", "CTNNB1"))
bp_summary_df <- hits_pancan %>% filter(sl_partner_gene %in% confounders)
# removing long names for better aesthetics
bp_summary_df$sl_partner_gene <- unlist(lapply(bp_summary_df$sl_partner_gene,
                                               function(x){strsplit(x, ",")[[1]][1]}))

bp_summary_df  <- bp_summary_df %>% arrange(sl_partner_gene)

a <- ggplot(bp_summary_df, aes(x= factor(driver_gene, levels = unique(driver_gene)), y=sl_partner_gene, size = -log10(mut_pvalue))) +
  geom_point(fill = "#055C95",color="#055C95", alpha=0.45) +
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() +
  xlab("Driver genes") +
  ylab("SL partner genes") +
  labs(size = "-log10(p-value)")+
  theme(panel.grid.major = element_line(size = 0.1),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size=0.5,color="#525252"),
        axis.line = element_line(colour="#525252"),
        axis.text.y=element_text(angle=0, vjust=0.5, hjust=1,size=11,colour="#525252"),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=11,colour="#525252"),
        axis.title.x=element_text(angle=0,size=13,face="bold",vjust=0,colour="#525252"),
        axis.title.y=element_text(angle=90, size=13,face="bold",vjust=0,colour="#525252"),
        strip.text.x = element_text(size = 9, colour = "#525252", face = "bold"),
        strip.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(vjust = 0.5, size=11,colour="#525252"),
        legend.text = element_text(size=10.5,colour="#525252"),
        legend.key.size = unit(0.25, "cm"))

row_1 <- plot_grid(a, NULL,labels = c("A", "B"),  rel_widths = c(0.7,0.3), nrow = 1, label_size = 13)
row_2 <- plot_grid(NULL, labels = c("C"),  nrow = 1, label_size = 13)
fin_plot <- plot_grid(row_1, row_2,
                      nrow = 2,
                      ncol = 1,
                      rel_heights = c(1,0.4),
                      # align = "v",
                      axis = 'l')

ggsave(fin_plot, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/finalPlot_supp1",Sys.Date(),".pdf"),
       width = 14, height = 11)

# Hits for the Venn diagram of FP and F
rm(list = ls())

# Loading cancer specific hits
load("~/Downloads/Slidr_Results_new/ContCN/ProcessedData.Rdata")

path_results   <- "~/Downloads/Slidr_Results_new/"
# Testing all pairs
# all_hits_liver <- slidr::identifySLHits(canc_data = all_data$liver, path_results = path_results, WT_pval_thresh = 0, qval_thresh = Inf)
all_hits_liver <- read.delim("~/Downloads/Slidr_Results_new/Hit_List/SL_hits_liver.txt", stringsAsFactors = FALSE)
all_hits_liver$fdr <- p.adjust(all_hits_liver$mut_pvalue, "fdr")
all_hits_liver <- all_hits_liver %>% dplyr::filter(WT_pvalue >= 0.1 & fdr <= 0.1)

# Hits from our method
sub_hits_liver <- read.delim("~/Downloads/Slidr_Results_new/ContCN/Hit_List/SL_hits_liver.txt", stringsAsFactors = FALSE)
sub_hits_liver <- sub_hits_liver %>% dplyr::filter(WT_pvalue >= 0.1)
