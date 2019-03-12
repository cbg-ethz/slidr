# Plotting figure for pan-cancer hits

library(cowplot)
library(ggplot2)
library(viridis)
library(dplyr)

# Load the summary data for the plots
# load("~/Downloads/Slidr_Results/Pan_cancer_8pc/summary_data.Rdata")
load("~/Downloads/Slidr_Results/Pan_cancer_8pc/summary_data.Rdata")

# Load the hit list
hits_pancan <- read.delim("~/Downloads/Slidr_Results/Pan_cancer_8pc/Hit_List/SL_hits_pan_cancer.txt",stringsAsFactors = F)

# Stacked barplot
barplot_df  <- summary_df %>% dplyr::distinct(driver_gene, n_samples, canc_type, tot_samples)
barplot_df$canc_type <- factor(barplot_df$canc_type,
                               labels = c("Autonomic ganglia", "Bone", "Breast", "CNS", "Endometrium", "Blood", "Kidney", "Large intestine", "Liver", "Lung", "Oesophagus",
                                          "Ovary", "Pancreas", "Pleura", "Skin", "Soft tissue", "Stomach", "UADT", "Urinary tract"))

colour_sites <- c("#CC99BB", "#AA4488", "#771155", "#77AADD",
                  "#4477AA", "#114477", "#77CCCC", "#44AAAA", "#117777", "#88CCAA",
                  "#44AA77", "#117744", "#DDDD77", "#AAAA44", "#777711", "#DDAA77",
                  "#AA7744", "#774411", "#DD7788")

a <- ggplot(barplot_df, aes(fill = canc_type, x = reorder(driver_gene,-tot_samples), y = n_samples)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Number of mutated cell lines") +
  xlab("Driver genes") +
  labs(fill = "Primary site")+
  theme(panel.grid= element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size=0.5,color="#525252"),
        axis.line = element_line(colour="#525252"),
        axis.text.y=element_text(angle=0, vjust=0.5, hjust=1,size=10,colour="#525252"),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=10,colour="#525252"),
        axis.title.x=element_text(angle=0,size=12,face="bold",vjust=0,colour="#525252"),
        axis.title.y=element_text(angle=90, size=12,face="bold",vjust=-15,colour="#525252"),
        strip.text.x = element_text(size = 10, colour = "#525252", face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(vjust = 1, size=12,colour="#525252"),
        legend.text = element_text(size=10,colour="#525252"),
        legend.key.size = unit(0.5, "cm"))+
  scale_fill_manual(values = colour_sites)+
  scale_y_continuous(breaks = seq(0, 280, by = 40))

# Loading the causal results for plotting confounders and heatmap
load("~/Downloads/Slidr_Results/Pan_cancer_8pc/causal_res.RData")
# Loading the summary dataframe after removing confounders
load("~/Downloads/Slidr_Results/Pan_cancer_8pc/filt_summary_data.Rdata")

# Confounders plot for the top causally inferred hits
confounders   <- best_causal_hits$sl_partner_gene #c("ADSL","MAT2A","PRMT5","RPA3")
bp_summary_df <- hits_pancan %>% filter(sl_partner_gene %in% confounders)
# removing long names for better aesthetics
bp_summary_df$sl_partner_gene <- unlist(lapply(bp_summary_df$sl_partner_gene,
                                        function(x){strsplit(x, ",")[[1]][1]}))
b <- ggplot(bp_summary_df, aes(x=driver_gene, y=sl_partner_gene, size = -log(mut_pvalue, base = 10))) +
  geom_point(fill = "#055C95",color="#055C95", alpha=0.45) +
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() +
  xlab("Driver genes") +
  ylab("SL partner genes") +
  labs(size = "-log(p-value)")+
  theme(panel.grid.major = element_line(size = 0.1),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size=0.5,color="#525252"),
        axis.line = element_line(colour="#525252"),
        axis.text.y=element_text(angle=0, vjust=0.5, hjust=1,size=10,colour="#525252"),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=10,colour="#525252"),
        axis.title.x=element_text(angle=0,size=12,face="bold",vjust=0,colour="#525252"),
        axis.title.y=element_text(angle=90, size=12,face="bold",vjust=-12,colour="#525252"),
        strip.text.x = element_text(size = 9, colour = "#525252", face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(vjust = 0.5, size=12,colour="#525252"),
        legend.text = element_text(size=10,colour="#525252"),
        legend.key.size = unit(0.6, "cm"))

# Heatmap of important hits after causal filtering step
# removing long names for better aesthetics
summary_df$sl_partner_gene <- unlist(lapply(summary_df$sl_partner_gene,
                                               function(x){strsplit(x, ",")[[1]][1]}))

d <- ggplot(summary_df, aes(x = sl_partner_gene, y = canc_type)) +
  geom_tile(aes(fill = -log(as.numeric(as.character(mut_pvalue))))) +
  facet_grid(.~ factor(driver_gene), scales = "free_x", space = "free_x") +
  scale_fill_viridis(option="D", begin = 0, end = 1, alpha = 0.9) +
  theme_bw() +
  ylab("Primary sites") +
  xlab("SL partner genes") +
  labs(fill = "-log(p-value)")+
  theme(panel.grid= element_blank(),
        axis.ticks = element_line(size=0.5,color="#525252"),
        axis.line = element_line(colour="#525252"),
        axis.text.y=element_text(angle=0, vjust=0.5, hjust=1,size=10,colour="#525252"),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1,size=10,colour="#525252"),
        axis.title.x=element_text(angle=0,size=12,face="bold",vjust=0,colour="#525252"),
        axis.title.y=element_text(angle=90, size=12,face="bold",vjust=-2,colour="#525252"),
        strip.text.x = element_text(size = 9, angle = 90, colour = "#525252", face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(vjust = 1, size=12,colour="#525252"),
        legend.text = element_text(size=10,colour="#525252"),
        legend.key.size = unit(0.6, "cm"))+
  scale_y_discrete(labels = c("Autonomic ganglia", "Bone", "Breast", "CNS", "Endometrium",
                              "Blood", "Kidney", "Large intestine", "Liver", "Lung", "Oesophagus",
                              "Ovary", "Pancreas", "Skin", "Soft tissue", "Stomach", "UADT", "Urinary tract")) # No Pleura after causal filtering

# aligning the y axis of all plots
al_plots <- align_plots(a,b,d, align = "v", axis = "l" )
row_1 <- plot_grid(al_plots[[1]], labels = c("A"),  nrow = 1, label_size = 13)
row_2 <- plot_grid(al_plots[[2]], NULL,labels = c("B", "C"),  rel_widths = c(0.75,0.25), nrow = 1, label_size = 13)
row_3 <- plot_grid(al_plots[[3]], labels = c("D"),  nrow = 1, label_size = 13)
fin_plot <- plot_grid(row_1, row_2, row_3,
                      nrow = 3,
                      ncol = 1,
                      rel_heights = c(1,0.85,1),
                      align = "v",
                      axis = 'l')

ggsave(fin_plot, filename = paste0("/Volumes/bsse_group_beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/finalPlot_res1",Sys.Date(),".pdf"),
       width = 20, height = 20)

