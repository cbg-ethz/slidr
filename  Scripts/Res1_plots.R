# Load the summary data for the plots
load("~/Downloads/Slidr_Results/Pan_cancer_8pc/summary_data.Rdata")

# Load the hit list
hits_pancan <- read.delim("~/Downloads/Slidr_Results/Pan_cancer_8pc/Hit_List/SL_hits_pan_cancer.txt",stringsAsFactors = F)

# Stacked barplot
barplot_df  <- summary_df %>% dplyr::distinct(driver_gene, n_samples, canc_type, tot_samples)
barplot_df$canc_type <- factor(barplot_df$canc_type,
                               labels = c("Autonomic ganglia", "Bone", "Breast", "CNS", "Endometrium", "Blood", "Kidney", "Large intestine", "Liver", "Lung", "Oesophagus",
                                          "Ovary", "Pancreas", "Pleura", "Skin", "Soft tissue", "Stomach", "UADT", "Urinary tract"))

coloursSites <- c("#CC99BB", "#AA4488", "#771155", "#77AADD",
                  "#4477AA", "#114477", "#77CCCC", "#44AAAA", "#117777", "#88CCAA",
                  "#44AA77", "#117744", "#DDDD77", "#AAAA44", "#777711", "#DDAA77",
                  "#AA7744", "#774411", "#DD7788")

c <- ggplot(barplot_df, aes(fill = canc_type, x = reorder(driver_gene,-tot_samples), y = n_samples)) +
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
  scale_fill_manual(values = coloursSites)+
  scale_y_continuous(breaks = seq(0, 280, by = 40))

# Heatmap of important hits
# top_sl_pairs    <- paste(hits_pancan$driver_gene[1:100], hits_pancan$sl_partner_gene[1:100], sep = "_")
top_drivers   <- c("ACVR2A", "APC", "BRAF", "COL5A1","CTNNB1", "DMRTA1", "KRAS", "MTAP", "NOTCH1", "SYNE1") #"RPL22", COL24A1
hm_summary_df <- summary_df %>% filter(driver_gene %in% top_drivers)

d <- ggplot(hm_summary_df, aes(x = sl_partner_gene, y = canc_type)) +
  geom_tile(aes(fill = -log(as.numeric(as.character(mut_pvalue))))) +
  facet_grid(.~ factor(driver_gene), scales = "free_x", space = "free_x") +
  scale_fill_viridis(option="D", begin = 0.1, end = 0.9, alpha = 0.9) +
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
        strip.text.x = element_text(size = 9, colour = "#525252", face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(vjust = 1, size=12,colour="#525252"),
        legend.text = element_text(size=10,colour="#525252"),
        legend.key.size = unit(0.6, "cm"))+
  scale_y_discrete(labels = c("Autonomic ganglia", "Bone", "Breast", "CNS", "Endometrium",
                              "Blood", "Kidney", "Large intestine", "Liver", "Lung", "Oesophagus",
                              "Ovary", "Pancreas", "Pleura", "Skin", "Soft tissue", "Stomach", "UADT", "Urinary tract"))

# Confounders plot
confounders   <- c("ADSL","MAT2A","PRMT5","RPA3")
bp_summary_df <- hits_pancan %>% filter(sl_partner_gene %in% confounders)
e <- ggplot(bp_summary_df, aes(x=driver_gene, y=sl_partner_gene, size = -log(mut_pvalue, base = 10))) +
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

al_plots <- align_plots(c,d,e, align = "v", axis = "l" )
row_1 <- plot_grid(NULL, NULL, labels = c("A","B"), nrow = 1, rel_widths = c(0.6,0.4), label_size = 13)
row_2 <- plot_grid(al_plots[[1]], labels = c("C"),  nrow = 1, label_size = 13)
row_3 <- plot_grid(al_plots[[2]], labels = c("D"),  nrow = 1, label_size = 13)
row_4 <- plot_grid(al_plots[[3]],NULL, NULL, labels = c("E","","F"), nrow = 1, rel_widths = c(0.6,0.02,0.38),label_size = 13)
fin_plot <- plot_grid(row_1, row_2, row_3, row_4,
                      nrow = 4,
                      ncol = 1,
                      rel_heights = c(0.85,1.1,1, 0.65),
                      align = "v",
                      axis = 'l')

ggsave(fin_plot, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/finalPlot_res1",Sys.Date(),".pdf"),
       width = 16, height = 22)

