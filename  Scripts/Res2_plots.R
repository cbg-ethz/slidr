# Plotting figure for cancer specific hits.

library(ggplot2)
library(dplyr)
library(cowplot)
library(viridis)
library(circlize)
library(ggplotify)

load("~/Downloads/Slidr_Results_new/ContCN/ProcessedData.Rdata")

# # Reading the hit list for 17 cancers
# hit_files   <- list.files("~/Downloads/Slidr_Results/ContCN/Hit_List/", full.names = TRUE)
# hits        <- lapply(hit_files, function(x){read.delim(x, stringsAsFactors = FALSE, sep = "\t")})
# names(hits) <- lapply(hit_files, function(x){sub("SL_hits_", "",strsplit(x, "\\/|\\.")[[1]][9])})
hits        <- lapply(hits, function(x){ x %>% dplyr::filter(WT_pvalue >= 0.1)})
hits        <- hits[Primary_sites]

# A vector of all the drivers
all_drivers  <- unique(unlist(lapply(hits, function(x){x$driver_gene})))
all_partners <- unique(unlist(lapply(hits, function(x){x$sl_partner_gene})))

# Create a data-frame of mutation frequencies of each driver in each cancer
mut_freq_df <- do.call(rbind.data.frame, sapply(all_data,
                                                function(x){
                                                  cbind(rep(x$primary_site, nrow(x$mutations)),
                                                        rownames(x$mutations),
                                                        rowMeans(x$mutations))}))
colnames(mut_freq_df) <- c("canc_type", "driver_gene", "freq")

# filtering out drivers without an SL partner
mut_freq_df <- mut_freq_df %>%
                dplyr::filter(driver_gene %in% all_drivers)
#mut_freq_df$canc_type <- factor(mut_freq_df$canc_type, levels = sort(levels(mut_freq_df$canc_type)))

#Filtering out drivers specific only to 1 cancer type
mut_freq_df <- mut_freq_df %>%
                dplyr::group_by(driver_gene) %>%
                dplyr::filter(n() > 1)

# Plotting driver frequencies across cancer types
a <- ggplot(mut_freq_df, aes(x = driver_gene, y = factor(canc_type)))+
      geom_tile(aes(fill = as.numeric(as.character(freq))), color = "white") +
      scale_fill_viridis(option="A", begin = 0, end = 0.95, alpha = 0.9, direction = -1) +
      theme_bw() +
      ylab("Primary sites") +
      xlab("Driver genes") +
      labs(fill = "Frequency")+
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
            legend.key.size = unit(0.6, "cm")) +
  scale_y_discrete(labels = c("Pancreas", "Lung", "Liver", "Breast", "Skin", "Endometrium", "Thyroid", "Stomach",
                              "Bone", "Ovary", "Blood", "Kidney", "CNS", "Oesophagus", "Urinary tract", "Large intestine"))


# Chord diagram
cs_hits <- read.delim("~/Downloads/Slidr_Results_new/CanSpecific_literature.txt", stringsAsFactors = FALSE)
cs_hits$sl_partner_gene <- sapply(cs_hits$sl_partner_gene, function(x){strsplit(x, ",")[[1]][1]})
cs_mat <- reshape2::acast(cs_hits, driver_gene~sl_partner_gene, value.var="mut_qvalue")
cs_mat[is.na(cs_mat)] <- 0
cs_mat[cs_mat>0]      <- 1
# colnames(cs_mat) <- sapply(colnames(cs_mat), function(x){strsplit(x, ",")[[1]][1]})

color_sites <- c("#CC99BB", "#771155", "#77AADD", "#44AA77",
                  "#114477", "#DD7788", "#44AAAA", "#774411",
                  "#DDDD77", "#147C47", "#DDAA77")

color_sites <- viridis(length(unique(cs_hits$Type)), option = "D", alpha = 0.85, begin = 0.1, end = 0.96)
names(color_sites) <- sort(unique(cs_hits$Type))

color_df <- cbind.data.frame(cs_hits$driver_gene, cs_hits$sl_partner_gene, color_sites[cs_hits$Type])

circos.clear()
circos.par(start.degree = 270)
# Basic chord diagram without the labels
chordDiagram(cs_mat,
             annotationTrack = "grid",
             annotationTrackHeight = 0.03,
             preAllocateTracks = 1,
             grid.col = "#525252",
             grid.border = "#525252",
             transparency = 0.5,
             col = color_df)
# Adding vertical labels track
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y){
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim),
              ylim[1] + .1,
              sector.name,
              facing = "clockwise",
              niceFacing = TRUE,
              adj = c(0.05, 0.5),
              cex = 0.8,
              col = "#525252")
  # circos.axis(h = "top",
  #             labels.cex = 1e-7,
  #             major.tick = FALSE,
  #             sector.index = sector.name,
  #             track.index = 2)
}, bg.border = NA)

# Converting it into gtable for using plot_grid
b <- plot_to_gtable(recordPlot())
b <- as.ggplot(b)
b <- b +
      theme_nothing() +
      theme(aspect.ratio = 1)

# Adding legend: Plotting a density plot to get the legend
lgd <- get_legend(ggplot(cs_hits, aes(mut_qvalue, color = Type, fill = Type)) +
                    geom_density(alpha = .85) +
                    theme(legend.position = "right",
                          legend.title = element_text(vjust = 1, size=12,colour="#525252"),
                          legend.text = element_text(size=10,colour="#525252"),
                          legend.key.size = unit(0.5, "cm"),
                          legend.key = element_rect(color="white")) +
                    scale_fill_manual(values = color_sites, name = "Primary site")+
                    scale_color_manual(values = color_sites, name = "Primary site"))

row_1 <- plot_grid(a, labels = c("A"),  nrow = 1, label_size = 13)
row_2 <- plot_grid(b, lgd, NULL, labels = c("B","","C"), rel_widths = c(0.3,0.005,0.2), nrow = 1, label_size = 13)
fin_plot <- plot_grid(row_1, row_2,
                      nrow = 2,
                      ncol = 1,
                      rel_heights = c(1,1.2))

ggsave(fin_plot, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/finalPlot_res2",Sys.Date(),".pdf"),
       width = 12, height = 12.5)
