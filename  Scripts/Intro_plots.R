library(viridis)
library(ggplot2)
library(gridBase)
library(dplyr)
library(cowplot)

# Circular bar plot

load("/Volumes//beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Plots/2018/13Aug2018/Pan_cancer_8pc/pancan.Rdata")

sites_vec   <- sapply(pc_celllines, function(x){paste(strsplit(x, "_")[[1]][-1], collapse = " ")})
sites_tab   <- table(sites_vec)
sites_tab   <- sites_tab[sites_tab > 2]
sites_df    <- as.data.frame(sites_tab)
colnames(sites_df) <- c("individual", "value")
sites_df$id <- 1:nrow(sites_df)


label_data <- sites_df
label_data$sites <- c("Autonomic ganglia", "Bone", "Breast", "CNS", "Endometrium",
                      "Blood", "Kidney", "Large intestine", "Liver", "Lung", "Oesophagus",
                      "Ovary", "Pancreas", "Pleura", "Skin", "Soft tissue", "Stomach", "Thyroid",
                      "UADT", "Urinary tract")

# calculate the ANGLE of the labels
number_of_bar = nrow(label_data)
angle = 90 - 360 * (label_data$id - 0.5)/number_of_bar

# calculate the alignment of labels: right or left
label_data$hjust <- ifelse( angle < -90, 1, 0)
# flip angle BY to make them readable
label_data$angle <- ifelse(angle < -90, angle+180, angle)

b <- ggplot(sites_df, aes(x=individual, y=value)) +
  geom_bar(stat="identity", fill=alpha("#055C95", 0.6)) +
  ylim(-30,80) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")) +
  coord_polar(start = 0) +
  geom_text(data = label_data,
            aes(x = id,
                y = value + 12,
                label = sites,
                hjust = hjust),
            color = "#525252",
            alpha = 0.8,
            fontface = "bold",
            size = 3.5,
            angle = label_data$angle,
            inherit.aes = FALSE) +
  geom_text(data = label_data,
            aes(x = id,
                y = value + 2 ,
                label = value,
                hjust = hjust),
            color = "#525252",
            alpha = 0.8,
            #fontface = "bold",
            size = 3.5,
            angle = label_data$angle,
            inherit.aes = FALSE)

ggsave(b, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/primarySites",Sys.Date(),".pdf"),
       width = 5.5, height = 5.5)


rm(list = setdiff(ls(), "b"))
set.seed(496)

n_celllines <- 10
n_perturbs  <- 4000
n_mutations <- 4

celllines   <- paste(rep("CL", n_celllines), 1:n_celllines, sep = " ")
perturbs    <- paste(rep("PG", n_perturbs), 1:n_perturbs, sep = " ")
mutations   <- paste(rep("MG", n_mutations), 1:n_mutations, sep = " ")

mut_mat                      <- matrix(0,n_mutations, n_celllines)
dimnames(mut_mat)            <- list(mutations,celllines)
mut_mat["MG 1",c("CL 2", "CL 3", "CL 4", "CL 5")] <- 1
mut_mat["MG 2",c("CL 7", "CL 9", "CL 10")]        <- 1
mut_mat["MG 3",c("CL 1", "CL 6", "CL 5", "CL 8")] <- 1
mut_mat["MG 4",c("CL 1", "CL 3", "CL 10")]        <- 1
mut_df                       <- melt(mut_mat)
mut_df$bordercol             <- ifelse(mut_df$Var1 == "MG 1", "0", "1")

p <- ggplot(mut_df, aes(Var1, Var2)) +
  geom_tile(aes(fill = as.factor(value), width=0.97, height=0.97), size = 0.5, alpha=0.75) +
  #geom_tile(aes(fill = as.factor(value), color = as.factor(bordercol), width=0.98, height=0.98), size = 0.5, alpha=0.75) +
  ylab("Cell Lines") +
  xlab("Mutated Genes") +
  labs(fill = "Genotype") +
  theme_bw()+
  theme(plot.margin = unit(c(1, 0.5, 1, 0.75), "lines"),
        panel.grid.major= element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y=element_text(angle=0, hjust=1,size=11,colour="#525252"),
        axis.text.x=element_text(angle=50, hjust=1,size=11,colour="#525252"),
        axis.title.x=element_text(angle=0,size=13,face="bold",vjust=1,colour="#525252"),
        axis.title.y=element_text(angle=90, size=13,face="bold",vjust=0,colour="#525252"),
        legend.position = "right",
        legend.title = element_text(hjust=1,angle=0, size=12,colour="#525252"),
        legend.text = element_text(size=11,colour="#525252"),
        legend.key.size = unit(0.5, "cm"))+
  #scale_fill_gradient(low="white",high="#055C95") +
  scale_fill_manual(values = c("#DEDEDE","#055C95"), labels = c("WT", "Mut")) +
  scale_color_manual(values = c("#525252","white"), guide = FALSE)


p <- ggdraw(p) + draw_label("*", x = 0.32, y = 0.96, colour = "#098CCB", size = 15)
# p <- ggdraw(p) + draw_label("*", x = 0.29, y = 0.96, colour = "#098CCB", size = 15)
# ggsave(p, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/mutations",Sys.Date(),".pdf"),
#        width = 3.8, height = 6)


per_mat                    <- matrix(rnorm(n_celllines * n_perturbs, 0, 0.4), n_perturbs, n_celllines)
dimnames(per_mat)          <- list(perturbs,celllines)
mut1_celllines             <- celllines[which(mut_mat[1,] == 1)]
per_mat[8, mut1_celllines] <- rnorm(length(mut1_celllines), -2., 0.25)
per_df                     <- melt(per_mat[1:15,])
per_df$bordercol           <- ifelse(per_df$Var1 == "PG 8", "0", "1")

q <- ggplot(per_df, aes(Var1, Var2)) +
  geom_tile(aes(fill = value, width=0.97, height=0.97), size = 0.5, alpha=0.75) +
  #geom_tile(aes(fill = value, color = as.factor(bordercol), width=0.97, height=0.97), size = 0.5, alpha=0.75) +
  ylab("Cell Lines") +
  xlab("Perturbed genes") +
  theme_bw()+
  labs(fill = "Viability") +
  theme(plot.margin = unit(c(1, 0.75, 1, 0.5), "lines"),
        panel.grid.major= element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y=element_text(angle=0, hjust=1,size=11,colour="#525252"),
        axis.text.x=element_text(angle=50, hjust=1,size=11,colour="#525252"),
        axis.title.x=element_text(angle=0,size=13,face="bold",vjust=1,colour="#525252"),
        axis.title.y=element_text(angle=90, size=13,face="bold",vjust=0,colour="#525252"),
        legend.position = "right",
        legend.title = element_text(vjust=0.8,size=12,colour="#525252"),
        legend.text = element_text(size=11,colour="#525252"),
        legend.key.size = unit(0.6, "cm"))+
  scale_fill_viridis(option="viridis",begin = 0.2,end = 0.85,alpha = 0.75) +
  scale_color_manual(values = c("#525252","white"), guide = FALSE)

# pdf(file = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/perturbations",Sys.Date(),".pdf"), width = 10, height = 6)
# q
# grid.text("*", x = unit(0.49, "npc"), y = unit(0.965, "npc"), gp=gpar(fontsize=22, col="#8F55E1"))
# dev.off()

q <- ggdraw(q) + draw_label("*", x = 0.485, y = 0.96, colour = "#7868B2", size = 15)
#ggsave(q, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/perturbations",Sys.Date(),".pdf"),
#       width = 10, height = 6)


# Violin plots
violin_df         <- melt(per_mat)
#violin_df$class   <- sapply(violin_df$Var2, function(x){max(mut_df$value[mut_df$Var2 == x])})
violin_df$class   <- ifelse(violin_df$Var2 %in% mut1_celllines, "Mut", "WT")
subset_df         <- violin_df %>% dplyr::filter(Var1 == "PG 8")
subset_df$outlier <- ifelse(subset_df$Var2 %in% mut1_celllines, "Mut", "WT")
subset_df$names   <- paste0(subset_df$Var1, " (", subset_df$Var2,")")
mut1              <- mutations[1]

v <- ggplot(violin_df, aes(x=factor(class), y=value)) +
  geom_violin(aes(fill=factor(class),colour=factor(class)),alpha=0.5,position=position_dodge(0.9), size=0.9)+
  geom_boxplot(aes(colour=factor(class)),width=0.1,outlier.shape=NA, size=0.7) +
  ggtitle("MG 1 - PG 8 SL pair") +
  ylab("Viabilities of all perturbed genes") +
  xlab(paste0(mut1, " genotype"))+
  ylim(-3,2) +
  theme_bw()+
  theme(panel.grid= element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_line(size=0.5,color="#525252"),
        axis.line = element_line(colour="#525252"),
        axis.text.y=element_text(angle=0, vjust=0.3, hjust=0.5,size=11,colour="#525252"),
        axis.text.x=element_text(angle=0, vjust=0.1, hjust=0.5,size=11,colour="#525252"),
        axis.title.x=element_text(angle=0,size=13,face="bold",vjust=0,colour="#525252"),
        axis.title.y=element_text(angle=90, size=13,face="bold",vjust=1,colour="#525252"),
        plot.title = element_text(angle=0, size=14,face="bold",hjust = 0.5,colour="#525252"),
        legend.position="none")+
  scale_color_manual(name='method', values=c("#055C95", "#A8AFB0" ))+
  scale_fill_manual(name='method', values=c("white", "white")) +
  geom_point(data = subset_df, aes(x=factor(outlier), y=value),
             color = '#B30D31',
             shape = 16,
             size = 2.5) +
  geom_text(data = subset_df, aes(x=factor(outlier), y=value, label=names),
            hjust=-0.15,
            vjust=0.5,
            color = "#525252",
            alpha = 0.8,
            fontface = "bold",
            size = 3.25)

ggsave(v, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/violinPlot",Sys.Date(),".pdf"),
       width = 7, height = 7)


row_1 <- plot_grid(NULL, b, labels = c("A","B"), nrow = 1, label_size = 13)
row_2 <- plot_grid(NULL, p, q, labels = c("C","",""), rel_widths = c(0.05,0.45,1), nrow = 1, label_size = 13)
row_3 <- plot_grid(NULL,v,NULL, labels = c("","",""), rel_widths = c(0.5,0.6,0.5), nrow = 1, label_size = 13)
fin_plot <- plot_grid(row_1, row_2, NULL,row_3, nrow = 4, rel_heights = c(1,0.75,0.06,0.75))

ggsave(fin_plot, filename = paste0("/Volumes/beerenwinkel/ssumana/Documents/ETH/CRISPR/SLIDR/Figures/finalPlot",Sys.Date(),".pdf"),
       width = 11, height = 15)

