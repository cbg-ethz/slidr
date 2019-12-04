library(ggplot2)
library(sysfonts)
font_add_google("Muli")
library(ggimage)
library(hexSticker)


load("~/Downloads/Slidr_Results/sticker_data.Rdata")
# Generate boxplot
p <- ggplot(temp_data,aes(factor(mut_status),viabilities)) +
  geom_boxplot(fill = "white", color = "grey30",
               alpha=0.4,
               size = 0.1,
               width = 0.3,
               outlier.shape=NA) +
  geom_jitter(aes(colour = CN_type, fill = mut_type),
              position = position_jitter(width = 0.2, height = 0.05),
              shape = 21,
              size = 0.5,
              stroke = 0.3) +
  ylab(paste("Gene B")) +
  xlab("Gene A")+
  scale_color_manual(name = "CNA", values=subset_CNA_colors, labels = names(subset_CNA_colors))+
  scale_fill_manual(name="Mutation type",values=subset_colors, labels = gsub("_"," ", names(subset_colors))) +
  theme_transparent() +
  ylim(-2,1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.2, color="#DCDCDC"),
        axis.ticks = element_line(size = 0.1,color="#DCDCDC"),
        axis.text.y=element_text(angle=0, vjust=0.3, hjust=0.5,size=4,colour="#DCDCDC"),
        axis.text.x=element_text(angle=0, vjust=0.1, hjust=0,size=4,colour="#DCDCDC"),
        axis.title.x=element_text(angle=0, size=5,vjust=0,colour="#DCDCDC"),
        axis.title.y=element_text(angle=90, size=5,vjust=0,colour="#DCDCDC"),
        legend.position="none")

# generate sticker
sticker(p,
        s_x = 0.9,
        s_y = 0.8,
        s_w = 1.1,
        s_h = 1.1,
        h_size = 0.5,
        h_fill = "#121214",
        h_color = "#42454C",
        package= "slidr",
        p_color = "#DCDCDC",
        p_size = 7,
        p_x = 1,
        p_y = 1.5,
        p_family = "Muli",
        spotlight = TRUE,
        l_x = 1,
        l_y = 1.8,
        l_alpha = 0.13,
        l_width = 10,
        l_height = 5)
