### Script for Barplot/boxplot for bulkseq

library(ggpubr)
library("tidyverse")
library("readxl")
library(openxlsx)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(stringr)
library(ggpattern)
library(ragg)
library(ggsignif)
library(ggsci)
library(ggforce)

counts <- read.csv("counts_normalizedCPM.csv")
counts.DEG <- read.csv("DEGs_fdr005_fc1.csv")

colnames(counts.DEG)[1] <- "gene"

counts.DEG$pvalue[counts.DEG$meta.q < 0.05] = "*"
counts.DEG$pvalue[counts.DEG$meta.q < 0.001] = "**"
counts.DEG$pvalue[counts.DEG$meta.q < 0.0001] = "***"
counts.DEG$pvalue[counts.DEG$meta.q < 0.00001] = "****"

counts.merge <- merge(counts.DEG, y = counts, by= "gene", all.x = TRUE, all.y = FALSE)
## delete useless column
counts.merge <- counts.merge[,-c(2:9,11:14)]

selected.genes <- read.xlsx("_Gene List To Graph.xlsx", sheet = "Sheet2")
selected.genes.split <- split.data.frame(selected.genes, as.factor(selected.genes$Pathway))


for (i in 1 : length(selected.genes.split)) {

## filter relevant genes from counts
geneList <- counts.merge %>% filter(gene %in% as.vector(selected.genes.split[[i]]$Genes))

## melt colums and create new group column
geneList <- melt(geneList,id.vars = c("gene", "pvalue"))

## create group column
geneList$group[grepl("Young", geneList$variable, fixed = TRUE)] = "young"
geneList$group[grepl("Old", geneList$variable, fixed = TRUE)] = "old"

geneList$group <- factor(geneList$group, levels = c("young", "old"))

agg_png(paste("boxplot_",unique(selected.genes.split[[i]]$Pathway),".png", sep = ""), width = 12, height = 2.6, res = 300, units = "in", background = "transparent")
#cairo_pdf(paste("boxplot_",unique(selected.genes.split[[i]]$Pathway),".pdf", sep = ""), width = 12, height = 2.6, pointsize = 12, fallback_resolution = 300)

print(ggplot(geneList, aes(x = group, y = value, fill = group, pattern = group, color = group)) +
  #geom_violin(aes(color = group), alpha= 0, trim=FALSE) +
  geom_boxplot_pattern(position = position_dodge(0.9),
                       #alpha = .5,
                       lwd = .25,
                       pattern_fill = "white",
                       pattern_color = 'white',
                       pattern_angle = 45,
                       pattern_size = 0.2,
                       pattern_density = 0.1,
                       pattern_spacing = 0.05,
                       pattern_key_scale_factor = 0.5) +
  geom_jitter(width=0.15, show.legend = FALSE) +
  geom_signif(data = geneList, aes(annotations = pvalue), comparisons=list(c("young", "old")), manual = TRUE,
              vjust = 0.4, tip_length = 0, size = 0.25, linetype = "solid") +
  #facet_wrap(~ gene, scales = "free", ncol = 10) +
  ggforce::facet_row(~ gene, scales = 'free', space = "fixed") +
  #ggtitle("IGF-1 signature top10 upregulated genes") +
  ylab("logCPM mRNA expression") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  #scale_fill_manual(values = c("white","white","#5c99d6","#5c99d6")) +
  #scale_fill_npg(alpha = .5) +
  scale_fill_manual(values = c("#c8def9", "#ff7171")) +
  scale_color_manual(values = c("#000000","#000000")) +
  scale_pattern_manual(values = c(young = "none", old = "none")) +
  #scale_linetype_manual(values=c("dashed", "solid", "dashed", "solid","solid","solid")) +
  #guides(pattern = guide_legend(override.aes = list(pattern_alpha = 1))) +
  theme_classic(
    base_family = "Arial"
  ) + theme(
    axis.line = element_line(colour = 'black', size = .25),
    axis.title.x=element_blank(),
    legend.title = element_blank(),
    legend.text = element_text( size = 14),
    #axis.text.x=element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.text.x=element_blank(),
    #axis.text.y=element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold.italic', size = 11),
    axis.title.y = element_text(size = 14),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent'),
    aspect.ratio = 3/1
  ))

invisible(dev.off())
}
