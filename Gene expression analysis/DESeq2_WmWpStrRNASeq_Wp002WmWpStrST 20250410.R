#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")

### DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

### EnhacedVolcano (https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) and (https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)
#devtools::install_github('kevinblighe/EnhancedVolcano')
#OR
#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#BiocManager::install('EnhancedVolcano', force = TRUE)

### Pheatmap (https://cran.r-project.org/web/packages/pheatmap/index.html)
#install.packages("pheatmap")

### tidyverse (https://www.tidyverse.org/)
#install.packages("tidyverse")

### RColorBrewer (https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf)
install.packages("RColorBrewer")

install.packages("readr")
install.packages("factoextra")

library(tximport)
library(readr)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap) 
library(factoextra)
library(viridis)

setwd("Y:/hxue/W.paniculata_W.multiflora_RNA-Seq")



#txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

files_WpWm <- c(
  "quants_Wp002WmWpStrST/EB_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/EB_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/LB_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/LB_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/FL_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/FL_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRc_quant/quant.sf",   
  "quants_Wp002WmWpStrST/Wp_MiSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARc_quant/quant.sf",
"quants_Wp002WmWpStrST/L_ES_A_04_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_B_05_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_C_06_quant/quant.sf",
"quants_Wp002WmWpStrST/L_MS_A_16_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_B_17_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_C_18_quant/quant.sf",
"quants_Wp002WmWpStrST/L_LS_A_28_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_B_29_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_C_30_quant/quant.sf",
"quants_Wp002WmWpStrST/L_EA_A_10_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_B_11_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_C_12_quant/quant.sf",
"quants_Wp002WmWpStrST/L_MA_A_22_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_B_23_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_C_24_quant/quant.sf",
"quants_Wp002WmWpStrST/L_LA_A_34_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_B_35_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_C_36_quant/quant.sf",
"quants_Wp002WmWpStrST/R_ES_A_01_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_B_02_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_C_03_quant/quant.sf",
"quants_Wp002WmWpStrST/R_MS_A_13_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_B_14_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_C_15_quant/quant.sf",
"quants_Wp002WmWpStrST/R_LS_A_25_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_B_26_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_C_27_quant/quant.sf",
"quants_Wp002WmWpStrST/R_EA_A_07_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_B_08_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_C_09_quant/quant.sf",
"quants_Wp002WmWpStrST/R_MA_A_19_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_B_20_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_C_21_quant/quant.sf",
"quants_Wp002WmWpStrST/R_LA_A_31_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_B_32_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_C_33_quant/quant.sf")

Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

txi_WpWm <- tximport(files_WpWm, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)

sampleID_WpWm = c(
  "OWpEaSLa","OWpEaSLb","OWpEaSLc",
  "OWpEaSRa","OWpEaSRb","OWpEaSRc",
  "OWpLaSLa","OWpLaSLb","OWpLaSLc",
  "OWpLaSRa","OWpLaSRb","OWpLaSRc",
  "OWpFLSLa","OWpFLSLb","OWpFLSLc",
  "OWpFLSRa","OWpFLSRb","OWpFLSRc",
  "WpEaSLa","WpEaSLb","WpEaSLc",
  "WpMiSLa","WpMiSLb","WpMiSLc",
  "WpLaSLa","WpLaSLb","WpLaSLc",
  "WpEaALa","WpEaALb","WpEaALc",
  "WpMiALa","WpMiALb","WpMiALc",
  "WpLaALa","WpLaALb","WpLaALc",
  "WpEaSRa","WpEaSRb","WpEaSRc",
  "WpMiSRa","WpMiSRb","WpMiSRc",
  "WpLaSRa","WpLaSRb","WpLaSRc",
  "WpEaARa","WpEaARb","WpEaARc",
  "WpMiARa","WpMiARb","WpMiARc",
  "WpLaARa","WpLaARb","WpLaARc",
  "WmEaSLa","WmEaSLb","WmEaSLc",
  "WmMiSLa","WmMiSLb","WmMiSLc",
  "WmLaSLa","WmLaSLb","WmLaSLc",
  "WmEaALa","WmEaALb","WmEaALc",
  "WmMiALa","WmMiALb","WmMiALc",
  "WmLaALa","WmLaALb","WmLaALc",
  "WmEaSRa","WmEaSRb","WmEaSRc",
  "WmMiSRa","WmMiSRb","WmMiSRc",
  "WmLaSRa","WmLaSRb","WmLaSRc",
  "WmEaARa","WmEaARb","WmEaARc",
  "WmMiARa","WmMiARb","WmMiARc",
  "WmLaARa","WmLaARb","WmLaARc"
)
type_WpWm = c(
  "WpLSEaO","WpLSEaO","WpLSEaO",
  "WpRSEaO","WpRSEaO","WpRSEaO",
  "WpLSLaO","WpLSLaO","WpLSLaO",
  "WpRSLaO","WpRSLaO","WpRSLaO",
  "WpLSOFO","WpLSOFO","WpLSOFO",
  "WpRSOFO","WpRSOFO","WpRSOFO",
  "WpLSEa","WpLSEa","WpLSEa",
  "WpLSMi","WpLSMi","WpLSMi",
  "WpLSLs","WpLSLa","WpLSLa",
  "WpLAEa","WpLAEa","WpLAEa",
  "WpLAMi","WpLAMi","WpLAMi",
  "WpLALa","WpLALa","WpLALa",
  "WpRSEa","WpRSEa","WpRSEa",
  "WpRSMi","WpRSMi","WpRSMi",
  "WpRSLa","WpRSLa","WpRSLa",
  "WpRAEa","WpRAEa","WpRAEa",
  "WpRAMi","WpRAMi","WpRAMi",
  "WpRALa","WpRALa","WpRALa",
  "WmLSEa","WmLSEa","WmLSEa",
  "WmLSMi","WmLSMi","WmLSMi",
  "WmLSLa","WmLSLa","WmLSLa",
  "WmLAEa","WmLAEa","WmLAEa",
  "WmLAMi","WmLAMi","WmLAMi",
  "WmLALa","WmLALa","WmLALa",
  "WmRSEa","WmRSEa","WmRSEa",
  "WmRSMi","WmRSMi","WmRSMi",
  "WmRSLa","WmRSLa","WmRSLa",
  "WmRAEa","WmRAEa","WmRAEa",
  "WmRAMi","WmRAMi","WmRAMi",
  "WmRALa","WmRALa","WmRALa"
  )



colData_WpWm <- data.frame(sampleID_WpWm, type_WpWm)
dds_WpWm <- DESeqDataSetFromTximport(txi_WpWm, colData_WpWm, design= ~ type_WpWm)

dds_WpWm <- estimateSizeFactors(dds_WpWm)
normalized_counts_Wp <- counts(dds_Wp, normalized=TRUE)
colnames(normalized_counts_Wp) <- sampleID_Wp
write.csv(normalized_counts_Wp, file = "Wp_normalized_gene_count_matrix.csv")

##PCA
vsd_WpWm <- vst(dds_WpWm, blind=FALSE)
z <- DESeq2::plotPCA(vsd_WpWm, intgroup = "type_WpWm", ntop = 500, returnData = FALSE)
z$data$species_organ = c(
  "1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS",
  "1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS",
  "1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS",
  "2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA",
  "1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS",
  "2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA",
  "3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS",
  "4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA",
  "3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS",
  "4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA"
)
z$data$stage_morph_batch = as.factor(c(
  "1ELO","1ELO","1ELO",
  "1ERO","1ERO","1ERO",
  "3LLO","3LLO","3LLO",
  "3LRO","3LRO","3LRO",
  "4FLO","4FLO","4FLO",
  "4FRO","4FRO","4FRO",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR"
))
#nudge <- position_nudge(y = 1)
#z + geom_text(aes(label = sampleID_Wp), position = nudge) + ggtitle("W. paniculata + W. multiflora RNA-seq")
#z + ggtitle("W. paniculata W. multiflora RNA-seq")

color_palette <- c(turbo(12))

ggplot(z$data, aes(x=PC1, y=PC2, color=stage_morph_batch, shape=species_organ, )) + 
  geom_point(size=2, stroke = 1.5) +
  scale_shape_manual(values = c("1WpS" = 16, "2WpA"= 17, "3WmS" = 1, "4WmA" = 2), name = "Species and organ") +
  scale_color_manual(values = color_palette, name = "Stage, morph, and batch") +
  theme_bw()+
  xlab("PC1: 28% variance") +
  ylab("PC2: 19% variance") +
  ggtitle("W. paniculata + W. multiflora RNA-seq") +
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=13)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=13)) 



##count matrix
dds_Wp <- estimateSizeFactors(dds_Wp)
normalized_counts_Wp <- counts(dds_Wp, normalized=TRUE)
colnames(normalized_counts_Wp) <- sampleID_Wp
write.csv(normalized_counts_Wp, file = "Wp_normalized_gene_count_matrix.csv")

res.pca_Wp <- prcomp(normalized_counts_Wp)

fviz_pca_ind(res.pca_Wp,
#  col.ind = "cos2", # Color by the quality of representation
#  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE     # Avoid text overlapping
)

plot(res.pca_Wp$rotation[,1],res.pca_Wp$rotation[,2])

png(file="Wp.Wp002WmWpStrST.817(miR156-R1).png", width=400, height=400)
plotCounts(dds_WpWm, gene="Wp002WmWpStrST.817", transform = FALSE, intgroup = "type_WpWm", returnData = FALSE)
dev.off()

png(file="Wp.Wp002WmWpStrST.812(YUC-R1).png", width=400, height=400)
plotCounts(dds_WpWm, gene="Wp002WmWpStrST.812", transform = FALSE, intgroup = "type_WpWm", returnData = FALSE)
dev.off()

png(file="Wp.Wp002WmWpStrST.813(anti-YUC-R1).png", width=400, height=400)
plotCounts(dds_WpWm, gene="Wp002WmWpStrST.813", transform = FALSE, intgroup = "type_WpWm", returnData = FALSE)
dev.off()

png(file="Wp.Wp002WmWpStrST.53482(YUC-P1).png", width=2400, height=800)
plotCounts(dds_WpWm, gene="Wp002WmWpStrST.53482", transform = FALSE, intgroup = "type_WpWm", returnData = FALSE)
dev.off()

png(file="Wp.Wp002WmWpStrST.12320(YUC-P2).png", width=2400, height=800)
plotCounts(dds_WpWm, gene="Wp002WmWpStrST.12320", transform = FALSE, intgroup = "type_WpWm", returnData = FALSE)
dev.off()

png(file="WpWm.Wp002WmWpStrST.30480(QRT2-1).png", width=2400, height=800)
plotCounts(dds_WpWm, gene="Wp002WmWpStrST.30480", transform = FALSE, intgroup = "type_WpWm", returnData = FALSE)
dev.off()


files_Wp <- c(
  "quants_Wp002WmWpStrST/EB_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/EB_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/LB_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/LB_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/FL_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/FL_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRc_quant/quant.sf",   
  "quants_Wp002WmWpStrST/Wp_MiSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARc_quant/quant.sf")

Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

txi_Wp <- tximport(files_Wp, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)

sampleID_Wp = c(
  "OWpEaSLa","OWpEaSLb","OWpEaSLc",
  "OWpEaSRa","OWpEaSRb","OWpEaSRc",
  "OWpLaSLa","OWpLaSLb","OWpLaSLc",
  "OWpLaSRa","OWpLaSRb","OWpLaSRc",
  "OWpFLSLa","OWpFLSLb","OWpFLSLc",
  "OWpFLSRa","OWpFLSRb","OWpFLSRc",
  "WpEaSLa","WpEaSLb","WpEaSLc",
  "WpMiSLa","WpMiSLb","WpMiSLc",
  "WpLaSLa","WpLaSLb","WpLaSLc",
  "WpEaALa","WpEaALb","WpEaALc",
  "WpMiALa","WpMiALb","WpMiALc",
  "WpLaALa","WpLaALb","WpLaALc",
  "WpEaSRa","WpEaSRb","WpEaSRc",
  "WpMiSRa","WpMiSRb","WpMiSRc",
  "WpLaSRa","WpLaSRb","WpLaSRc",
  "WpEaARa","WpEaARb","WpEaARc",
  "WpMiARa","WpMiARb","WpMiARc",
  "WpLaARa","WpLaARb","WpLaARc"
)
type_Wp = c(
  "WpEaSLO","WpEaSLO","WpEaSLO",
  "WpEaSRO","WpEaSRO","WpEaSRO",
  "WpLaSLO","WpLaSLO","WpLaSLO",
  "WpLaSRO","WpLaSRO","WpLaSRO",
  "WpOFSLO","WpOFSLO","WpOFSLO",
  "WpOFSRO","WpOFSRO","WpOFSRO",
  "WpEaSL","WpEaSL","WpEaSL",
  "WpMiSL","WpMiSL","WpMiSL",
  "WpLaSL","WpLaSL","WpLaSL",
  "WpEaAL","WpEaAL","WpEaAL",
  "WpMiAL","WpMiAL","WpMiAL",
  "WpLaAL","WpLaAL","WpLaAL",
  "WpEaSR","WpEaSR","WpEaSR",
  "WpMiSR","WpMiSR","WpMiSR",
  "WpLaSR","WpLaSR","WpLaSR",
  "WpEaAR","WpEaAR","WpEaAR",
  "WpMiAR","WpMiAR","WpMiAR",
  "WpLaAR","WpLaAR","WpLaAR"
)

colData_Wp <- data.frame(sampleID_Wp, type_Wp)
dds_Wp <- DESeqDataSetFromTximport(txi_Wp, colData_Wp, design= ~ type_Wp)

##PCA
vsd_Wp <- vst(dds_Wp, blind=FALSE)
z <- DESeq2::plotPCA(vsd_Wp, intgroup = "type_Wp", ntop = 500, returnData = FALSE)
#nudge <- position_nudge(y = 1)
#z + geom_text(aes(label = sampleID_Wp), position = nudge) + ggtitle("W. paniculata + W. multiflora RNA-seq")
#z + ggtitle("W. paniculata RNA-seq") + geom_point(data=z$data, )

z$data$organ = c(
  "Style","Style","Style","Style","Style","Style","Style","Style","Style",
  "Style","Style","Style","Style","Style","Style","Style","Style","Style",
  "Style","Style","Style","Style","Style","Style","Style","Style","Style",
  "Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen",
  "Style","Style","Style","Style","Style","Style","Style","Style","Style",
  "Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen"
)
z$data$stage_morph_batch = as.factor(c(
  "1ELO","1ELO","1ELO",
  "1ERO","1ERO","1ERO",
  "3LLO","3LLO","3LLO",
  "3LRO","3LRO","3LRO",
  "4FLO","4FLO","4FLO",
  "4FRO","4FRO","4FRO",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR"
))
#nudge <- position_nudge(y = 1)
#z + geom_text(aes(label = sampleID_Wp), position = nudge) + ggtitle("W. paniculata + W. multiflora RNA-seq")
#z + ggtitle("W. paniculata W. multiflora RNA-seq")

color_palette <- c(turbo(12))

ggplot(z$data, aes(x=PC1, y=PC2, color=stage_morph_batch, shape=organ, )) + 
  geom_point(size=2, stroke = 1.5) +
  scale_shape_manual(values = c("Style" = 16, "Stamen"= 17), name = "Organ") +
  scale_color_manual(values = color_palette, name = "Stage, morph, and batch") +
  theme_bw()+
  xlab("PC1: 38% variance") +
  ylab("PC2: 20% variance") +
  ggtitle("W. paniculata RNA-seq") +
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=13)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=13)) 

files_WpStrWm <- c(
  "quants_Wp002WmWpStrST/Wp_EaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRc_quant/quant.sf",   
  "quants_Wp002WmWpStrST/Wp_MiSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_ES_A_04_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_B_05_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_C_06_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_MS_A_16_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_B_17_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_C_18_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_LS_A_28_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_B_29_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_C_30_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_EA_A_10_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_B_11_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_C_12_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_MA_A_22_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_B_23_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_C_24_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_LA_A_34_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_B_35_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_C_36_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_ES_A_01_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_B_02_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_C_03_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_MS_A_13_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_B_14_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_C_15_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_LS_A_25_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_B_26_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_C_27_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_EA_A_07_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_B_08_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_C_09_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_MA_A_19_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_B_20_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_C_21_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_LA_A_31_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_B_32_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_C_33_quant/quant.sf")

Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

txi_WpStrWm <- tximport(files_WpStrWm, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)

sampleID_WpStrWm = c(
  "WpEaSLa","WpEaSLb","WpEaSLc",
  "WpMiSLa","WpMiSLb","WpMiSLc",
  "WpLaSLa","WpLaSLb","WpLaSLc",
  "WpEaALa","WpEaALb","WpEaALc",
  "WpMiALa","WpMiALb","WpMiALc",
  "WpLaALa","WpLaALb","WpLaALc",
  "WpEaSRa","WpEaSRb","WpEaSRc",
  "WpMiSRa","WpMiSRb","WpMiSRc",
  "WpLaSRa","WpLaSRb","WpLaSRc",
  "WpEaARa","WpEaARb","WpEaARc",
  "WpMiARa","WpMiARb","WpMiARc",
  "WpLaARa","WpLaARb","WpLaARc",
  "WmEaSLa","WmEaSLb","WmEaSLc",
  "WmMiSLa","WmMiSLb","WmMiSLc",
  "WmLaSLa","WmLaSLb","WmLaSLc",
  "WmEaALa","WmEaALb","WmEaALc",
  "WmMiALa","WmMiALb","WmMiALc",
  "WmLaALa","WmLaALb","WmLaALc",
  "WmEaSRa","WmEaSRb","WmEaSRc",
  "WmMiSRa","WmMiSRb","WmMiSRc",
  "WmLaSRa","WmLaSRb","WmLaSRc",
  "WmEaARa","WmEaARb","WmEaARc",
  "WmMiARa","WmMiARb","WmMiARc",
  "WmLaARa","WmLaARb","WmLaARc"
)
type_WpStrWm = c(
  "WpEaSL","WpEaSL","WpEaSL",
  "WpMiSL","WpMiSL","WpMiSL",
  "WpLaSL","WpLaSL","WpLaSL",
  "WpEaAL","WpEaAL","WpEaAL",
  "WpMiAL","WpMiAL","WpMiAL",
  "WpLaAL","WpLaAL","WpLaAL",
  "WpEaSR","WpEaSR","WpEaSR",
  "WpMiSR","WpMiSR","WpMiSR",
  "WpLaSR","WpLaSR","WpLaSR",
  "WpEaAR","WpEaAR","WpEaAR",
  "WpMiAR","WpMiAR","WpMiAR",
  "WpLaAR","WpLaAR","WpLaAR",
  "WmEaSL","WmEaSL","WmEaSL",
  "WmMiSL","WmMiSL","WmMiSL",
  "WmLaSL","WmLaSL","WmLaSL",
  "WmEaAL","WmEaAL","WmEaAL",
  "WmMiAL","WmMiAL","WmMiAL",
  "WmLaAL","WmLaAL","WmLaAL",
  "WmEaSR","WmEaSR","WmEaSR",
  "WmMiSR","WmMiSR","WmMiSR",
  "WmLaSR","WmLaSR","WmLaSR",
  "WmEaAR","WmEaAR","WmEaAR",
  "WmMiAR","WmMiAR","WmMiAR",
  "WmLaAR","WmLaAR","WmLaAR"
)



colData_WpStrWm <- data.frame(sampleID_WpStrWm, type_WpStrWm)
dds_WpStrWm <- DESeqDataSetFromTximport(txi_WpStrWm, colData_WpStrWm, design= ~ type_WpStrWm)

##PCA
vsd_WpStrWm <- vst(dds_WpStrWm, blind=FALSE)
z <- DESeq2::plotPCA(vsd_WpStrWm, intgroup = "type_WpStrWm", ntop = 500, returnData = FALSE)
z$data$species_organ = c(
  "1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS",
  "2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA",
  "1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS","1WpS",
  "2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA","2WpA",
  "3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS",
  "4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA",
  "3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS","3WmS",
  "4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA","4WmA"
)
z$data$stage_morph = as.factor(c(
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR"
))
nudge <- position_nudge(y = 1)
z + geom_text(aes(label = sampleID_WpStrWm), position = nudge) + ggtitle("W. paniculata + W. multiflora stranded RNA-seq")
z + ggtitle("W. paniculata + W. multiflora RNA-seq")

color_palette <- c(turbo(6))

ggplot(z$data, aes(x=PC1, y=PC2, color=stage_morph, shape=species_organ, )) + 
  geom_point(size=2, stroke = 1.5) +
  scale_shape_manual(values = c("1WpS" = 16, "2WpA"= 17, "3WmS" = 1, "4WmA" = 2), name = "Species and organ") +
  scale_color_manual(values = color_palette, name = "Stage, morph, and batch") +
  theme_bw()+
  xlab("PC1: 36% variance") +
  ylab("PC2: 17% variance") +
  ggtitle("Stranded W. paniculata + W. multiflora RNA-seq") +
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=13)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=13)) 




files_WpStr <- c(
  "quants_Wp002WmWpStrST/Wp_EaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRc_quant/quant.sf",   
  "quants_Wp002WmWpStrST/Wp_MiSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARc_quant/quant.sf")

Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

txi_WpStr <- tximport(files_WpStr, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)

sampleID_WpStr = c(
  "WpEaSLa","WpEaSLb","WpEaSLc",
  "WpMiSLa","WpMiSLb","WpMiSLc",
  "WpLaSLa","WpLaSLb","WpLaSLc",
  "WpEaALa","WpEaALb","WpEaALc",
  "WpMiALa","WpMiALb","WpMiALc",
  "WpLaALa","WpLaALb","WpLaALc",
  "WpEaSRa","WpEaSRb","WpEaSRc",
  "WpMiSRa","WpMiSRb","WpMiSRc",
  "WpLaSRa","WpLaSRb","WpLaSRc",
  "WpEaARa","WpEaARb","WpEaARc",
  "WpMiARa","WpMiARb","WpMiARc",
  "WpLaARa","WpLaARb","WpLaARc"
)
type_WpStr = c(
  "WpLe1S1Ea","WpLe1S1Ea","WpLe1S1Ea",
  "WpLe1S2Mi","WpLe1S2Mi","WpLe1S2Mi",
  "WpLe1S3La","WpLe1S3La","WpLe1S3La",
  "WpLe2A1Ea","WpLe2A1Ea","WpLe2A1Ea",
  "WpLe2A2Mi","WpLe2A2Mi","WpLe2A2Mi",
  "WpLe2A3LA","WpLe2A3LA","WpLe2A3LA",
  "WpRi1S1Ea","WpRi1S1Ea","WpRi1S1Ea",
  "WpRi1S2Mi","WpRi1S2Mi","WpRi1S2Mi",
  "WpRi1S3La","WpRi1S3La","WpRi1S3La",
  "WpRi2A1Ea","WpRi2A1Ea","WpRi2A1Ea",
  "WpRi2A2Mi","WpRi2A2Mi","WpRi2A2Mi",
  "WpRi2A3LA","WpRi2A3LA","WpRi2A3LA"
)

colData_WpStr <- data.frame(sampleID_WpStr, type_WpStr)
dds_WpStr <- DESeqDataSetFromTximport(txi_WpStr, colData_WpStr, design= ~ type_WpStr)

##PCA
vsd_WpStr <- vst(dds_WpStr, blind=FALSE)
z <- DESeq2::plotPCA(vsd_WpStr, intgroup = "type_WpStr", ntop = 500, returnData = FALSE)
nudge <- position_nudge(y = 1)
z + geom_text(aes(label = sampleID_WpStr), position = nudge) + ggtitle("Stranded W. paniculata RNA-seq")


z$data$organ = c(
  "Style","Style","Style","Style","Style","Style","Style","Style","Style",
  "Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen",
  "Style","Style","Style","Style","Style","Style","Style","Style","Style",
  "Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen"
)
z$data$stage_morph = as.factor(c(
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR"
))
#nudge <- position_nudge(y = 1)
#z + geom_text(aes(label = sampleID_Wp), position = nudge) + ggtitle("W. paniculata + W. multiflora RNA-seq")
#z + ggtitle("W. paniculata W. multiflora RNA-seq")

color_palette <- c(turbo(6))

ggplot(z$data, aes(x=PC1, y=PC2, color=stage_morph, shape=organ, )) + 
  geom_point(size=2, stroke = 1.5) +
  scale_shape_manual(values = c("Style" = 16, "Stamen"= 17), name = "Organ") +
  scale_color_manual(values = color_palette, name = "Stage and morph") +
  theme_bw()+
  xlab("PC1: 35% variance") +
  ylab("PC2: 23% variance") +
  ggtitle("W. paniculata stranded RNA-seq") +
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=13)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=13)) 




files_Wm <- c(
  "quants_Wp002WmWpStrST/L_ES_A_04_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_B_05_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_C_06_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_MS_A_16_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_B_17_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_C_18_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_LS_A_28_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_B_29_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_C_30_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_EA_A_10_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_B_11_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_C_12_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_MA_A_22_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_B_23_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_C_24_quant/quant.sf",
  "quants_Wp002WmWpStrST/L_LA_A_34_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_B_35_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_C_36_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_ES_A_01_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_B_02_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_C_03_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_MS_A_13_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_B_14_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_C_15_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_LS_A_25_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_B_26_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_C_27_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_EA_A_07_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_B_08_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_C_09_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_MA_A_19_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_B_20_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_C_21_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_LA_A_31_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_B_32_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_C_33_quant/quant.sf")

Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

txi_Wm <- tximport(files_Wm, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)

sampleID_Wm = c(
  "WmEaSLa","WmEaSLb","WmEaSLc",
  "WmMiSLa","WmMiSLb","WmMiSLc",
  "WmLaSLa","WmLaSLb","WmLaSLc",
  "WmEaALa","WmEaALb","WmEaALc",
  "WmMiALa","WmMiALb","WmMiALc",
  "WmLaALa","WmLaALb","WmLaALc",
  "WmEaSRa","WmEaSRb","WmEaSRc",
  "WmMiSRa","WmMiSRb","WmMiSRc",
  "WmLaSRa","WmLaSRb","WmLaSRc",
  "WmEaARa","WmEaARb","WmEaARc",
  "WmMiARa","WmMiARb","WmMiARc",
  "WmLaARa","WmLaARb","WmLaARc"
)


type_Wm = c(
  "WmLe1S1Ea","WmLe1S1Ea","WmLe1S1Ea",
  "WmLe1S2Mi","WmLe1S2Mi","WmLe1S2Mi",
  "WmLe1S3La","WmLe1S3La","WmLe1S3La",
  "WmLe2A1Ea","WmLe2A1Ea","WmLe2A1Ea",
  "WmLe2A2Mi","WmLe2A2Mi","WmLe2A2Mi",
  "WmLe2A3LA","WmLe2A3LA","WmLe2A3LA",
  "WmRi1S1Ea","WmRi1S1Ea","WmRi1S1Ea",
  "WmRi1S2Mi","WmRi1S2Mi","WmRi1S2Mi",
  "WmRi1S3La","WmRi1S3La","WmRi1S3La",
  "WmRi2A1Ea","WmRi2A1Ea","WmRi2A1Ea",
  "WmRi2A2Mi","WmRi2A2Mi","WmRi2A2Mi",
  "WmRi2A3LA","WmRi2A3LA","WmRi2A3LA"
)



colData_Wm <- data.frame(sampleID_Wm, type_Wm)
dds_Wm <- DESeqDataSetFromTximport(txi_Wm, colData_Wm, design= ~ type_Wm)

##PCA
vsd_Wm <- vst(dds_Wm, blind=FALSE)
z <- DESeq2::plotPCA(vsd_Wm, intgroup = "type_Wm", ntop = 500, returnData = FALSE)
z$data$organ = c(
    "Style","Style","Style","Style","Style","Style","Style","Style","Style",
    "Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen",
    "Style","Style","Style","Style","Style","Style","Style","Style","Style",
    "Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen","Stamen"
)
z$data$stage_morph = as.factor(c(
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1EL","1EL","1EL",
  "2ML","2ML","2ML",
  "3LL","3LL","3LL",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR",
  "1ER","1ER","1ER",
  "2MR","2MR","2MR",
  "3LR","3LR","3LR"
))
nudge <- position_nudge(y = 1)
z + geom_text(aes(label = sampleID_Wm), position = nudge) + ggtitle("W. paniculata + W. multiflora stranded RNA-seq")
z + ggtitle("W. multiflora RNA-seq")

color_palette <- c(turbo(6))

ggplot(z$data, aes(x=PC1, y=PC2, color=stage_morph, shape=organ, )) + 
  geom_point(size=2, stroke = 1.5) +
  scale_shape_manual(values = c("Style" = 1, "Stamen"= 2), name = "Organ") +
  scale_color_manual(values = color_palette, name = "Stage and morph") +
  theme_bw()+
  xlab("PC1: 33% variance") +
  ylab("PC2: 27% variance") +
  ggtitle("W. multiflora RNA-seq") +
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=13)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=13)) 




files_WpO <- c(
  "quants_Wp002WmWpStrST/EB_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/EB_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/EB_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/LB_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/LB_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/LB_R_Z_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/FL_L_A_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_L_B_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_L_C_trimmed_quant/quant.sf",
  "quants_Wp002WmWpStrST/FL_R_X_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_R_Y_trimmed_quant/quant.sf","quants_Wp002WmWpStrST/FL_R_Z_trimmed_quant/quant.sf"
)
Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

txi_WpO <- tximport(files_WpO, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)

sampleID_WpO = c(
  "OWpEaSLa","OWpEaSLb","OWpEaSLc",
  "OWpEaSRa","OWpEaSRb","OWpEaSRc",
  "OWpLaSLa","OWpLaSLb","OWpLaSLc",
  "OWpLaSRa","OWpLaSRb","OWpLaSRc",
  "OWpFLSLa","OWpFLSLb","OWpFLSLc",
  "OWpFLSRa","OWpFLSRb","OWpFLSRc"
)
type_WpO = c(
  "WpEaSLO","WpEaSLO","WpEaSLO",
  "WpEaSRO","WpEaSRO","WpEaSRO",
  "WpLaSLO","WpLaSLO","WpLaSLO",
  "WpLaSRO","WpLaSRO","WpLaSRO",
  "WpOFSLO","WpOFSLO","WpOFSLO",
  "WpOFSRO","WpOFSRO","WpOFSRO"
)

colData_WpO <- data.frame(sampleID_WpO, type_WpO)
dds_WpO <- DESeqDataSetFromTximport(txi_WpO, colData_WpO, design= ~ type_WpO)

##PCA
vsd_WpO <- vst(dds_WpO, blind=FALSE)
z <- DESeq2::plotPCA(vsd_WpO, intgroup = "type_WpO", ntop = 500, returnData = FALSE)
#nudge <- position_nudge(y = 1)
#z + geom_text(aes(label = sampleID_Wp), position = nudge) + ggtitle("W. paniculata + W. multiflora RNA-seq")
z + ggtitle("W. paniculata RNA-seq") + geom_point(data=z$data, )

z$data$stage_morph = as.factor(c(
  "1EL","1EL","1EL",
  "1ER","1ER","1ER",
  "3LL","3LL","3LL",
  "3LR","3LR","3LR",
  "4FL","4FL","4FL",
  "4FR","4FR","4FR"
))
#nudge <- position_nudge(y = 1)
#z + geom_text(aes(label = sampleID_Wp), position = nudge) + ggtitle("W. paniculata + W. multiflora RNA-seq")
z + ggtitle("W. paniculata unstranded RNA-seq")

color_palette <- c(turbo(6))

ggplot(z$data, aes(x=PC1, y=PC2, color=stage_morph, )) + 
  geom_point(size=2, stroke = 1.5) +
  scale_color_manual(values = color_palette, name = "Stage and morph") +
  theme_bw()+
  xlab("PC1: 41% variance") +
  ylab("PC2: 24% variance") +
  ggtitle("W. paniculata unstranded RNA-seq") +
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=13)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=13)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=13)) 


##DGE analysis for WpStr


Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

files_WpES <- c(
  "quants_Wp002WmWpStrST/Wp_EaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRc_quant/quant.sf"
)
txi_WpES <- tximport(files_WpES, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WmES = c(
  "WpEaSLa","WpEaSLb","WpEaSLc",
  "WpEaSRa","WpEaSRb","WpEaSRc"
)
colData_WpES <- data.frame(sampleID_WpES, type_LR)
dds_WpES <- DESeqDataSetFromTximport(txi_WpES, colData_WmES, design= ~ type_LR)
dds_WpES <- DESeq(dds_WpES)
WmES_DESeq2 <- results(dds_WmES)
write.table(WmES_DESeq2, file = "WmES_DESeq2")


files_WpMA <- c(
  "quants_Wp002WmWpStrST/Wp_MiALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARc_quant/quant.sf"
)
txi_WpMA <- tximport(files_WpMA, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WpMA = c(
  "WpMiALa","WpMiALb","WpMiALc",
  "WpMiARa","WpMiARb","WpMiARc"
)
colData_WpMA <- data.frame(sampleID_WpMA, type_LR)
dds_WpMA <- DESeqDataSetFromTximport(txi_WpMA, colData_WpMA, design= ~ type_LR)
dds_WpMA <- DESeq(dds_WpMA)
WpMA_DESeq2 <- results(dds_WpMA)
write.table(WpMA_DESeq2, file = "WpStrMA_Wp002WmWpStrST_DESeq2")




files_WpStrES <- c(
  "quants_Wp002WmWpStrST/Wp_EaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf"
)
txi_WpStrES <- tximport(files_WpStrES, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WpStrES = c(
  "WpEaSLa","WpEaSLb","WpEaSLc",
  "WpEaSRa","WpEaSRb","WpEaSRc"
)
colData_WpStrES <- data.frame(sampleID_WpStrES, type_LR)
dds_WpStrES <- DESeqDataSetFromTximport(txi_WpStrES, colData_WpStrES, design= ~ type_LR)
dds_WpStrES <- DESeq(dds_WpStrES)
WpStrES_DESeq2 <- results(dds_WpStrES)
write.table(WpStrES_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WpStrES_Wp002WmWpStrST_DESeq2")


files_WpStrEA <- c(
  "quants_Wp002WmWpStrST/Wp_EaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALc_quant/quant.sf"
)
txi_WpStrEA <- tximport(files_WpStrEA, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WpStrEA = c(
  "WpEaALa","WpEaALb","WpEaALc",
  "WpEaARa","WpEaARb","WpEaARc"
)
colData_WpStrEA <- data.frame(sampleID_WpStrEA, type_LR)
dds_WpStrEA <- DESeqDataSetFromTximport(txi_WpStrEA, colData_WpStrEA, design= ~ type_LR)
dds_WpStrEA <- DESeq(dds_WpStrEA)
WpStrEA_DESeq2 <- results(dds_WpStrEA)
write.table(WpStrEA_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WpStrEA_Wp002WmWpStrST_DESeq2")


files_WpStrMS <- c(
  "quants_Wp002WmWpStrST/Wp_MiSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLc_quant/quant.sf"
)
txi_WpStrMS <- tximport(files_WpStrMS, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WpStrMS = c(
  "WpMiSLa","WpMiSLb","WpMiSLc",
  "WpMiSRa","WpMiSRb","WpMiSRc"
)
colData_WpStrMS <- data.frame(sampleID_WpStrMS, type_LR)
dds_WpStrMS <- DESeqDataSetFromTximport(txi_WpStrMS, colData_WpStrMS, design= ~ type_LR)
dds_WpStrMS <- DESeq(dds_WpStrMS)
WpStrMS_DESeq2 <- results(dds_WpStrMS)
write.table(WpStrMS_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WpStrMS_Wp002WmWpStrST_DESeq2")

files_WpStrLS <- c(
  "quants_Wp002WmWpStrST/Wp_LaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLc_quant/quant.sf"
)
txi_WpStrLS <- tximport(files_WpStrLS, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WpStrLS = c(
  "WpLaSLa","WpLaSLb","WpLaSLc",
  "WpLaSRa","WpLaSRb","WpLaSRc"
)
colData_WpStrLS <- data.frame(sampleID_WpStrLS, type_LR)
dds_WpStrLS <- DESeqDataSetFromTximport(txi_WpStrLS, colData_WpStrLS, design= ~ type_LR)
dds_WpStrLS <- DESeq(dds_WpStrLS)
WpStrLS_DESeq2 <- results(dds_WpStrLS)
write.table(WpStrLS_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WpStrLS_Wp002WmWpStrST_DESeq2")

files_WpStrLA <- c(
  "quants_Wp002WmWpStrST/Wp_LaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALc_quant/quant.sf"
)
txi_WpStrLA <- tximport(files_WpStrLA, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WpStrLA = c(
  "WpLaALa","WpLaALb","WpLaALc",
  "WpLaARa","WpLaARb","WpLaARc"
)
colData_WpStrLA <- data.frame(sampleID_WpStrLA, type_LR)
dds_WpStrLA <- DESeqDataSetFromTximport(txi_WpStrLA, colData_WpStrLA, design= ~ type_LR)
dds_WpStrLA <- DESeq(dds_WpStrLA)
WpStrLA_DESeq2 <- results(dds_WpStrLA)
write.table(WpStrLA_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WpStrLA_Wp002WmWpStrST_DESeq2")


files_WpStr <- c(
  "quants_Wp002WmWpStrST/Wp_EaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSLa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSLc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaALa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaALc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaSRc_quant/quant.sf",   
  "quants_Wp002WmWpStrST/Wp_MiSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaSRa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaSRc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_EaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_EaARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_MiARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_MiARc_quant/quant.sf",
  "quants_Wp002WmWpStrST/Wp_LaARa_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARb_quant/quant.sf","quants_Wp002WmWpStrST/Wp_LaARc_quant/quant.sf")

Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

txi_WpStr <- tximport(files_WpStr, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)



##DGE analysis for Wm

Wp002WmWpStrST.tx2gene <- read.table("Wp002WmWpStrST.tx2gene")

type_LR = c("L","L","L","R","R","R")

files_WmES <- c(
  "quants_Wp002WmWpStrST/L_ES_A_04_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_B_05_quant/quant.sf","quants_Wp002WmWpStrST/L_ES_C_06_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_ES_A_01_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_B_02_quant/quant.sf","quants_Wp002WmWpStrST/R_ES_C_03_quant/quant.sf"
)
txi_WmES <- tximport(files_WmES, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WmES = c(
  "WmEaSLa","WmEaSLb","WmEaSLc",
  "WmEaSRa","WmEaSRb","WmEaSRc"
)
colData_WmES <- data.frame(sampleID_WmES, type_LR)
dds_WmES <- DESeqDataSetFromTximport(txi_WmES, colData_WmES, design= ~ type_LR)
dds_WmES <- DESeq(dds_WmES)
WmES_DESeq2 <- results(dds_WmES)
write.table(WmES_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WmES_Wp002WmWpStrST_DESeq2")

files_WmEA <- c(
  "quants_Wp002WmWpStrST/L_EA_A_10_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_B_11_quant/quant.sf","quants_Wp002WmWpStrST/L_EA_C_12_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_EA_A_07_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_B_08_quant/quant.sf","quants_Wp002WmWpStrST/R_EA_C_09_quant/quant.sf"
)
txi_WmEA <- tximport(files_WmEA, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WmEA = c(
  "WmEaALa","WmEaALb","WmEaALc",
  "WmEaARa","WmEaARb","WmEaARc"
)
colData_WmEA <- data.frame(sampleID_WmEA, type_LR)
dds_WmEA <- DESeqDataSetFromTximport(txi_WmEA, colData_WmEA, design= ~ type_LR)
dds_WmEA <- DESeq(dds_WmEA)
WmEA_DESeq2 <- results(dds_WmEA)
write.table(WmEA_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WmEA_Wp002WmWpStrST_DESeq2")


files_WmMS <- c(
  "quants_Wp002WmWpStrST/L_MS_A_16_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_B_17_quant/quant.sf","quants_Wp002WmWpStrST/L_MS_C_18_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_MS_A_13_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_B_14_quant/quant.sf","quants_Wp002WmWpStrST/R_MS_C_15_quant/quant.sf"
)
txi_WmMS <- tximport(files_WmMS, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WmMS = c(
  "WmMiSLa","WmMiSLb","WmMiSLc",
  "WmMiSRa","WmMiSRb","WmMiSRc"
)
colData_WmMS <- data.frame(sampleID_WmMS, type_LR)
dds_WmMS <- DESeqDataSetFromTximport(txi_WmMS, colData_WmMS, design= ~ type_LR)
dds_WmMS <- DESeq(dds_WmMS)
WmMS_DESeq2 <- results(dds_WmMS)
write.table(WmMS_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WmMS_Wp002WmWpStrST_DESeq2")

files_WmMA <- c(
  "quants_Wp002WmWpStrST/L_MA_A_22_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_B_23_quant/quant.sf","quants_Wp002WmWpStrST/L_MA_C_24_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_MA_A_19_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_B_20_quant/quant.sf","quants_Wp002WmWpStrST/R_MA_C_21_quant/quant.sf"
)
txi_WmMA <- tximport(files_WmMA, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WmMA = c(
  "WmMiALa","WmMiALb","WmMiALc",
  "WmMiARa","WmMiARb","WmMiARc"
)
colData_WmMA <- data.frame(sampleID_WmMA, type_LR)
dds_WmMA <- DESeqDataSetFromTximport(txi_WmMA, colData_WmMA, design= ~ type_LR)
dds_WmMA <- DESeq(dds_WmMA)
WmMA_DESeq2 <- results(dds_WmMA)
write.table(WmMA_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WmMA_Wp002WmWpStrST_DESeq2")

files_WmLS <- c(
  "quants_Wp002WmWpStrST/L_LS_A_28_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_B_29_quant/quant.sf","quants_Wp002WmWpStrST/L_LS_C_30_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_LS_A_25_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_B_26_quant/quant.sf","quants_Wp002WmWpStrST/R_LS_C_27_quant/quant.sf"
  
)
txi_WmLS <- tximport(files_WmLS, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WmLS = c(
  "WmLaSLa","WmLaSLb","WmLaSLc",
  "WmLaSRa","WmLaSRb","WmLaSRc"
)
colData_WmLS <- data.frame(sampleID_WmLS, type_LR)
dds_WmLS <- DESeqDataSetFromTximport(txi_WmLS, colData_WmLS, design= ~ type_LR)
dds_WmLS <- DESeq(dds_WmLS)
WmLS_DESeq2 <- results(dds_WmLS)
write.table(WmLS_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WmLS_Wp002WmWpStrST_DESeq2")

files_WmLA <- c(
  "quants_Wp002WmWpStrST/L_LA_A_34_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_B_35_quant/quant.sf","quants_Wp002WmWpStrST/L_LA_C_36_quant/quant.sf",
  "quants_Wp002WmWpStrST/R_LA_A_31_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_B_32_quant/quant.sf","quants_Wp002WmWpStrST/R_LA_C_33_quant/quant.sf"
)
txi_WmLA <- tximport(files_WmLA, type = "salmon", tx2gene = Wp002WmWpStrST.tx2gene)
sampleID_WmLA = c(
  "WmLaALa","WmLaALb","WmLaALc",
  "WmLaARa","WmLaARb","WmLaARc"
)
colData_WmLA <- data.frame(sampleID_WmLA, type_LR)
dds_WmLA <- DESeqDataSetFromTximport(txi_WmLA, colData_WmLA, design= ~ type_LR)
dds_WmLA <- DESeq(dds_WmLA)
WmLA_DESeq2 <- results(dds_WmLA)
write.table(WmLA_DESeq2, file = "DESeq2 results with Mercator4 withSAURs/WmLA_Wp002WmWpStrST_DESeq2")

png(file="WpStr.Wp002WmWpStrST.30480(QRT2-1).png", width=400, height=400)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.30480", main="Wp002WmWpStrST.30480(QRT2-1)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.50488(QRT2-2).png", width=400, height=400)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.50488", main="Wp002WmWpStrST.50488(QRT2-2)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.50507(QRT2-3).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.50507", main="Wp002WmWpStrST.50507(QRT2-3)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.56810(QRT2-4).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.56810", main="Wp002WmWpStrST.56810(QRT2-4)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.13087(BGAL1-1).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.13087", main="Wp002WmWpStrST.13087(BGAL1-1)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.13090(BGAL1-2).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.13090", main="Wp002WmWpStrST.13090(BGAL1-2)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.24613(BGAL1-3).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.24613", main="Wp002WmWpStrST.24613(BGAL1-3)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.26866(BGAL1-4).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.26866", main="Wp002WmWpStrST.26866(BGAL1-4)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.35261(BGAL1-5).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.35261", main="Wp002WmWpStrST.35261(BGAL1-5)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.39690(BGAL1-6).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.39690", main="Wp002WmWpStrST.39690(BGAL1-6)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.45788(BGAL1-7).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.45788", main="Wp002WmWpStrST.45788(BGAL1-7)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.52551(BGAL1-8).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.52551", main="Wp002WmWpStrST.52551(BGAL1-8)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="WpStr.Wp002WmWpStrST.54164(BGAL1-9).png", width=2400, height=800)
plotCounts(dds_WpStr, gene="Wp002WmWpStrST.54164", main="Wp002WmWpStrST.54164(BGAL1-9)", transform = FALSE, intgroup = "type_WpStr", returnData = FALSE)
dev.off()

png(file="Wm.Wp002WmWpStrST.30480(QRT2-1).png", width=400, height=400)
Wp002WmWpStrST.30480 <- plotCounts(dds_Wm, gene="Wp002WmWpStrST.30480", main="Wp002WmWpStrST.30480(QRT2-1)", transform = FALSE, intgroup = "type_Wm", returnData = FALSE)
dev.off()

png(file="Wm.Wp002WmWpStrST.50488(QRT2-2).png", width=400, height=400)
plotCounts(dds_Wm, gene="Wp002WmWpStrST.50488", main="Wp002WmWpStrST.50488(QRT2-2)", transform = FALSE, intgroup = "type_Wm", returnData = FALSE)
dev.off()

png(file="Wm.Wp002WmWpStrST.50507(QRT2-3).png", width=400, height=400)
plotCounts(dds_Wm, gene="Wp002WmWpStrST.50507", main="Wp002WmWpStrST.50507(QRT2-3)", transform = FALSE, intgroup = "type_Wm", returnData = FALSE)
dev.off()

png(file="Wm.Wp002WmWpStrST.56810(QRT2-4).png", width=400, height=400)
plotCounts(dds_Wm, gene="Wp002WmWpStrST.56810", main="Wp002WmWpStrST.56810(QRT2-4)", transform = FALSE, intgroup = "type_Wm", returnData = FALSE)
dev.off()

