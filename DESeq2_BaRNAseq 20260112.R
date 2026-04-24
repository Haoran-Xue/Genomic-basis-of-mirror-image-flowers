if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")

### DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

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

setwd("D:/Michael Lenhard - HFSP Shared folder/Universität Potsdam/Michael Lenhard - HFSP Shared folder/Wachendorfia/Barberetta RNA-seq/")
setwd("G:/Universität Potsdam/Michael Lenhard - HFSP Shared folder/Wachendorfia/Barberetta RNA-seq")

files_all <- c(
  "quants/Ba2024_LeAnLa1_val_quant/quant.sf",
  "quants/Ba2024_LeAnLa2_val_quant/quant.sf",
  "quants/Ba2024_LeAnLa3_val_quant/quant.sf",
  "quants/Ba2024_LeAnMi1_val_quant/quant.sf",
  "quants/Ba2024_LeAnMi2_val_quant/quant.sf",
  "quants/Ba2024_LeAnMi4_val_quant/quant.sf",
  "quants/Ba2024_LeStLa1_val_quant/quant.sf",
  "quants/Ba2024_LeStLa4_val_quant/quant.sf",
  "quants/Ba2024_LeStLa5_val_quant/quant.sf",
  "quants/Ba2024_LeStMi1_val_quant/quant.sf",
  "quants/Ba2024_LeStMi2_val_quant/quant.sf",
  "quants/Ba2024_LeStMi4_val_quant/quant.sf",
  "quants/Ba2024_RiAnLa1_val_quant/quant.sf",
  "quants/Ba2024_RiAnLa2_val_quant/quant.sf",
  "quants/Ba2024_RiAnLa4_val_quant/quant.sf",
  "quants/Ba2024_RiAnMi1_val_quant/quant.sf",
  "quants/Ba2024_RiAnMi2_val_quant/quant.sf",
  "quants/Ba2024_RiAnMi4_val_quant/quant.sf",
  "quants/Ba2024_RiStLa1_val_quant/quant.sf",
  "quants/Ba2024_RiStLa4_val_quant/quant.sf",
  "quants/Ba2024_RiStLa5_val_quant/quant.sf",
  "quants/Ba2024_RiStMi1_val_quant/quant.sf",
  "quants/Ba2024_RiStMi3_val_quant/quant.sf",
  "quants/Ba2024_RiStMi4_val_quant/quant.sf",
  "quants/Ba2025_L_EA1_val_quant/quant.sf",
  "quants/Ba2025_L_EA2_val_quant/quant.sf",
  "quants/Ba2025_L_EA3_val_quant/quant.sf",
  "quants/Ba2025_L_EA4_val_quant/quant.sf",
  "quants/Ba2025_L_ES1_val_quant/quant.sf",
  "quants/Ba2025_L_ES2_val_quant/quant.sf",
  "quants/Ba2025_L_ES3_val_quant/quant.sf",
  "quants/Ba2025_L_ES4_val_quant/quant.sf",
  "quants/Ba2025_L_LA1_val_quant/quant.sf",
  "quants/Ba2025_L_LA2_val_quant/quant.sf",
  "quants/Ba2025_L_LA3_val_quant/quant.sf",
  "quants/Ba2025_L_LA4_val_quant/quant.sf",
  "quants/Ba2025_L_LS1_val_quant/quant.sf",
  "quants/Ba2025_L_LS2_val_quant/quant.sf",
  "quants/Ba2025_L_LS3_val_quant/quant.sf",
  "quants/Ba2025_L_LS4_val_quant/quant.sf",
  "quants/Ba2025_L_MA1_val_quant/quant.sf",
  "quants/Ba2025_L_MA2_val_quant/quant.sf",
  "quants/Ba2025_L_MA3_val_quant/quant.sf",
  "quants/Ba2025_L_MA4_val_quant/quant.sf",
  "quants/Ba2025_L_MS1_val_quant/quant.sf",
  "quants/Ba2025_L_MS2_val_quant/quant.sf",
  "quants/Ba2025_L_MS3_val_quant/quant.sf",
  "quants/Ba2025_L_MS4_val_quant/quant.sf",
  "quants/Ba2025_R_EA1_val_quant/quant.sf",
  "quants/Ba2025_R_EA2_val_quant/quant.sf",
  "quants/Ba2025_R_EA3_val_quant/quant.sf",
  "quants/Ba2025_R_EA4_val_quant/quant.sf",
  "quants/Ba2025_R_ES1_val_quant/quant.sf",
  "quants/Ba2025_R_ES2_val_quant/quant.sf",
  "quants/Ba2025_R_ES3_val_quant/quant.sf",
  "quants/Ba2025_R_ES4_val_quant/quant.sf",
  "quants/Ba2025_R_LA1_val_quant/quant.sf",
  "quants/Ba2025_R_LA2_val_quant/quant.sf",
  "quants/Ba2025_R_LA3_val_quant/quant.sf",
  "quants/Ba2025_R_LS1_val_quant/quant.sf",
  "quants/Ba2025_R_LS2_val_quant/quant.sf",
  "quants/Ba2025_R_LS3_val_quant/quant.sf",
  "quants/Ba2025_R_LS4_val_quant/quant.sf",
  "quants/Ba2025_R_MA1_val_quant/quant.sf",
  "quants/Ba2025_R_MA2_val_quant/quant.sf",
  "quants/Ba2025_R_MA3_val_quant/quant.sf",
  "quants/Ba2025_R_MA4_val_quant/quant.sf",
  "quants/Ba2025_R_MS1_val_quant/quant.sf",
  "quants/Ba2025_R_MS2_val_quant/quant.sf",
  "quants/Ba2025_R_MS3_val_quant/quant.sf",
  "quants/Ba2025_R_MS4_val_quant/quant.sf",
  "quants/F01_R_val_quant/quant.sf",
  "quants/F02_R_val_quant/quant.sf",
  "quants/F03_R_val_quant/quant.sf",
  "quants/F04_R_val_quant/quant.sf",
  "quants/F10_R_val_quant/quant.sf",
  "quants/F06_L_val_quant/quant.sf",
  "quants/F07_L_val_quant/quant.sf",
  "quants/F08_L_val_quant/quant.sf",
  "quants/F09_L_val_quant/quant.sf",
  "quants/F11_L_val_quant/quant.sf",
  "quants/S01_R_val_quant/quant.sf",
  "quants/S02_R_val_quant/quant.sf",
  "quants/S03_R_val_quant/quant.sf",
  "quants/S04_R_val_quant/quant.sf",
  "quants/S10_R_val_quant/quant.sf",
  "quants/S05_L_val_quant/quant.sf",
  "quants/S06_L_val_quant/quant.sf",
  "quants/S08_L_val_quant/quant.sf",
  "quants/S09_L_val_quant/quant.sf",
  "quants/S11_L_val_quant/quant.sf")

BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")

txi <- tximport(files_all, type = "salmon", tx2gene = BaN.ST.tx2gene)

sampleID = c(
  "Ba2024_LeAnLa1","Ba2024_LeAnLa2","Ba2024_LeAnLa3","Ba2024_LeAnMi1","Ba2024_LeAnMi2","Ba2024_LeAnMi4",
  "Ba2024_LeStLa1","Ba2024_LeStLa4","Ba2024_LeStLa5","Ba2024_LeStMi1","Ba2024_LeStMi2","Ba2024_LeStMi4",
  "Ba2024_RiAnLa1","Ba2024_RiAnLa2","Ba2024_RiAnLa4","Ba2024_RiAnMi1","Ba2024_RiAnMi2","Ba2024_RiAnMi4",
  "Ba2024_RiStLa1","Ba2024_RiStLa4","Ba2024_RiStLa5","Ba2024_RiStMi1","Ba2024_RiStMi3","Ba2024_RiStMi4",
  "Ba2025_L_EA1","Ba2025_L_EA2","Ba2025_L_EA3","Ba2025_L_EA4","Ba2025_L_ES1","Ba2025_L_ES2","Ba2025_L_ES3","Ba2025_L_ES4",
  "Ba2025_L_LA1","Ba2025_L_LA2","Ba2025_L_LA3","Ba2025_L_LA4","Ba2025_L_LS1","Ba2025_L_LS2","Ba2025_L_LS3","Ba2025_L_LS4",
  "Ba2025_L_MA1","Ba2025_L_MA2","Ba2025_L_MA3","Ba2025_L_MA4","Ba2025_L_MS1","Ba2025_L_MS2","Ba2025_L_MS3","Ba2025_L_MS4",
  "Ba2025_R_EA1","Ba2025_R_EA2","Ba2025_R_EA3","Ba2025_R_EA4","Ba2025_R_ES1","Ba2025_R_ES2","Ba2025_R_ES3","Ba2025_R_ES4",
  "Ba2025_R_LA1","Ba2025_R_LA2","Ba2025_R_LA3","Ba2025_R_LS1","Ba2025_R_LS2","Ba2025_R_LS3","Ba2025_R_LS4",
  "Ba2025_R_MA1","Ba2025_R_MA2","Ba2025_R_MA3","Ba2025_R_MA4","Ba2025_R_MS1","Ba2025_R_MS2","Ba2025_R_MS3","Ba2025_R_MS4",
  "Ba2026_F01_R","Ba2026_F02_R","Ba2026_F03_R","Ba2026_F04_R","Ba2026_F10_R","Ba2026_F06_L","Ba2026_F07_L","Ba2026_F08_L","Ba2026_F09_L","Ba2026_F11_L",
  "Ba2026_S01_R","Ba2026_S02_R","Ba2026_S03_R","Ba2026_S04_R","Ba2026_S10_R","Ba2026_S05_L","Ba2026_S06_L","Ba2026_S08_L","Ba2026_S09_L","Ba2026_S11_L"
)

type = c(
  "LeSmnLa","LeSmnLa","LeSmnLa","LeSmnMi","LeSmnMi","LeSmnMi",
  "LeStyLa","LeStyLa","LeStyLa","LeStyMi","LeStyMi","LeStyMi",
  "RiSmnLa","RiSmnLa","RiSmnLa","RiSmnMi","RiSmnMi","RiSmnMi",
  "RiStyLa","RiStyLa","RiStyLa","RiStyMi","RiStyMi","RiStyMi",
  "LeSmnEa","LeSmnE","LeSmnE","LeSmnE","LeStyEa","LeStyEa","LeStyEa","LeStyEa",
  "LeSmnLa","LeSmnL","LeSmnL","LeSmnL","LeStyLa","LeStyLa","LeStyLa","LeStyLa",
  "LeSmnM","LeSmnM","LeSmnM","LeSmnM","LeStyMi","LeStyMi","LeStyMi","LeStyMi",
  "RiSmnEa","RiSmnE","RiSmnE","RiSmnE","RiStyEa","RiStyEa","RiStyEa","RiStyEa",
  "RiSmnLa","RiSmnLa","RiSmnLa","RiStyLa","RiStyLa","RiStyLa","RiStyLa",
  "RiSmnMi","RiSmnMi","RiSmnMi","RiSmnMi","RiStyMi","RiStyMi","RiStyMi","RiStyMi",
  "RiFilVE","RiFilVE","RiFilVE","RiFilVE","RiFilVE","LeFilVE","LeFilVE","LeFilVE","LeFilVE","LeFilVE",
  "RiStyVE","RiStyVE","RiStyVE","RiStyVE","RiStyVE","LeStyVE","LeStyVE","LeStyVE","LeStyVE","LeStyVE"
  )



colData <- data.frame(sampleID, type)
dds <- DESeqDataSetFromTximport(txi, colData, design= ~ type)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
colnames(normalized_counts) <- sampleID
write.csv(normalized_counts, file = "Ba_normalized_gene_count_matrix.csv")


png(file="BaNST40750.YUC-R.png", width=1800, height=1500)
plotCounts(dds, gene="BaNST40750", transform = FALSE, main = "BaNST40750.YUC-R", intgroup = "type", returnData = FALSE)
dev.off()


png(file="BaNST17896.YUC-P1.png", width=1800, height=1500)
plotCounts(dds, gene="BaNST17896", transform = FALSE, main = "BaNST17896.YUC-P1", intgroup = "type", returnData = FALSE)
dev.off()

##PCA
vst <- vst(dds, blind=FALSE)

vst_matrix <- getVarianceStabilizedData(dds)

z <- DESeq2::plotPCA(vst, intgroup = "sampleID", ntop = 500, returnData = FALSE)

nudge <- position_nudge(y = 1)

png(file="c RNA-seq.png", width=1200, height=600)
z + geom_text(aes(label = sampleID), position = nudge) + ggtitle("sampleID RNA-seq") +
  theme_bw()
dev.off()



files_2026 <- c(
  "quants/F01_R_val_quant/quant.sf",
  "quants/F02_R_val_quant/quant.sf",
  "quants/F03_R_val_quant/quant.sf",
  "quants/F04_R_val_quant/quant.sf",
  "quants/F10_R_val_quant/quant.sf",
  "quants/F06_L_val_quant/quant.sf",
  "quants/F07_L_val_quant/quant.sf",
  "quants/F08_L_val_quant/quant.sf",
  "quants/F09_L_val_quant/quant.sf",
  "quants/F11_L_val_quant/quant.sf",
  "quants/S01_R_val_quant/quant.sf",
  "quants/S02_R_val_quant/quant.sf",
  "quants/S03_R_val_quant/quant.sf",
  "quants/S04_R_val_quant/quant.sf",
  "quants/S10_R_val_quant/quant.sf",
  "quants/S05_L_val_quant/quant.sf",
  "quants/S06_L_val_quant/quant.sf",
  "quants/S08_L_val_quant/quant.sf",
  "quants/S09_L_val_quant/quant.sf",
  "quants/S11_L_val_quant/quant.sf")

BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")

txi_2026 <- tximport(files_2026, type = "salmon", tx2gene = BaN.ST.tx2gene)

sampleID_2026 = c(
  "F01_R","F02_R","F03_R","F04_R","F10_R","F06_L","F07_L","F08_L","F09_L","F11_L",
  "S01_R","S02_R","S03_R","S04_R","S10_R","S05_L","S06_L","S08_L","S09_L","S11_L"
)

type_2026 = c(
  "RiFilVE","RiFilVE","RiFilVE","RiFilVE","RiFilVE","LeFilVE","LeFilVE","LeFilVE","LeFilVE","LeFilVE",
  "RiStyVE","RiStyVE","RiStyVE","RiStyVE","RiStyVE","LeStyVE","LeStyVE","LeStyVE","LeStyVE","LeStyVE"
)



colData_2026 <- data.frame(sampleID_2026, type_2026)
dds_2026 <- DESeqDataSetFromTximport(txi_2026, colData_2026, design= ~ type_2026)

dds_2026 <- estimateSizeFactors(dds_2026)
normalized_counts_2026 <- counts(dds_2026, normalized=TRUE)
colnames(normalized_counts_2026) <- sampleID_2026
write.csv(normalized_counts_2026, file = "Ba2026_normalized_gene_count_matrix.csv")

##PCA
vst_2026 <- vst(dds_2026, blind=FALSE)

vst_matrix <- getVarianceStabilizedData(dds)

z <- DESeq2::plotPCA(vst_2026, intgroup = "sampleID_2026", ntop = 500, returnData = FALSE)

nudge <- position_nudge(y = 1)

png(file="Ba2026 RNA-seq.png", width=1200, height=600)
z + geom_text(aes(label = sampleID_2026), position = nudge) + ggtitle("Ba2026 RNA-seq") +
  theme_bw()
dev.off()




filesVS <- c(
    "quants/S05_L_val_quant/quant.sf",
    "quants/S06_L_val_quant/quant.sf",
    "quants/S08_L_val_quant/quant.sf",
    "quants/S09_L_val_quant/quant.sf",
    "quants/S11_L_val_quant/quant.sf",
    "quants/S01_R_val_quant/quant.sf",
    "quants/S02_R_val_quant/quant.sf",
    "quants/S03_R_val_quant/quant.sf",
    "quants/S04_R_val_quant/quant.sf",
    "quants/S10_R_val_quant/quant.sf")

txiVS <- tximport(filesVS, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDVS = c(
  "S05_L","S06_L","S08_L","S09_L","S11_L","S01_R","S02_R","S03_R","S04_R","S10_R"
)
typeVLR = c("L","L","L","L","L","R","R","R","R","R")

colDataVS <- data.frame(sampleIDVS, typeVLR)
ddsVS <- DESeqDataSetFromTximport(txiVS, colDataVS, design= ~ typeVLR)
ddsVS <- DESeq(ddsVS)
BaVS_DESeq2 <- results(ddsVS)
write.table(BaVS_DESeq2, file = "BaVS_DESeq2.txt")


filesVF <- c(
  "quants/F06_L_val_quant/quant.sf",
  "quants/F07_L_val_quant/quant.sf",
  "quants/F08_L_val_quant/quant.sf",
  "quants/F09_L_val_quant/quant.sf",
  "quants/F11_L_val_quant/quant.sf",
  "quants/F01_R_val_quant/quant.sf",
  "quants/F02_R_val_quant/quant.sf",
  "quants/F03_R_val_quant/quant.sf",
  "quants/F04_R_val_quant/quant.sf",
  "quants/F10_R_val_quant/quant.sf")

txiVF <- tximport(filesVF, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDVF = c(
  "F06_L","F07_L","F08_L","F09_L","F11_L","F01_R","F02_R","F03_R","F04_R","F10_R"
)
typeVLR = c("L","L","L","L","L","R","R","R","R","R")

colDataVF <- data.frame(sampleIDVF, typeVLR)
ddsVF <- DESeqDataSetFromTximport(txiVF, colDataVF, design= ~ typeVLR)
ddsVF <- DESeq(ddsVF)
BaVF_DESeq2 <- results(ddsVF)
write.table(BaVF_DESeq2, file = "BaVF_DESeq2.txt", sep = "\t")


png(file=".png", width=900, height=750)

plotCounts(ddsVF, gene="BaNST40754", transform = FALSE, main = "YUC-R", intgroup = "typeVLR", returnData = FALSE)

plotCounts(ddsVS, gene="BaNST40754", transform = FALSE, main = "YUC-R", intgroup = "typeVLR", returnData = FALSE)


png(file="BaNST12945.SPL12.png", width=900, height=750)
plotCounts(ddsVS, gene="BaNST12945", transform = FALSE, main = "BaNST12945.SPL12", intgroup = "typeVLR", returnData = FALSE)



dev.off()







filesES <- c(
  "quants/Ba2025_L_ES1_val_quant/quant.sf",
  "quants/Ba2025_L_ES2_val_quant/quant.sf",
  "quants/Ba2025_L_ES3_val_quant/quant.sf",
  "quants/Ba2025_L_ES4_val_quant/quant.sf",
  "quants/Ba2025_R_ES1_val_quant/quant.sf",
  "quants/Ba2025_R_ES2_val_quant/quant.sf",
  "quants/Ba2025_R_ES3_val_quant/quant.sf",
  "quants/Ba2025_R_ES4_val_quant/quant.sf")

BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")
txiES <- tximport(filesES, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDES = c(
  "ES1_L","ES2_L","ES3_L","ES4_L","ES1_R","ES2_R","ES3_R","ES4_R"
)

typeELR = c("L","L","L","L","R","R","R","R")

colDataES <- data.frame(sampleIDES, typeELR)
ddsES <- DESeqDataSetFromTximport(txiES, colDataES, design= ~ typeELR)
ddsES <- DESeq(ddsES)
BaES_DESeq2 <- results(ddsES)
write.table(BaES_DESeq2, file = "BaES_DESeq2.txt")



filesMS <- c(
  "quants/Ba2024_LeStMi1_val_quant/quant.sf",
  "quants/Ba2024_LeStMi2_val_quant/quant.sf",
  "quants/Ba2024_LeStMi4_val_quant/quant.sf",
  "quants/Ba2025_L_MS1_val_quant/quant.sf",
  "quants/Ba2025_L_MS2_val_quant/quant.sf",
  "quants/Ba2025_L_MS3_val_quant/quant.sf",
  "quants/Ba2025_L_MS4_val_quant/quant.sf",
  "quants/Ba2024_RiStMi1_val_quant/quant.sf",
  "quants/Ba2024_RiStMi3_val_quant/quant.sf",
  "quants/Ba2024_RiStMi4_val_quant/quant.sf",
  "quants/Ba2025_R_MS1_val_quant/quant.sf",
  "quants/Ba2025_R_MS2_val_quant/quant.sf",
  "quants/Ba2025_R_MS3_val_quant/quant.sf",
  "quants/Ba2025_R_MS4_val_quant/quant.sf")

BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")
txiMS <- tximport(filesMS, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDMS = c(
  "MS1_L","MS2_L","MS3_L","MS4_L","MS5_L","MS6_L","MS7_L","MS1_R","MS2_R","MS3_R","MS4_R","MS5_R","MS6_R","MS7_R"
)

typeMLLR = c("L","L","L","L","L","L","L","R","R","R","R","R","R","R")

colDataMS <- data.frame(sampleIDMS, typeMLLR)
ddsMS <- DESeqDataSetFromTximport(txiMS, colDataMS, design= ~ typeMLLR)
ddsMS <- DESeq(ddsMS)
BaMS_DESeq2 <- results(ddsMS)
write.table(BaMS_DESeq2, file = "BaMS_DESeq2.txt")




filesLS <- c(
  "quants/Ba2024_LeStLa1_val_quant/quant.sf",
  "quants/Ba2024_LeStLa4_val_quant/quant.sf",
  "quants/Ba2024_LeStLa5_val_quant/quant.sf",
  "quants/Ba2025_L_LS1_val_quant/quant.sf",
  "quants/Ba2025_L_LS2_val_quant/quant.sf",
  "quants/Ba2025_L_LS3_val_quant/quant.sf",
  "quants/Ba2025_L_LS4_val_quant/quant.sf",
  "quants/Ba2024_RiStLa1_val_quant/quant.sf",
  "quants/Ba2024_RiStLa4_val_quant/quant.sf",
  "quants/Ba2024_RiStLa5_val_quant/quant.sf",
  "quants/Ba2025_R_LS1_val_quant/quant.sf",
  "quants/Ba2025_R_LS2_val_quant/quant.sf",
  "quants/Ba2025_R_LS3_val_quant/quant.sf",
  "quants/Ba2025_R_LS4_val_quant/quant.sf")


BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")
txiLS <- tximport(filesLS, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDLS = c(
  "LS1_L","LS2_L","LS3_L","LS4_L","LS5_L","LS6_L","LS7_L","LS1_R","LS2_R","LS3_R","LS4_R","LS5_R","LS6_R","LS7_R"
)

typeMLLR = c("L","L","L","L","L","L","L","R","R","R","R","R","R","R")

colDataLS <- data.frame(sampleIDLS, typeMLLR)
ddsLS <- DESeqDataSetFromTximport(txiLS, colDataLS, design= ~ typeMLLR)
ddsLS <- DESeq(ddsLS)
BaLS_DESeq2 <- results(ddsLS)
write.table(BaLS_DESeq2, file = "BaLS_DESeq2.txt")




filesEA <- c(
  "quants/Ba2025_L_EA1_val_quant/quant.sf",
  "quants/Ba2025_L_EA2_val_quant/quant.sf",
  "quants/Ba2025_L_EA3_val_quant/quant.sf",
  "quants/Ba2025_L_EA4_val_quant/quant.sf",
  "quants/Ba2025_R_EA1_val_quant/quant.sf",
  "quants/Ba2025_R_EA2_val_quant/quant.sf",
  "quants/Ba2025_R_EA3_val_quant/quant.sf",
  "quants/Ba2025_R_EA4_val_quant/quant.sf")

BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")
txiEA <- tximport(filesEA, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDEA = c(
  "EA1_L","EA2_L","EA3_L","EA4_L","EA1_R","EA2_R","EA3_R","EA4_R"
)

typeELR = c("L","L","L","L","R","R","R","R")

colDataEA <- data.frame(sampleIDEA, typeELR)
ddsEA <- DESeqDataSetFromTximport(txiEA, colDataEA, design= ~ typeELR)
ddsEA <- DESeq(ddsEA)
BaEA_DESeq2 <- results(ddsEA)
write.table(BaEA_DESeq2, file = "BaEA_DESeq2.txt")



filesMA <- c(
  "quants/Ba2024_LeAnMi1_val_quant/quant.sf",
  "quants/Ba2024_LeAnMi2_val_quant/quant.sf",
  "quants/Ba2024_LeAnMi4_val_quant/quant.sf",
  "quants/Ba2025_L_MA1_val_quant/quant.sf",
  "quants/Ba2025_L_MA2_val_quant/quant.sf",
  "quants/Ba2025_L_MA3_val_quant/quant.sf",
  "quants/Ba2025_L_MA4_val_quant/quant.sf",
  "quants/Ba2024_RiAnMi1_val_quant/quant.sf",
  "quants/Ba2024_RiAnMi2_val_quant/quant.sf",
  "quants/Ba2024_RiAnMi4_val_quant/quant.sf",
  "quants/Ba2025_R_MA1_val_quant/quant.sf",
  "quants/Ba2025_R_MA2_val_quant/quant.sf",
  "quants/Ba2025_R_MA3_val_quant/quant.sf",
  "quants/Ba2025_R_MA4_val_quant/quant.sf")

BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")
txiMA <- tximport(filesMA, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDMA = c(
  "MA1_L","MA2_L","MA3_L","MA4_L","MA5_L","MA6_L","MA7_L","MA1_R","MA2_R","MA3_R","MA4_R","MA5_R","MA6_R","MA7_R"
)

typeMLLR = c("L","L","L","L","L","L","L","R","R","R","R","R","R","R")

colDataMA <- data.frame(sampleIDMA, typeMLLR)
ddsMA <- DESeqDataSetFromTximport(txiMA, colDataMA, design= ~ typeMLLR)
ddsMA <- DESeq(ddsMA)
BaMA_DESeq2 <- results(ddsMA)
write.table(BaMA_DESeq2, file = "BaMA_DESeq2.txt")




filesLA <- c(
  "quants/Ba2024_LeAnLa1_val_quant/quant.sf",
  "quants/Ba2024_LeAnLa2_val_quant/quant.sf",
  "quants/Ba2024_LeAnLa3_val_quant/quant.sf",
  "quants/Ba2025_L_LA1_val_quant/quant.sf",
  "quants/Ba2025_L_LA2_val_quant/quant.sf",
  "quants/Ba2025_L_LA3_val_quant/quant.sf",
  "quants/Ba2025_L_LA4_val_quant/quant.sf",
  "quants/Ba2024_RiAnLa1_val_quant/quant.sf",
  "quants/Ba2024_RiAnLa2_val_quant/quant.sf",
  "quants/Ba2024_RiAnLa4_val_quant/quant.sf",
  "quants/Ba2025_R_LA1_val_quant/quant.sf",
  "quants/Ba2025_R_LA2_val_quant/quant.sf",
  "quants/Ba2025_R_LA3_val_quant/quant.sf")


BaN.ST.tx2gene <- read.table("BaN.ST.tx2gene")
txiLA <- tximport(filesLA, type = "salmon", tx2gene = BaN.ST.tx2gene)
sampleIDLA = c(
  "LA1_L","LA2_L","LA3_L","LA4_L","LA5_L","LA6_L","LA7_L","LA1_R","LA2_R","LA3_R","LA4_R","LA5_R","LA6_R"
)

typeLALR = c("L","L","L","L","L","L","L","R","R","R","R","R","R")

colDataLA <- data.frame(sampleIDLA, typeLALR)
ddsLA <- DESeqDataSetFromTximport(txiLA, colDataLA, design= ~ typeLALR)
ddsLA <- DESeq(ddsLA)
BaLA_DESeq2 <- results(ddsLA)
write.table(BaLA_DESeq2, file = "BaLA_DESeq2.txt")







res.pca.vst_stage_type <- prcomp(vst_stage_type_matrix)

library(ggfortify)
library(factoextra)

fviz_pca_ind(res.pca.vst_stage_type,
             #  col.ind = "cos2", # Color by the quality of representation
             #  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

autoplot(pc1, data = iris, colour = "Species")



# Biplot
biplot(res.pca.vst_stage_type)

# Convert to data.frame
rotation_df <- as.data.frame(res.pca.vst_stage_type$rotation)

stage_abc <- as.factor(c(
  rep(c("E","M","V"), each = 3, times = 3)
))


# Add your factor column
anther_type_abc <- as.factor(c(
  rep("C", 9),
  rep("F", 9),
  rep("C", 9)
))

# Check result
str(rotation_df)

# Color by stage
stage_cols <- c("E" = "blue", "M" = "purple", "V" = "green")

# Shape by anther type
# pch = 19 = solid circle, pch = 1 = hollow circle
anther_pch <- c("C" = 19, "F" = 1)





library(ggplot2)
library(rlang)

# --- Compute variance explained ---
pc_sdev <- res.pca.vst_stage_type$sdev
variance_explained <- pc_sdev^2 / sum(pc_sdev^2)
percent_var <- 100 * variance_explained

# --- Combined factor: species + anther type ---
rotation_df$anther_type_full <- NA
rotation_df$anther_type_full[grepl("^Caf", rotation_df$sample)] <- "Caf_connected"
rotation_df$anther_type_full[grepl("^Caa", rotation_df$sample) & anther_type_abc == "C"] <- "Caa_connected"
rotation_df$anther_type_full[grepl("^Caa", rotation_df$sample) & anther_type_abc == "F"] <- "Caa_free"

rotation_df$anther_type_full <- factor(rotation_df$anther_type_full,
                                       levels = c("Caf_connected",
                                                  "Caa_connected",
                                                  "Caa_free"))

# Shapes
anther_shapes <- c("Caf_connected" = 15,   # square solid
                   "Caa_connected" = 16,   # circle solid
                   "Caa_free" = 1)         # circle hollow

# Pretty labels
anther_labels <- c(
  expression(italic("C. alba flavescens") ~ "connected"),
  expression(italic("C. alba alba") ~ "connected"),
  expression(italic("C. alba alba") ~ "free")
)

# Stage colors
stage_cols <- c("E" = "blue", "M" = "purple", "V" = "green")

# --- PDF with plots ---
pdf("Cyanella_anther_RNA-seq_PCA_plots.pdf", width = 7, height = 5)

n_pcs <- ncol(rotation_df) - 4  # subtract metadata columns

for (i in seq(1, n_pcs, by = 2)) {
  
  pc_x <- paste0("PC", i)
  pc_y_index <- ifelse(i + 1 <= n_pcs, i + 1, i)
  pc_y <- paste0("PC", pc_y_index)
  
  # axis labels with % variance
  xlab_str <- paste0(pc_x, " (", round(percent_var[i], 1), "%)")
  ylab_str <- paste0(pc_y, " (", round(percent_var[pc_y_index], 1), "%)")
  
  p <- ggplot(rotation_df, aes(x = !!sym(pc_x), y = !!sym(pc_y),
                               color = stage_abc,
                               shape = anther_type_full,
                               label = sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.7, size = 3, show.legend = FALSE) +
    scale_color_manual(values = stage_cols) +
    scale_shape_manual(values = anther_shapes,
                       labels = anther_labels) +
    labs(x = xlab_str, y = ylab_str,
         color = "Stage", shape = "Subspecies and anther type") +
    theme_bw() +
    theme(legend.position = "right")
  
  print(p)
}

dev.off()



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










color_palette <- c(turbo(12))


png(file="Ba2025.MIR156-R.png", width=900, height=750)
plotCounts(dds, gene="MIR156-R", transform = FALSE, intgroup = "type", returnData = FALSE)
dev.off()


png(file="Ba2025.YUC-R.png", width=900, height=750)
plotCounts(dds, gene="YUC-R", transform = FALSE, intgroup = "type", returnData = FALSE)
dev.off()

BaST11047
plotCounts(dds, gene="BaST11047", transform = FALSE, intgroup = "type", returnData = FALSE)

##DGE analysis for Ba2025

filesES <- c(
  "quants/L_ES1_quant/quant.sf","quants/L_ES2_quant/quant.sf","quants/L_ES3_quant/quant.sf","quants/L_ES4_quant/quant.sf",
  "quants/R_ES1_quant/quant.sf","quants/R_ES2_quant/quant.sf","quants/R_ES3_quant/quant.sf","quants/R_ES4_quant/quant.sf"
)
txiES <- tximport(filesES, type = "salmon", tx2gene = Ba2025.ST.tx2gene)
sampleIDES = c(
  "LES1","LES2","LES3","LES4",
  "RES1","RES2","RES3","RES4"
)
type_LR = c("L","L","L","L","R","R","R","R")

colDataES <- data.frame(sampleIDES, type_LR)
ddsES <- DESeqDataSetFromTximport(txiES, colDataES, design= ~ type_LR)
ddsES <- DESeq(ddsES)
Ba2025ES_DESeq2 <- results(ddsES)
write.table(Ba2025ES_DESeq2, file = "Ba2025ES_DESeq2.txt")


filesEA <- c(
  "quants/L_EA1_quant/quant.sf","quants/L_EA2_quant/quant.sf","quants/L_EA3_quant/quant.sf","quants/L_EA4_quant/quant.sf",
  "quants/R_EA1_quant/quant.sf","quants/R_EA2_quant/quant.sf","quants/R_EA3_quant/quant.sf","quants/R_EA4_quant/quant.sf"
)
txiEA <- tximport(filesEA, type = "salmon", tx2gene = Ba2025.ST.tx2gene)
sampleIDEA = c(
  "LEA1","LEA2","LEA3","LEA4",
  "REA1","REA2","REA3","REA4"
)
colDataEA <- data.frame(sampleIDEA, type_LR)
ddsEA <- DESeqDataSetFromTximport(txiEA, colDataEA, design= ~ type_LR)
ddsEA <- DESeq(ddsEA)
Ba2025EA_DESeq2 <- results(ddsEA)
write.table(Ba2025EA_DESeq2, file = "Ba2025EA_DESeq2.txt")


filesMS <- c(
  "quants/L_MS1_quant/quant.sf","quants/L_MS2_quant/quant.sf","quants/L_MS3_quant/quant.sf","quants/L_MS4_quant/quant.sf",
  "quants/R_MS1_quant/quant.sf","quants/R_MS2_quant/quant.sf","quants/R_MS3_quant/quant.sf","quants/R_MS4_quant/quant.sf"
)
txiMS <- tximport(filesMS, type = "salmon", tx2gene = Ba2025.ST.tx2gene)
sampleIDMS = c(
  "LMS1","LMS2","LMS3","LMS4",
  "RMS1","RMS2","RMS3","RMS4"
)
colDataMS <- data.frame(sampleIDMS, type_LR)
ddsMS <- DESeqDataSetFromTximport(txiMS, colDataMS, design= ~ type_LR)
ddsMS <- DESeq(ddsMS)
Ba2025MS_DESeq2 <- results(ddsMS)
write.table(Ba2025MS_DESeq2, file = "Ba2025MS_DESeq2.txt")



filesMA <- c(
  "quants/L_MA1_quant/quant.sf","quants/L_MA2_quant/quant.sf","quants/L_MA3_quant/quant.sf","quants/L_MA4_quant/quant.sf",
  "quants/R_MA1_quant/quant.sf","quants/R_MA2_quant/quant.sf","quants/R_MA3_quant/quant.sf","quants/R_MA4_quant/quant.sf"
)
txiMA <- tximport(filesMA, type = "salmon", tx2gene = Ba2025.ST.tx2gene)
sampleIDMA = c(
  "LMA1","LMA2","LMA3","LMA4",
  "RMA1","RMA2","RMA3","RMA4"
)
colDataMA <- data.frame(sampleIDMA, type_LR)
ddsMA <- DESeqDataSetFromTximport(txiMA, colDataMA, design= ~ type_LR)
ddsMA <- DESeq(ddsMA)
Ba2025MA_DESeq2 <- results(ddsMA)
write.table(Ba2025MA_DESeq2, file = "Ba2025MA_DESeq2.txt")


filesLS <- c(
  "quants/L_LS1_quant/quant.sf","quants/L_LS2_quant/quant.sf","quants/L_LS3_quant/quant.sf","quants/L_LS4_quant/quant.sf",
  "quants/R_LS1_quant/quant.sf","quants/R_LS2_quant/quant.sf","quants/R_LS3_quant/quant.sf","quants/R_LS4_quant/quant.sf"
)
txiLS <- tximport(filesLS, type = "salmon", tx2gene = Ba2025.ST.tx2gene)
sampleIDLS = c(
  "LLS1","LLS2","LLS3","LLS4",
  "RLS1","RLS2","RLS3","RLS4"
)
colDataLS <- data.frame(sampleIDLS, type_LR)
ddsLS <- DESeqDataSetFromTximport(txiLS, colDataLS, design= ~ type_LR)
ddsLS <- DESeq(ddsLS)
Ba2025LS_DESeq2 <- results(ddsLS)
write.table(Ba2025LS_DESeq2, file = "Ba2025LS_DESeq2.txt")

filesLA <- c(
  "quants/L_LA1_quant/quant.sf","quants/L_LA2_quant/quant.sf","quants/L_LA3_quant/quant.sf","quants/L_LA4_quant/quant.sf",
  "quants/R_LA1_quant/quant.sf","quants/R_LA2_quant/quant.sf","quants/R_LA3_quant/quant.sf"
)
txiLA <- tximport(filesLA, type = "salmon", tx2gene = Ba2025.ST.tx2gene)
sampleIDLA = c(
  "LLA1","LLA2","LLA3","LLA4",
  "RLA1","RLA2","RLA3"
)
type_LR2 = c("L","L","L","L","R","R","R")

colDataLA <- data.frame(sampleIDLA, type_LR2)
ddsLA <- DESeqDataSetFromTximport(txiLA, colDataLA, design= ~ type_LR2)
ddsLA <- DESeq(ddsLA)
Ba2025LA_DESeq2 <- results(ddsLA)
write.table(Ba2025LA_DESeq2, file = "Ba2025LA_DESeq2.txt")

# Subset to -5 to 5 and remove NA
filtered_df <- Ba2025ES_DESeq2 %>%
  filter(!is.na(baseMean), !is.na(log2FoldChange),
         baseMean > 100, log2FoldChange < -2)


 x <- x[!is.na(x) & x <= -2]

hist(x,
     breaks = seq(-5, 5, by = 0.2),
     col = "skyblue",
     border = "white",
     main = "Histogram of log2FoldChange - Ba2025ES",
     xlab = "log2FoldChange",
     ylab = "Frequency")

