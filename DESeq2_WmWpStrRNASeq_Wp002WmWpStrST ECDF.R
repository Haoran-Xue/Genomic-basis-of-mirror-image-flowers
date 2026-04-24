library(ggplot2)
library(dplyr)
library(tidyr)



# Colors for bins
color_palette_bins <- c(
  "IAA" = "#00abf0",
  "ERF" = "#eb287b",
  "other bins" = "gray"
)

# Line types for groups
linetype_palette_group <- c("Wp" = "solid", "Wm" = "dashed")

# Store median differences
all_diffs <- list()

# oragn type
o_vec <- c("ES", "EA", "MS", "MA", "LS", "LA")

for (o in o_vec) {
  WpStr_df_name <- paste0("WpStr", o, "_Mercator")
  
  if (!exists(WpStr_df_name) | !exists(Wm_df_name)) {
    message("Skipping ", o, ": missing one or both groups")
    next
  }
  
  WpStr_df <- get(WpStr_df_name)
  Wm_df <- get(Wm_df_name)
  
  # Plotting
  pdf(file = paste0("Wp_Wm_", o, "_ECDF_IAA+ERF_no35.2.pdf"), width = 16, height = 8)
  p <- ggplot() +
    labs(y = "F(log2FC)", x = "log2FC") +
    xlim(-5, 5) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  get_bin_filter <- function(df, bin_pattern) {
    df[grep(bin_pattern, df$bin), ]
  }
  
  exclude_pattern <- "^15\\.5\\.7\\.6|,15\\.5\\.7\\.6|^11\\.2\\.2\\.2|,11\\.2\\.2\\.2|^35\\.2|,35\\.2"
  
  # Wp
  p <- p + stat_ecdf(data = get_bin_filter(WpStr_df, "^15\\.5\\.7\\.6|,15\\.5\\.7\\.6"),
                     aes(log2FC, color = "ERF", linetype = "Wp"), size = 1)
  p <- p + stat_ecdf(data = get_bin_filter(WpStr_df, "^11\\.2\\.2\\.2|,11\\.2\\.2\\.2"),
                     aes(log2FC, color = "IAA", linetype = "Wp"), size = 1)
  p <- p + stat_ecdf(data = WpStr_df[!grepl(exclude_pattern, WpStr_df$bin), ],
                     aes(log2FC, color = "other bins", linetype = "Wp"), size = 1)
  
  # Wm
  p <- p + stat_ecdf(data = get_bin_filter(Wm_df, "^15\\.5\\.7\\.6|,15\\.5\\.7\\.6"),
                     aes(log2FC, color = "ERF", linetype = "Wm"), size = 1)
  p <- p + stat_ecdf(data = get_bin_filter(Wm_df, "^11\\.2\\.2\\.2|,11\\.2\\.2\\.2"),
                     aes(log2FC, color = "IAA", linetype = "Wm"), size = 1)
  p <- p + stat_ecdf(data = Wm_df[!grepl(exclude_pattern, Wm_df$bin), ],
                     aes(log2FC, color = "other bins", linetype = "Wm"), size = 1)
  
  group_labels <- c(
    Wp = expression(italic("W. paniculata")),
    Wm = expression(italic("W. multiflora"))
  )
  
  p <- p + scale_color_manual(
    name = "Bin",
    values = color_palette_bins,
    breaks = c("ERF", "IAA", "other bins")
  ) +
    scale_linetype_manual(
      name = "Species",
      values = linetype_palette_group,
      breaks = c("Wp", "Wm"),
      labels = group_labels
    ) +
    guides(
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2, keywidth = 3)
    ) +
    ggtitle(paste0("ECDF of log2FC for ", o, " (Wp vs Wm)"))
  
  print(p)
  dev.off()
  message("Saved plot for ", o)
  
  # --- Median Difference Calculation ---
  for (pair in list(list(df = WpStr_df, species = "Wp"), list(df = Wm_df, species = "Wm"))) {
    df <- pair$df
    species <- pair$species
    
    df <- df %>%
      mutate(
        group = case_when(
          grepl("^15\\.5\\.7\\.6|,15\\.5\\.7\\.6", bin) ~ "ERF",
          grepl("^11\\.2\\.2\\.2|,11\\.2\\.2\\.2", bin) ~ "IAA",
          !grepl(exclude_pattern, bin) ~ "other",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(group)) %>%
      group_by(group) %>%
      summarise(median_log2FC = median(log2FC, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = median_log2FC) %>%
      mutate(
        stage_organ = o,
        species = species,
        diff_ERF = if ("ERF" %in% names(.) && "other" %in% names(.)) ERF - other else NA,
        diff_IAA = if ("IAA" %in% names(.) && "other" %in% names(.)) IAA - other else NA
      ) %>%
      select(stage_organ, species, diff_ERF, diff_IAA)
    
    all_diffs[[length(all_diffs) + 1]] <- df
  }
}

# Final table with all median differences
median_diff_df <- bind_rows(all_diffs)
print(median_diff_df)





library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)

setwd("Y:/hxue/W.paniculata_W.multiflora_RNA-Seq/DESeq2 results with Mercator4 withSAURs")

WpStrMA_Mercator <- read.table("WpStrMA_Wp002WmWpStrST_DESeq2.noNA.Mercator4_SAURs.bin.merged.txt",
                          sep = "\t", header=TRUE, comment.char="", quote="", fill=T, encoding="UTF-8")
WmMA_Mercator <- read.table("WmMA_DESeq2.noNA.Mercator.withSAURs.bin.merged.txt",
                            sep = "\t", header=TRUE, comment.char="", quote="", fill=T, encoding="UTF-8")


color_palette <- c(viridis(4, option = "H"))

WmMA_Mercator$baseMean <- as.numeric(WmMA_Mercator$baseMean)
WmMA_Mercator$log2FC <- as.numeric(WmMA_Mercator$log2FC)

pdf(file="log2FC_vs_log10baseMean.pdf", width=16, height=8)

p<-ggplot(WmMA_Mercator, aes(x=log10(baseMean), y=log2FC)) +
  geom_point(size=0.5, stroke = 0.5, color = "grey") +
  geom_point(data = WmMA_Mercator[grep("^15\\.5\\.7\\.1\\.|,15\\.5\\.7\\.1\\.", WmMA_Mercator$bin), ], size=1, stroke = 0.5, aes(colour = "NAC transcription factor")) +
  geom_point(data = WmMA_Mercator[grep("^15\\.5\\.7\\.5\\.1\\.|,15\\.5\\.7\\.5\\.1\\.", WmMA_Mercator$bin), ], size=1, stroke = 0.5, aes(colour = "WRKY transcription factor")) +
  geom_point(data = WmMA_Mercator[grep("^15\\.5\\.7\\.6\\.11\\.|,15\\.5\\.7\\.6\\.11\\.", WmMA_Mercator$bin), ], size=1, stroke = 0.5, aes(colour = "ERF-IX transcription factor")) +
  geom_point(data = WmMA_Mercator[grep("^11\\.2\\.2\\.2|,11\\.2\\.2\\.2", WmMA_Mercator$bin), ], size=1, stroke = 0.5, aes(colour = "Aux/IAA transcriptional repressor")) +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=15)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=15)) +
  theme(axis.title.y = element_text(angle = 00, vjust = 0.5, size=15)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=15)) +
  geom_line(aes(y=0), color = "black") +
  scale_color_manual(name='Mercator bins',
                       breaks = c("NAC transcription factor",
                                  "WRKY transcription factor",
                                  "ERF-IX transcription factor",
                                  "Aux/IAA transcriptional repressor"),
                       values = c("NAC transcription factor" = color_palette[1],
                                  "WRKY transcription factor" = color_palette[2],
                                  "ERF-IX transcription factor" = color_palette[3],
                                  "Aux/IAA transcriptional repressor" = color_palette[4])) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) 
p
dev.off()

setEPS()
postscript("log2FC_vs_log10baseMeanr.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
p
dev.off()

pdf(file="ECDF.pdf", width=16, height=8)  
ggplot(WmMA_Mercator, aes(log2FC)) +
  stat_ecdf(data = WmMA_Mercator[grep("^15\\.5\\.7\\.1|,15\\.5\\.7\\.1", WmMA_Mercator$bin), ], size=0.5, aes(color = "NAC transcription factor")) +
  labs(y = "F(log2FC)", x="log2FC") +
  xlim(-5,5)+
  stat_ecdf(data = WmMA_Mercator[grep("^15\\.5\\.7\\.5\\.1|,15\\.5\\.7\\.5\\.1", WmMA_Mercator$bin), ], aes(color = "WRKY transcription factor")) +
  stat_ecdf(data = WmMA_Mercator[grep("^15\\.5\\.7\\.6\\.1\\.11|,15\\.5\\.7\\.6\\.1\\.11", WmMA_Mercator$bin), ], aes(color = "ERF-IX transcription factor")) +
  stat_ecdf(data = WmMA_Mercator[grep("^11\\.2\\.2\\.2|,11\\.2\\.2\\.2", WmMA_Mercator$bin), ], aes(color = "Aux/IAA transcriptional repressor")) +
  stat_ecdf(data = WmMA_Mercator[!grepl("^21\\.3\\.1\\.2\\.|,21\\.3\\.1\\.2\\.|^21\\.4\\.1\\.1\\.2\\.|,21\\.4\\.1\\.1\\.2\\.|
                                     ^11\\.11\\.2\\.4\\.1|,11\\.11\\.2\\.4\\.1|11\\.2\\.2\\.6|,11\\.2\\.2\\.6|
                                     11\\.2\\.2\\.2|,11\\.2\\.2\\.2|20\\.1\\.1\\.|,20\\.1\\.1\\.", WmMA_Mercator$bin), ],
            aes(color = "other bins"), geom = "step") +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=15)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=15)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=15)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=15)) +
  scale_color_manual(name='Mercator bins',
                     breaks = c("NAC transcription factor",
                                "WRKY transcription factor",
                                "ERF-IX transcription factor",
                                "Aux/IAA transcriptional repressor",
                                "other bins"),
                     values = c("NAC transcription factor" = color_palette[1],
                                "WRKY transcription factor" = color_palette[2],
                                "ERF-IX transcription factor" = color_palette[3],
                                "Aux/IAA transcriptional repressor" = color_palette[4],
                                "other bins" = "grey")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) 

dev.off()


WpStrMA_Mercator <- read.table("DESeq2 results with Mercator4 withSAURs/WpStrMA_Wp002WmWpStrST_DESeq2.noNA.Mercator4_SAURs.bin.merged.txt",
                               sep = "\t", header=TRUE, comment.char="", quote="", fill=T, encoding="UTF-8")

pdf(file="WpStrMA_ECDF_2.pdf", width=10, height=6)  
ggplot(WpStrMA_Mercator, aes(log2FC)) +
  stat_ecdf(data = WpStrMA_Mercator[grep("^11\\.2\\.2\\.2|,11\\.2\\.2\\.2", WpStrMA_Mercator$bin), ], size=0.5, aes(color = "Aux/IAA")) +
  labs(y = "ECDF(log2FC)", x="log2FC") +
  xlim(-5,5)+
  stat_ecdf(data = WpStrMA_Mercator[!grepl("^11\\.2\\.2\\.2|,11\\.2\\.2\\.2", WpStrMA_Mercator$bin), ],
            aes(color = "other bins"), geom = "step") +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=20)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=20)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=20)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=20)) +
  scale_color_manual(name='Functional bins',
                     breaks = c("Aux/IAA",
                                "other bins"),
                     values = c("Aux/IAA" = "red",
                                "other bins" = "grey")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.title=element_text(size=22)) +
  theme(legend.text=element_text(size=22)) 

dev.off()

pdf(file="WpStrMA_ECDF_35.pdf", width=10, height=6)  
ggplot(WpStrMA_Mercator, aes(log2FC)) +
  stat_ecdf(data = WpStrMA_Mercator[grep("^35|,35", WpStrMA_Mercator$bin), ], size=0.5, aes(color = "No Mercator4 annotation")) +
  labs(y = "ECDF(log2FC)", x="log2FC") +
  xlim(-5,5)+
  stat_ecdf(data = WpStrMA_Mercator[!grepl("^35|,35", WpStrMA_Mercator$bin), ],
            aes(color = "other bins"), geom = "step") +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=20)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=20)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=20)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=20)) +
  scale_color_manual(name='Functional bins',
                     breaks = c("No Mercator4 annotation",
                                "other bins"),
                     values = c("No Mercator4 annotation" = "red",
                                "other bins" = "grey")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.title=element_text(size=22)) +
  theme(legend.text=element_text(size=22)) 

dev.off()


WpStrMA_Mercator <- read.table("DESeq2 results with Mercator4 withSAURs/WpStrMA_Wp002WmWpStrST_DESeq2.noNA.Mercator4_SAURs.bin.merged.txt",
                               sep = "\t", header=TRUE, comment.char="", quote="", fill=T, encoding="UTF-8")

pdf(file="WpStrMA_ECDF_2.pdf", width=10, height=6)  
ggplot(WpStrMA_Mercator, aes(log2FC)) +
  stat_ecdf(data = WpStrMA_Mercator[grep("^11\\.2\\.2\\.2|,11\\.2\\.2\\.2", WpStrMA_Mercator$bin), ], size=0.5, aes(color = "Aux/IAA")) +
  labs(y = "ECDF(log2FC)", x="log2FC") +
  xlim(-5,5)+
  stat_ecdf(data = WpStrMA_Mercator[!grepl("^11\\.2\\.2\\.2|,11\\.2\\.2\\.2", WpStrMA_Mercator$bin), ],
            aes(color = "other bins"), geom = "step") +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 00, hjust = 1, size=20)) +
  theme(axis.text.x = element_text(angle = 00, hjust = 1, size=20)) +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, size=20)) +
  theme(axis.title.x = element_text(angle = 00, hjust = 0.5, size=20)) +
  scale_color_manual(name='Functional bins',
                     breaks = c("Aux/IAA",
                                "other bins"),
                     values = c("Aux/IAA" = "red",
                                "other bins" = "grey")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.title=element_text(size=22)) +
  theme(legend.text=element_text(size=22)) 

dev.off()
