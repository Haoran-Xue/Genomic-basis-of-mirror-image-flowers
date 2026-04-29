install.packages("ggplot2")
install.packages("ggpubr")

library(ggpubr)
library(ggplot2)
library(tidyr)
library(dplyr)



setwd("D:/Michael Lenhard - HFSP Shared folder/Universität Potsdam/Michael Lenhard - HFSP Shared folder/Wachendorfia ms1/Figure 5 R supergene")
#setwd("G:/Universität Potsdam/Michael Lenhard - HFSP Shared folder/Wachendorfia ms1/Figure 5 R supergene")

#repeat_content <- data.frame(read.csv(file = "repeat_content_comparison.csv", header = TRUE, row.names = NULL))

df <- data.frame(
  Wp.E = c(0.8666, 0.9159, 0.4403, 0.8885, 0.7751, 0.9672, 0.8358, 0.9992, 0.9698, 0.8971,
           0.5839, 0.5330, 0.8635, 0.7875, 0.6582, 0.9128, 0.8087, 0.9496, 0.7872, 0.6781,NA,NA,NA,NA,NA),
  Wp.N = c(0.7707, 0.8790, 0.9765, 0.7663, 0.9871, 0.7363, 0.8776, 0.9251, 0.7683, 0.9223,
           0.6790, 0.8557, 0.9298, 0.6234, 0.6511, 0.6930, 0.6557, 0.7664, 0.8232, 0.9583,NA,NA,NA,NA,NA),
  Wt.E = c(0.9801, 0.4553, 0.9532, 0.9548, 0.8844, 0.9310, 0.9512, 0.8755, 0.9984, 0.9065,
           0.7057, 0.8157, 0.8025, 0.8918, 0.8693, 0.9458, 0.8202, 0.9089, 0.7971, 0.9181,NA,NA,NA,NA,NA),
  Wt.N = c(0.6415, 0.9134, 0.9584, 0.9820, 0.9498, 0.9203, 0.9696, 0.9941, 0.9812, 0.9875,
           0.9999, 0.9910, 0.9989, 0.9744, 0.9286, 0.9663, 0.9251, 0.9677, 0.7688, 0.9999,NA,NA,NA,NA,NA),
#  Ba.E = c(0.9963, 0.9741, 0.5956, 0.9300, 0.9560, 0.8014, 0.8042, 0.6819, 0.8618, 1.0000,
#           NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#  Ba.N = c(0.8515, 0.9826, 0.6606, 0.8003, 0.9280, 0.6690, 0.8098, 0.8627, 0.9982, 0.7821,
#           NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  BaN.E = c(0.6833, 0.7214, 0.9946, 0.8498, 0.9358, 0.9961, 0.9965, 0.9856, 0.9732, 0.9911,
           0.5494, 0.7934, 0.9993, 0.957, 0.8937, 0.8266, 0.8361, 0.9017, 0.6922, 0.7295, 0.5899, 0.85, 0.7644, 0.8602, 0.9944),
  BaN.N = c(0.454, 0.6098, 0.997, 0.6935, 0.4042, 0.9978, 0.8489, 0.7139, 0.755, 0.8106, 
           0.8622, 0.831, 0.7906, 0.7587, 0.8039, 0.2898, 0.5767, 0.9993, 0.9992, 0.9991, 0.921, 0.8511, 0.3369, 0.0781, 0.7025)
  
)


df_long <- df %>%
  mutate(ID = row_number()) %>%  # keep row id if needed
  pivot_longer(cols = -ID, names_to = "Species_Condition", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  separate(Species_Condition, into = c("Species", "Condition"), sep = "\\.") %>%
  mutate(
    Species = factor(Species, levels = c("Wp", "Wt", "BaN")),
    Condition = factor(Condition, levels = c("E", "N"), labels = c("R-locus", "Neighbouring"))
  )

# Species labels for italics
species_labels <- c(
  "Wp" = "italic('W. paniculata')",
  "Wt" = "italic('W. thyrsiflora')",
  "BaN" = "italic('B. aurea') * ' (Ngeli)'"
)

pvals <- df_long %>%
  group_by(Species) %>%
  summarise(p = wilcox.test(Value ~ Condition)$p.value) %>%
  mutate(
    label = paste0("italic(p)*' = '*", signif(p, 3)),
    Species = factor(Species, levels = c("Wp", "Wt", "BaN"))
  )

# Plot
p <- ggplot(df_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~Species, labeller = as_labeller(species_labels, label_parsed)) +
  scale_fill_manual(values = c("R-locus" = "#4daf4a", "Neighbouring" = "#e41a1c")) +
  labs(
    y = "Proportion of repeats",
    x = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  geom_text(
    data = pvals,
    aes(x = 1.5, y = max(df_long$Value) * 1.05, label = label), 
    parse = TRUE,
    inherit.aes = FALSE,
    size = 4  
  )


print(p)

# Save to PDF
pdf(file = "R-locus_vs_neighbouring_repeat_content.pdf", width = 5, height = 3)
print(p)
dev.off()

