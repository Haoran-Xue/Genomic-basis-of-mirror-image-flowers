install.packages("ggplot2")
install.packages("cowplot")

library(ggplot2)
library(cowplot)


setwd("D:/Michael Lenhard - HFSP Shared folder/Universität Potsdam/Michael Lenhard - HFSP Shared folder/Wachendorfia/Wachendorfia_DNA-Seq_essentials")

#50 kb
W.RM.ptg000001l.50kb.norm.LtoR <- data.frame(read.table(file = "W_Wp002.RM.ptg000001l.50kb.regions.LtoR.norm.txt", header = TRUE, row.names = NULL))

W.4 <- subset(W.RM.ptg000001l.50kb.norm.LtoR, W.RM.ptg000001l.50kb.norm.LtoR$from >= 22750000 & W.RM.ptg000001l.50kb.norm.LtoR$from < 23950000)
write.csv(W.4, file = "Wachendorfia_pools_RM_LtoR.50kb.ptg000001l_22750000-23950000.csv", row.names = FALSE)

# Add position column (in MB)
W.4 <- W.4 %>%
  mutate(pos = (from + 25000) / 1e6)

# Extract flanking y-values for dashed segments
get_y_at_from <- function(data, from_val, col) {
  data %>%
    filter(from == from_val) %>%
    pull({{ col }})
}

x_start <- (22900000 + 25000)/1e6
x_end   <- (23050000 + 25000)/1e6

y_wb_start <- get_y_at_from(W.4, 22900000, Wb.LtoR)
y_wb_end   <- get_y_at_from(W.4, 23050000, Wb.LtoR)

y_wm_start <- get_y_at_from(W.4, 22900000, Wm.LtoR)
y_wm_end   <- get_y_at_from(W.4, 23050000, Wm.LtoR)

y_wp_start <- get_y_at_from(W.4, 22900000, Wp.LtoR)
y_wp_end   <- get_y_at_from(W.4, 23050000, Wp.LtoR)

y_wt_start <- get_y_at_from(W.4, 22900000, Wt.LtoR)
y_wt_end   <- get_y_at_from(W.4, 23050000, Wt.LtoR)



pdf(file="Wachendorfia_pools_RM_LtoR.50kb.ptg000001l_22750000-23950000.pdf", width=25, height=15)

pw.4 <- ggplot(W.4, aes(x=(from+25000)/1000000)) +
  ylim(0,1.75) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wb.LtoR, color="W. brachyandra"),size = 1) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wm.LtoR, color ="W. multiflora"),size = 1) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wp.LtoR, color="W. paniculata"), size = 1) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wt.LtoR, color="W. thyrsiflora"), size = 1) +
  theme_bw() +
  scale_color_manual(
    name = NULL,
    values = c("W. brachyandra" = "black", 
               "W. multiflora" = "red", 
               "W. paniculata" = "blue", 
               "W. thyrsiflora" = "turquoise"),
    labels = c(
      expression(italic("W. brachyandra")),
      expression(italic("W. multiflora")),
      expression(italic("W. paniculata")),
      expression(italic("W. thyrsiflora"))
    )
  ) +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_blank()
  ) +
  annotate("rect", xmin = 23273529/1e6, xmax = 23276494/1e6,
           ymin = 0, ymax = 1.75,
           fill = "#00abf0", alpha = 0.3) +
  annotate("rect", xmin = 23355346/1e6, xmax = 23357266/1e6,
           ymin = 0, ymax = 1.75,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (23273529 + 23276494)/2 / 1e6, y = 1.7, 
           label = expression(italic("YUC-R")), size = 7, fontface = "bold", color = "black") +
  annotate("text", x = (23355346 + 23357266)/2 / 1e6, y = 1.7, 
           label = expression(italic("miR156-R")), size = 7, fontface = "bold", color = "black")+
  geom_segment(aes(x = x_start, xend = x_end, y = y_wb_start, yend = y_wb_end),
               linetype = "dashed", color = "black", size = 1) +
  geom_segment(aes(x = x_start, xend = x_end, y = y_wm_start, yend = y_wm_end),
               linetype = "dashed", color = "red", size = 1) +
  geom_segment(aes(x = x_start, xend = x_end, y = y_wp_start, yend = y_wp_end),
               linetype = "dashed", color = "blue", size = 1) +
  geom_segment(aes(x = x_start, xend = x_end, y = y_wt_start, yend = y_wt_end),
               linetype = "dashed", color = "turquoise", size = 1)


plot(pw.4)
dev.off()


pdf(file="Wp_pools_RM_LandR.50kb.ptg000001l_22750000-23950000.pdf", width=25, height=15)

pw.5 <- ggplot(W.4, aes(x=(from+25000)/1000000)) +
  ylim(0,0.5) +
  ylab("Normalized coverage") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
#  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wb.L, color="W. brachyandra"),size = 1,linetype = "dashed") +
#  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wm.L, color ="W. multiflora"),size = 1,linetype = "dashed") +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = (WpAB.L+WpC.L)/2), color="black", size = 1,linetype = "dashed") +
#  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wt.L, color="W. thyrsiflora"), size = 1,linetype = "dashed") +
#  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wb.R, color="W. brachyandra"),size = 1) +
#  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wm.R, color ="W. multiflora"),size = 1) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = (WpAB.R+WpC.R)/2), color="black", size = 1) +
#  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wt.R, color="W. thyrsiflora"), size = 1) +
  theme_bw() +
#  scale_color_manual(
#    name = NULL,
#    values = c(
#      "W. brachyandra" = "black", 
#               "W. multiflora" = "red", 
#               "W. paniculata" = "black"),
#               "W. thyrsiflora" = "turquoise"),
#    labels = c(
#      expression(italic("W. brachyandra")),
#      expression(italic("W. multiflora")),
#      expression(italic("W. paniculata"))
#      expression(italic("W. thyrsiflora"))
#    )
#  ) +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  annotate("rect", xmin = 23273529/1e6, xmax = 23276494/1e6,
           ymin = 0, ymax = 0.5,
           fill = "#00abf0", alpha = 0.3) +
  annotate("rect", xmin = 23355346/1e6, xmax = 23357266/1e6,
           ymin = 0, ymax = 0.5,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (23273529 + 23276494)/2 / 1e6, y = 0.48, 
           label = expression(italic("YUC-R")), size = 7, fontface = "bold", color = "black") +
  annotate("text", x = (23355346 + 23357266)/2 / 1e6, y = 0.46, 
           label = expression(italic("miR156-R")), size = 7, fontface = "bold", color = "black")

plot(pw.5)
dev.off()



Wt_Wt302R.RM.ptg000020l.50kb.norm.LtoR <- data.frame(read.table(file = "Wt_Wt302R.RM.ptg000020l.50kb.LtoR.norm.txt", header = TRUE, row.names = NULL))


Wt.1 <- subset(Wt_Wt302R.RM.ptg000020l.50kb.norm.LtoR, Wt_Wt302R.RM.ptg000020l.50kb.norm.LtoR$from >= 10850000 & Wt_Wt302R.RM.ptg000020l.50kb.norm.LtoR$from < 12050000)
write.csv(Wt.1, file = "Wt_pools_RM_LtoR.50kb.h2tg000020l_10850000-12050000.csv", row.names = FALSE)


wt_x_start <- (11250000 + 25000)/1e6
wt_x_end   <- (11400000 + 25000)/1e6

wt_y_start <- get_y_at_from(Wt.1, 11250000, Wt.LtoR)
wt_y_end   <- get_y_at_from(Wt.1, 11400000, Wt.LtoR)

pdf(file="Wt_pools_RM_LtoR.50kb.ptg000020l_10850000-12050000.pdf", width=25, height=15)

pWt.1 <- ggplot(Wt.1, aes(x=(from+25000)/1000000)) +
  ylim(0,1.55) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = Wt.1, aes(x=(from+25000)/1000000, y = Wt.LtoR), color="black",size = 1) +
  theme_bw() +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")) +
  annotate("rect", xmin = 11413091/1e6, xmax = 11416222/1e6,
           ymin = 0, ymax = 1.55,
           fill = "#00abf0", alpha = 0.3) +
  annotate("rect", xmin = 11502790/1e6, xmax = 11504661/1e6,
           ymin = 0, ymax = 1.55,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (11413091 + 11416222)/2 / 1e6, y = 1.55, 
           label = expression(italic("YUC-R")), size = 7, fontface = "bold", color = "black") +
  annotate("text", x = (11502790 + 11504661)/2 / 1e6, y = 1.55, 
           label = expression(italic("miR156-R")), size = 7, fontface = "bold", color = "black")+
  geom_segment(aes(x = wt_x_start, xend = wt_x_end, y = wt_y_start, yend = wt_y_end),
               linetype = "dashed", color = "black", size = 1)

plot(pWt.1)

dev.off()


Wt.2 <- subset(Wt_Wt302R.RM.ptg000020l.50kb.norm.LtoR, Wt_Wt302R.RM.ptg000020l.50kb.norm.LtoR$from >= 10800000 & Wt_Wt302R.RM.ptg000020l.50kb.norm.LtoR$from < 12200000)
write.csv(Wt.2, file = "Wt_pools_RM_LtoR.50kb.ptg000020l_10800000-12200000.csv", row.names = FALSE)

pdf(file="Wt_pools_RM_LtoR.50kb.ptg000020l_10800000-12200000.pdf", width=25, height=15)

pWt.2 <- ggplot(Wt.2, aes(x=(from+25000)/1000000)) +
  ylim(0,2.1) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = Wt.2, aes(x=(from+25000)/1000000, y = Wt.LtoR), color="black",size = 1) +
  theme_bw() +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")) +
  annotate("rect", xmin = 11413091/1e6, xmax = 11416222/1e6,
           ymin = 0, ymax = 2.1,
           fill = "#00abf0", alpha = 0.3) +
  annotate("rect", xmin = 11502790/1e6, xmax = 11504661/1e6,
           ymin = 0, ymax = 2.1,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (11413091 + 11416222)/2 / 1e6, y = 2.05, 
           label = expression(italic("YUC-R")), size = 7, fontface = "bold", color = "black") +
  annotate("text", x = (11502790 + 11504661)/2 / 1e6, y = 2.05, 
           label = expression(italic("miR156-R")), size = 7, fontface = "bold", color = "black")+
  geom_segment(aes(x = wt_x_start, xend = wt_x_end, y = wt_y_start, yend = wt_y_end),
               linetype = "dashed", color = "black", size = 1)

plot(pWt.2)
dev.off()





pdf(file="Wt_pools_RM_LandR.50kb.ptg000020l_10800000-12200000.pdf", width=25, height=15)

pWt.3 <- ggplot(Wt.2, aes(x=(from+25000)/1000000)) +
  ylim(0,0.16) +
  ylab("Normalized coverage") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = Wt.2, aes(x=(from+25000)/1000000, y = Wt.R), color="black",size = 1) +
  geom_line(data = Wt.2, aes(x=(from+25000)/1000000, y = Wt.L), color="black",size = 1,linetype = "dashed") +
  theme_bw() +
  scale_color_manual(
    name = NULL,
    values = c("L" = "black", 
               "R" = "black"),
    labels = c("L", "R")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("L" = "dashed", 
               "R" = "solid"),
    labels = c("L", "R")
  ) +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")) +
  annotate("rect", xmin = 11413091/1e6, xmax = 11416222/1e6,
           ymin = 0, ymax = 0.16,
           fill = "#00abf0", alpha = 0.3) +
  annotate("rect", xmin = 11502790/1e6, xmax = 11504661/1e6,
           ymin = 0, ymax = 0.16,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (11413091 + 11416222)/2 / 1e6, y = 0.155, 
           label = expression(italic("YUC-R")), size = 7, fontface = "bold", color = "black") +
  annotate("text", x = (11502790 + 11504661)/2 / 1e6, y = 0.155, 
           label = expression(italic("miR156-R")), size = 7, fontface = "bold", color = "black")

plot(pWt.3)
dev.off()






B.RM.h2tg000059l.50kb.norm.LtoR <- data.frame(read.table(file = "Ba_Bahap2.RM.h2tg000059l.50kb.regions.LtoR.norm.txt", header = TRUE, row.names = NULL))


B.1 <- subset(B.RM.h2tg000059l.50kb.norm.LtoR, B.RM.h2tg000059l.50kb.norm.LtoR$from >= 600000 & B.RM.h2tg000059l.50kb.norm.LtoR$from < 1800000)
write.csv(B.1, file = "Barberetta_pools_RM_LtoR.50kb.:q!:q_600000-1800000.csv", row.names = FALSE)

Ba_x_start <- (1550000 + 25000)/1e6
Ba_x_end   <- (1650000 + 25000)/1e6

Ba_y_start <- get_y_at_from(B.1, 1550000, Ba.LtoR)
Ba_y_end   <- get_y_at_from(B.1, 1650000, Ba.LtoR)


pdf(file="Barberetta_pools_RM_LtoR.50kb.h2tg000059l_600000-1800000.pdf", width=25, height=15)

pB.1 <- ggplot(B.1, aes(x=(from+25000)/1000000)) +
  ylim(0,1.6) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = B.1, aes(x=(from+25000)/1000000, y = Ba.LtoR), color="black",size = 1) +
  theme_bw() +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")) +
  annotate("rect", xmin = 1173034/1e6, xmax = 1178435/1e6,
           ymin = 0, ymax = 1.6,
           fill = "#00abf0", alpha = 0.3) +
  annotate("rect", xmin = 1226802/1e6, xmax = 1228710/1e6,
           ymin = 0, ymax = 1.6,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (1173034 + 1178435)/2 / 1e6, y = 1.55, 
           label = expression(italic("YUC-R")), size = 7, fontface = "bold", color = "black") +
  annotate("text", x = (1226802 + 1228710)/2 / 1e6, y = 1.50, 
           label = expression(italic("miR156-R")), size = 7, fontface = "bold", color = "black")+
  geom_segment(aes(x = Ba_x_start, xend = Ba_x_end, y = Ba_y_start, yend = Ba_y_end),
               linetype = "dashed", color = "black", size = 1)

plot(pB.1)
dev.off()


pdf(file="Barberetta_pools_RM_LandR.50kb.h2tg000059l_600000-1800000.pdf", width=25, height=15)

pB.2 <- ggplot(B.1, aes(x=(from+25000)/1000000)) +
  ylim(0,0.5) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = B.1, aes(x=(from+25000)/1000000, y = Ba.R), color="black",size = 1) +
  geom_line(data = B.1, aes(x=(from+25000)/1000000, y = Ba.L), color="black",size = 1,linetype = "dashed") +
  theme_bw() +
  scale_color_manual(
    name = NULL,
    values = c("L" = "black", 
               "R" = "black"),
    labels = c("L", "R")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("L" = "dashed", 
               "R" = "solid"),
    labels = c("L", "R")
  ) +
  theme_bw() +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")) +
  annotate("rect", xmin = 1173034/1e6, xmax = 1178435/1e6,
           ymin = 0, ymax = 0.5,
           fill = "#00abf0", alpha = 0.3) +
  annotate("rect", xmin = 1226802/1e6, xmax = 1228710/1e6,
           ymin = 0, ymax = 0.5,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (1173034 + 1178435)/2 / 1e6, y = 0.49, 
           label = expression(italic("YUC-R")), size = 7, fontface = "bold", color = "black") +
  annotate("text", x = (1226802 + 1228710)/2 / 1e6, y = 0.48, 
           label = expression(italic("miR156-R")), size = 7, fontface = "bold", color = "black")+
  geom_segment(aes(x = Ba_x_start, xend = Ba_x_end, y = Ba_y_start, yend = Ba_y_end),
               linetype = "dashed", color = "black", size = 1)

plot(pB.2)
dev.off()






Ba_pools.bwa.BaN.RM.50kb.norm.LtoR.ptg000064l <- data.frame(read.table(file = "Ba_pools.bwa.BaN.RM.50kb.norm.LtoR.ptg000064l.txt", header = TRUE, row.names = NULL))


BaN.1 <- subset(Ba_pools.bwa.BaN.RM.50kb.norm.LtoR.ptg000064l, Ba_pools.bwa.BaN.RM.50kb.norm.LtoR.ptg000064l$from >= 4000000 & Ba_pools.bwa.BaN.RM.50kb.norm.LtoR.ptg000064l$from < 5200000)
write.csv(BaN.1, file = "Ba_pools_BaNgeli_reference_RM.50kb.LtoR.ptg000064l_4000000-5200000.csv", row.names = FALSE)

pdf(file="Ba_pools_BaNgeli_reference_RM.50kb.LtoR.ptg000064l_4000000-5200000.pdf", width=15, height=9)

## LtoR < 1.31

pBaN.1 <- ggplot(BaN.1, aes(x=(from+25000)/1000000)) +
  ylim(0,1.6) +
  ylab("L:R coverage ratio") +
  xlab("Position on contig 64 (Mb)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = BaN.1, aes(x=(from+25000)/1000000, y = LtoR.norm), color="#FFA500",linewidth = 3) +
  theme_bw() +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")) +
  
  annotate("rect", xmin = 4534256/1e6, xmax = 4540425/1e6,
           ymin = 0, ymax = 1.4,
           fill = "#00abf0", alpha = 0.3) +
  annotate("text", x = (4534256 - 50000) / 1e6, y = 1.35, 
           label = "bold(italic('YUC-R'))", size = 7, color = "#00abf0", parse = TRUE) +
  
  annotate("rect", xmin = 4618443/1e6, xmax = 4623320/1e6,
           ymin = 0, ymax = 1.4,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (4623320 + 65000) / 1e6, y = 1.35, 
           label = "bold(italic('MIR156-R'))", size = 7, color = "#eb287b", parse = TRUE) +
  
  annotate("rect", xmin = 4654805 / 1e6, xmax = 4658063 / 1e6,
           ymin = 0, ymax = 1.2,
           fill = "#00FF00", alpha = 0.3) +
  annotate("text", x = (4658063 + 25000) / 1e6, y = 1.15, 
           label = "bold(italic('Pol'))", size = 7,  color = "#00FF00", parse = TRUE)

plot(pBaN.1)
dev.off()


pdf(file="Ba_pools_BaNgeli_reference_RM.50kb.LandR.ptg000064l_4000000-5200000.pdf", width=15, height=9)

pBaN.2 <- ggplot(BaN.1, aes(x=(from+25000)/1000000)) +
  ylim(0,0.45) +
  ylab("L:R coverage ratio") +
  xlab("Position on contig 64 (Mb)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = BaN.1, aes(x=(from+25000)/1000000, y = R.norm), color="black",size = 3) +
  geom_line(data = BaN.1, aes(x=(from+25000)/1000000, y = L.norm), color="black",size = 3, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(
    name = NULL,
    values = c("L" = "black", 
               "R" = "black"),
    labels = c("L", "R")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("L" = "dashed", 
               "R" = "solid"),
    labels = c("L", "R")
  ) +
  theme_bw() +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) +
  theme(
    legend.position = c(0.05, 0.95),          # Upper left inside the plot
    legend.justification = c(0, 1),           # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black")) +
  
  annotate("rect", xmin = 4534256/1e6, xmax = 4540425/1e6,
           ymin = 0, ymax = 0.45,
           fill = "#00abf0", alpha = 0.3) +
  annotate("text", x = (4534256 - 50000) / 1e6, y = 0.44, 
           label = "bold(italic('YUC-R'))", size = 7, color = "#00abf0", parse = TRUE) +
  
  annotate("rect", xmin = 4618443/1e6, xmax = 4623320/1e6,
           ymin = 0, ymax = 0.45,
           fill = "#eb287b", alpha = 0.3) +
  annotate("text", x = (4623320 + 65000) / 1e6, y = 0.44, 
           label = "bold(italic('MIR156-R'))", size = 7, color = "#eb287b", parse = TRUE) +
  
  annotate("rect", xmin = 4654805 / 1e6, xmax = 4658063 / 1e6,
           ymin = 0, ymax = 0.4,
           fill = "#00FF00", alpha = 0.3) +
  annotate("text", x = (4658063 + 25000) / 1e6, y = 0.39, 
           label = "bold(italic('Pol'))", size = 7,  color = "#00FF00", parse = TRUE)

plot(pBaN.2)
dev.off()
