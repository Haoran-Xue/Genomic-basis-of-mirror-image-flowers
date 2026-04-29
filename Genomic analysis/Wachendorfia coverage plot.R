install.packages("ggplot2")
install.packages("cowplot")

library(ggplot2)
library(cowplot)


setwd("//141.89.108.211/genetics-shared/hxue/Wachendorfia paniculata genome assembly")

Wachendorfia.10kb.norm <- data.frame(read.table(file = "Wachendorfia.10kb.norm", header = TRUE, row.names = NULL))
Wachendorfia.100kb.norm <- data.frame(read.table(file = "Wachendorfia.100kb.norm", header = TRUE, row.names = NULL))


datasubset1 <- subset(Wachendorfia.10kb.norm, Wachendorfia.10kb.norm$contig=="ptg000001l" & Wachendorfia.10kb.norm$von >= 0 & Wachendorfia.10kb.norm$von <= 30630000)
png(file="Wachendorfia.10kb.ptg000001l.png", width=1000, height=1000)
p1 <- ggplot() +
  theme(legend.position = "bottom") +
  ylim(0,2) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = WT_left), linetype = "dashed", color="blue",size = 1) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = WT_right), color="blue",size = 1) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = (WpAB_L+WpC_L)/2), linetype = "dashed", color = "black",size = 1) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = (WpAB_R+WpC_R)/2), color = "black",size = 1) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = Wm_L), linetype = "dashed", color="red",size = 1) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = Wm_R), color="red",size = 1) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = WB_left), linetype = "dashed", color="turquoise",size = 1) +
  geom_line(data = datasubset1, aes(x=(from+500)/1000000, y = WB_right), color="turquoise",size = 1) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  ggtitle("ptg000001l") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("W. paniculata" = "black", "W. multiflora"="red", "W. thyrsiflora" = "blue", "W. brachyandra" = "turquoise"), name = "Species") +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 25)) 


plot(p1)

p1
dev.off()



datasubset4 <- subset(Wachendorfia.10kb.norm, Wachendorfia.10kb.norm$contig=="ptg000001l" & Wachendorfia.10kb.norm$von >= 23200000 & Wachendorfia.10kb.norm$von <= 23500000)

png(file="Wachendorfia.10kb.ptg000001l_23200000-23500000.png", width=1000, height=1000)


p4 <- ggplot(datasubset4, aes(x=(von+5000)/1000000)) +
  ylim(0,5) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = WT_left), linetype = "dashed", color="blue",size = 1) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = WT_right), color="blue",size = 1) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = (WpAB_L+WpC_L)/2), linetype = "dashed", color = "black",size = 1) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = (WpAB_R+WpC_R)/2), color = "black",size = 1) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = Wm_L), linetype = "dashed", color="red",size = 1) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = Wm_R), color="red",size = 1) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = WB_left), linetype = "dashed", color="turquoise",size = 1) +
  geom_line(data = datasubset4, aes(x=(von+5000)/1000000, y = WB_right), color="turquoise",size = 1) +
  ggtitle("ptg000001l") +
  theme_bw() +
  scale_color_manual(values = c("W. paniculata" = "black", "W. multiflora"="red", "W. thyrsiflora" = "blue", "W. brachyandra" = "turquoise"), name = "Species") +
  theme(text = element_text(size = 25)) 
plot(p4)

dev.off()



datasubset2.3 <- subset(Wachendorfia.100kb.norm, Wachendorfia.100kb.norm$contig=="ptg000001l" & Wachendorfia.100kb.norm$von >= 22800000 & Wachendorfia.100kb.norm$von <= 23800000)

png(file="Wachendorfia.100kb.ptg000001l_22800000-23810000.png", width=1000, height=1000)


p2.3 <- ggplot(datasubset2.3, aes(x=(from+500)/1000000)) +
  ylim(0,5) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = WT_left), linetype = "dashed", color="blue",size = 1) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = WT_right), color="blue",size = 1) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = (WpAB_L+WpC_L)/2), linetype = "dashed", color = "black",size = 1) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = (WpAB_R+WpC_R)/2), color = "black",size = 1) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = Wm_L), linetype = "dashed", color="red",size = 1) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = Wm_R), color="red",size = 1) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = WB_left), linetype = "dashed", color="turquoise",size = 1) +
  geom_line(data = datasubset2.3, aes(x=(from+500)/1000000, y = WB_right), color="turquoise",size = 1) +
  ggtitle("ptg000001l") +
  theme_bw() +
  scale_color_manual(values = c("W. paniculata" = "black", "W. multiflora"="red", "W. thyrsiflora" = "blue", "W. brachyandra" = "turquoise"), name = "Species") +
  theme(text = element_text(size = 25)) 
plot(p2.3)

dev.off()



datasubset2.4 <- subset(Wachendorfia.100kb.norm, Wachendorfia.100kb.norm$contig=="ptg000001l" & Wachendorfia.100kb.norm$von >= 23200000 & Wachendorfia.100kb.norm$von <= 23500000)

png(file="Wachendorfia.100kb.ptg000001l_23200000-23500000.png", width=1000, height=1000)

p2.4 <- ggplot(datasubset2.4, aes(x=(from+500)/1000000)) +
  ylim(0,5) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = WT_left), linetype = "dashed", color="blue",size = 1) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = WT_right), color="blue",size = 1) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = (WpAB_L+WpC_L)/2), linetype = "dashed", color = "black",size = 1) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = (WpAB_R+WpC_R)/2), color = "black",size = 1) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = Wm_L), linetype = "dashed", color="red",size = 1) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = Wm_R), color="red",size = 1) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = WB_left), linetype = "dashed", color="turquoise",size = 1) +
  geom_line(data = datasubset2.4, aes(x=(from+500)/1000000, y = WB_right), color="turquoise",size = 1) +
  ggtitle("ptg000001l") +
  theme_bw() +
  scale_color_manual(values = c("W. paniculata" = "black", "W. multiflora"="red", "W. thyrsiflora" = "blue", "W. brachyandra" = "turquoise"), name = "Species") +
  theme(text = element_text(size = 25)) 
plot(p2.4)

dev.off()


plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = c(substitute(paste(italic("W. paniculata"))), substitute(paste(italic("W. multiflora"))), substitute(paste(italic("W. thyrsiflora"))), substitute(paste(italic("W. brachyandra")))), pch=16, pt.cex=1.5, cex=1.5, bty='n', col = c("black", "red", "blue", "turquoise"))
mtext("Species", at=0.05, cex=2)


pdf(args[2], width = 20, height = 24)
plot_grid(top, mid1, mid2, mid3, mid4, mid5, mid6, mid7, mid8, mid9, mid10, mid11, mid12, mid13, bottom, ncol=1, align = "v")
dev.off()



setwd("Q:/Michael Lenhard - HFSP Shared folder/Universität Potsdam/Michael Lenhard - HFSP Shared folder/Wachendorfia/Wachendorfia_DNA-Seq_essentials")
W.q10.ptg000001l.1kb.norm.LtoR <- data.frame(read.table(file = "W.q10.ptg000001l.1kb.norm.LtoR", header = TRUE, row.names = NULL))
W.q10.ptg000001l.10kb.norm.LtoR <- data.frame(read.table(file = "W.q10.ptg000001l.10kb.norm.LtoR", header = TRUE, row.names = NULL))


W.1 <- subset(W.q10.ptg000001l.1kb.norm.LtoR, W.q10.ptg000001l.1kb.norm.LtoR$from >= 23000000 & W.q10.ptg000001l.1kb.norm.LtoR$from < 23600000)

png(file="Wachendorfia_pools_LtoR.1kb.ptg000001l_23000000-23600000.png", width=2000, height=1000)

pw.1 <- ggplot(W.1, aes(x=(from+500)/1000000)) +
  ylim(0,4) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = W.1, aes(x=(from+500)/1000000, y = Wb, color="W. brachyandra"),size = 1) +
  geom_line(data = W.1, aes(x=(from+500)/1000000, y = Wm, color ="W. multiflora"),size = 1) +
  geom_line(data = W.1, aes(x=(from+500)/1000000, y = Wp, color="W. paniculata"), size = 1) +
  geom_line(data = W.1, aes(x=(from+500)/1000000, y = Wt, color="W. thyrsiflora"), size = 1) +
  ggtitle("ptg000001l") +
  theme_bw() +
  scale_color_manual("", 
                     breaks = c("W. brachyandra", "W. multiflora", "W. paniculata", "W. thyrsiflora"),
                     values = c("black", "red", "blue", "green")) +
    theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) 
plot(pw.1)
dev.off()


W.2 <- subset(W.q10.ptg000001l.1kb.norm.LtoR, W.q10.ptg000001l.1kb.norm.LtoR$from >= 23150000 & W.q10.ptg000001l.1kb.norm.LtoR$from < 23450000)

pdf(file="Wachendorfia_pools_LtoR.1kb.ptg000001l_23150000-23450000.pdf", width=100, height=60)

ggplot(W.2, aes(x=(from+500)/1000000)) +
  ylim(0,2) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = W.2, aes(x=(from+500)/1000000, y = Wb, color="W. brachyandra"),size = 1) +
  geom_line(data = W.2, aes(x=(from+500)/1000000, y = Wm, color ="W. multiflora"),size = 1) +
  geom_line(data = W.2, aes(x=(from+500)/1000000, y = Wp, color="W. paniculata"), size = 1) +
  geom_line(data = W.2, aes(x=(from+500)/1000000, y = Wt, color="W. thyrsiflora"), size = 1) +
  ggtitle("ptg000001l") +
  theme_bw() +
  scale_color_manual("", 
                     breaks = c("W. brachyandra", "W. multiflora", "W. paniculata", "W. thyrsiflora"),
                     values = c("black", "red", "blue", "green")) +
    theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) 

dev.off()



W.3 <- subset(W.q10.ptg000001l.10kb.norm.LtoR, W.q10.ptg000001l.10kb.norm.LtoR$from >= 23000000 & W.q10.ptg000001l.10kb.norm.LtoR$from < 23600000)
write.csv(W.3, file = "Wachendorfia_pools_LtoR.10kb.ptg000001l_23000000-23600000.csv", row.names = FALSE)
pdf(file="Wachendorfia_pools_LtoR.10kb.ptg000001l_23000000-23600000.pdf", width=25, height=15)

pw.3 <- ggplot(W.3, aes(x=(from+5000)/1000000)) +
  ylim(0,2.5) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = W.3, aes(x=(from+5000)/1000000, y = Wb, color="W. brachyandra"),size = 1) +
  geom_line(data = W.3, aes(x=(from+5000)/1000000, y = Wm, color ="W. multiflora"),size = 1) +
  geom_line(data = W.3, aes(x=(from+5000)/1000000, y = Wp, color="W. paniculata"), size = 1) +
  geom_line(data = W.3, aes(x=(from+5000)/1000000, y = Wt, color="W. thyrsiflora"), size = 1) +
  ggtitle("ptg000001l") +
  theme_bw() +
  scale_color_manual("", 
                     breaks = c("W. brachyandra", "W. multiflora", "W. paniculata", "W. thyrsiflora"),
                     values = c("black", "red", "blue", "green")) +
    theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) 
plot(pw.3)
dev.off()

W.q10.ptg000001l.50kb.norm.LtoR <- data.frame(read.table(file = "W.q10.ptg000001l.50kb.norm.LtoR.txt", header = TRUE, row.names = NULL))

W.4 <- subset(W.q10.ptg000001l.50kb.norm.LtoR, W.q10.ptg000001l.50kb.norm.LtoR$from >= 22700000 & W.q10.ptg000001l.50kb.norm.LtoR$from < 23900000)
write.csv(W.4, file = "Wachendorfia_pools_LtoR.50kb.ptg000001l_22700000-23900000.csv", row.names = FALSE)
pdf(file="Wachendorfia_pools_LtoR.50kb.ptg000001l_22700000-23900000.pdf", width=25, height=15)

pw.4 <- ggplot(W.4, aes(x=(from+25000)/1000000)) +
  ylim(0,2.5) +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  theme(axis.text.x = element_text(size =20)) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wb, color="W. brachyandra"),size = 1) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wm, color ="W. multiflora"),size = 1) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wp, color="W. paniculata"), size = 1) +
  geom_line(data = W.4, aes(x=(from+25000)/1000000, y = Wt, color="W. thyrsiflora"), size = 1) +
  ggtitle("ptg000001l") +
  theme_bw() +
  scale_color_manual("", 
                     breaks = c("W. brachyandra", "W. multiflora", "W. paniculata", "W. thyrsiflora"),
                     values = c("black", "red", "blue", "green")) +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.text=element_text(size=30))+
  theme(text = element_text(size = 25)) 
plot(pw.4)
dev.off()
