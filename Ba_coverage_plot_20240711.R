install.packages("ggplot2")
install.packages("cowplot")

library(ggplot2)
library(cowplot)


setwd("//141.89.108.211/genetics-shared/hxue/Barberetta genome assembly and coverage analysis//")

Ba_Ba700mhap2.10kb.norm <- data.frame(read.table(file = "Ba_Ba700mhap2.unique.10kb.regions.norm", header = TRUE, row.names = NULL))
Ba_Ba700mhap2.50kb.norm <- data.frame(read.table(file = "Ba_Ba700mhap2.unique.50kb.regions.norm", header = TRUE, row.names = NULL))
Ba_Ba700mhap2.100kb.norm <- data.frame(read.table(file = "Ba_Ba700mhap2.unique.100kb.regions.norm", header = TRUE, row.names = NULL))


datasubset1 <- subset(Ba_Ba700mhap2.100kb.norm, Ba_Ba700mhap2.100kb.norm$contig=="h2tg000059l")
png(file="Ba_Ba700mhap2.100kb.h2tg000059l.png", width=1000, height=1000)
p1 <- ggplot() +
  theme(legend.position = "bottom") +
  ylim(0,2.5) +
  geom_line(data = datasubset1, aes(x=(bis+50000)/1000000, y = Ba_L/Ba_R), linetype = "dashed", color="blue",size = 1) +
  ylab("Normalized coverage") +
  xlab("MB") +
  ggtitle("h2tg000059l") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 25)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))

plot(p1)

png(file="Ba_Ba700mhap2.100kb.h2tg000059l.png", width=1000, height=1000)
plot(p1)
dev.off()

datasubset2 <- subset(Ba_Ba700mhap2.10kb.norm, Ba_Ba700mhap2.10kb.norm$contig=="h2tg000059l" & Ba_Ba700mhap2.10kb.norm$von >= 1000000 & Ba_Ba700mhap2.10kb.norm$von <= 1350000)
p2 <- ggplot() +
  theme(legend.position = "bottom") +
  ylim(0,1.5) +
  geom_line(data = datasubset2, aes(x=(von+5000)/1000000, y = Ba_L/Ba_R), color="blue", linewidth = 1) +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_hline(yintercept=1,linetype="dashed", color = "orange") +
  geom_vline(xintercept=c(1.173447,1.178435,1.227045,1.227309),linetype="dashed") +
  ylab("Left to right coverage ratio") +
  xlab("MB") +
  ggtitle("h2tg000059l") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 25))

plot(p2)

png(file="Ba_Ba700mhap2.10kb.h2tg000059l_1-1.35MB.png", width=1000, height=500)
plot(p2)
dev.off()

datasubset1 <- subset(Ba_Ba700mhap2.100kb.norm, Ba_Ba700mhap2.100kb.norm$contig=="h2tg000059l")
png(file="Ba_Ba700mhap2.100kb.h2tg000059l.png", width=1000, height=1000)
p1 <- ggplot() +
  theme(legend.position = "bottom") +
  ylim(0,2.5) +
  geom_line(data = datasubset1, aes(x=(bis+50000)/1000000, y = Ba_L/Ba_R), color="blue",size = 1) +
  geom_vline(xintercept=2,linetype="dashed") +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  ggtitle("h2tg000059l") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 25)) 

plot(p1)

png(file="Ba_Ba700mhap2.100kb.h2tg000059l.png", width=1000, height=1000)
plot(p1)
dev.off()

datasubset2 <- subset(Ba_Ba700mhap2.10kb.norm, Ba_Ba700mhap2.10kb.norm$contig=="ptg000001l" & Ba_Ba700mhap2.10kb.norm$bis >= 21000000 & Ba_Ba700mhap2.10kb.norm$bis <= 24000000)

png(file="Ba_Ba700mhap2.10kb.ptg000001l_21000000-24010000.png", width=1000, height=1000)

p1 <- ggplot() +
  theme(legend.position = "bottom") +
  ylim(0,2.5) +
  geom_line(data = datasubset1, aes(x=(bis+50000)/1000000, y = Ba_L/Ba_R), color="blue",size = 1) +
  geom_vline(xintercept=2,linetype="dashed") +
  ylab("L:R coverage ratio") +
  xlab("MB") +
  ggtitle("h2tg000059l") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 25)) 

plot(p1)

png(file="Ba_Ba700mhap2.100kb.h2tg000059l.png", width=1000, height=1000)
plot(p1)
dev.off()

datasubset3 <- subset(Ba_Ba700mhap2.10kb.norm, Ba_Ba700mhap2.10kb.norm$contig=="ptg000001l" & Ba_Ba700mhap2.10kb.norm$bis >= 22800000 & Ba_Ba700mhap2.10kb.norm$bis <= 23800000)

png(file="Ba_Ba700mhap2.10kb.ptg000001l_22800000-23810000.png", width=1000, height=1000)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = c(substitute(paste(italic("W. paniculata"))), substitute(paste(italic("W. multiflora"))), substitute(paste(italic("W. thyrsiflora"))), substitute(paste(italic("W. brachyandra")))), pch=16, pt.cex=1.5, cex=1.5, bty='n', col = c("black", "red", "blue", "turquoise"))
mtext("Species", at=0.05, cex=2)

pdf(args[2], width = 20, height = 24)
plot_grid(top, mid1, mid2, mid3, mid4, mid5, mid6, mid7, mid8, mid9, mid10, mid11, mid12, mid13, bottom, ncol=1, align = "v")
dev.off()

