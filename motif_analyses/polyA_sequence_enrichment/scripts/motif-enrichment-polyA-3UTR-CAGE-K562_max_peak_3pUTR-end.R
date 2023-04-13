library("ggplot2")
library("smoother")
library("cowplot")

smoothing_window <- 1

# polyA sequence enrichment

# max peak in CAGE 3'UTR (excluding 200bp into 3'UTRs)
utr3.2 <- read.table("CAGE-Hs-K562-ENCSR000CJN-rep1-rep2.max_peak.filtered.3UTR-flanked100.polyA.enrichment.txt", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
utr3.2 <- t(utr3.2)
colnames(utr3.2)[1] <- "coverage"
utr3.2 <- as.data.frame(utr3.2)
utr3.2$map <- c(-100:100)
utr3.2$smooth <- smth(utr3.2$coverage, window = smoothing_window, method = "gaussian")

# 3'UTR end
utr3.4 <- read.table("gencode.v27.basic.3UTR-End-flanked100.polyA.enrichment.txt", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
utr3.4 <- t(utr3.4)
colnames(utr3.4)[1] <- "coverage"
utr3.4 <- as.data.frame(utr3.4)
utr3.4$map <- c(-100:100)
utr3.4$smooth <- smth(utr3.4$coverage, window = smoothing_window, method = "gaussian")

# plot
pdf("polyAsequence_enrichment_relative_to_3pUTR_CAGE_peak_and_3'UTR-end.pdf")

tans <- 1.0
g_size <- 0.4
adj <- 0.4

ggRNAmap.smooth.3UTR <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 10, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(utr3.3$map, as.vector(utr3.3$smooth), colour="dominant 3'UTR CAGE"),size=g_size, alpha=tans) +
  geom_line(aes(utr3.4$map, as.vector(utr3.4$smooth), colour="3'UTR end"),size=g_size, alpha=tans) + 
  ggtitle("Canonical PolyA signal enrichment\n relative to 3'UTR end maximum 3'UTR CAGE peak") + 
  xlab("position relative to 3'UTR end and 3'UTR CAGE peak") + 
  ylab("normalised motif coverage") + 
  theme(text=element_text(size=10),axis.text=element_text(size=10), axis.title=element_text(size=10,face="plain")) + 
  scale_colour_manual(values=c("#42476E", "#439915", "#FF8100", "#D94D3D", "#F2CC17", "#CFB7A2", "#3F8782")) + 
  scale_x_continuous(limits = c(-95, 95)) 
ggRNAmap.smooth.3UTR

dev.off()