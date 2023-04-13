library("ggplot2")
library("smoother")
library("cowplot")

args<-commandArgs(TRUE)
smoothing_window <- 3 # size of the smoothing window

# GGG enrichment relative to its target
GGG.cov <- read.table(args[1], header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
GGG.cov <- t(GGG.cov)
colnames(GGG.cov)[1] <- "coverage"
GGG.cov <- as.data.frame(GGG.cov)
GGG.cov$map <- c(-100:100)
GGG.cov$smooth <- smth(GGG.cov$coverage, window = smoothing_window, method = "gaussian")

# plot coverage
pdf(args[2])

tans <- 1.0
g_size <- 0.4
adj <- 0.4

gggEnrichment <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 10, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(GGG.cov$map, as.vector(GGG.cov$smooth), colour="GGG"),size=g_size, alpha=tans) + 
  ggtitle("G-enrichment around max 3'UTR CAGE peaks") + 
  xlab("position relative to CAGE peak") + 
  ylab("normalised GGGG coverage") + 
  theme(text=element_text(size=10),axis.text=element_text(size=10), axis.title=element_text(size=10,face="plain")) + 
  gggEnrichment

dev.off()
