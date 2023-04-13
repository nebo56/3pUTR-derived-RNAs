library("ggplot2")
library(data.table)

args<-commandArgs(TRUE)

# For each coverage will first be seperated by strand, then summed positions in the region as metaplot, followed by normalisation and indexing relative to the peak

# RNAseq cDNA-ends: 
######################################
RNAseq.1 <- fread(args[1], sep='\t')
RNAseq.1.plus <- RNAseq.1[which(RNAseq.1$V6 == "+"),]
RNAseq.1.minus <- RNAseq.1[which(RNAseq.1$V6 == "-"),]
RNAseq.1.plus.sum <- aggregate(RNAseq.1.plus$V8, list(RNAseq.1.plus$V7), FUN=sum)
RNAseq.1.minus.sum <- aggregate(RNAseq.1.minus$V8, list(RNAseq.1.minus$V7), FUN=sum)
RNAseq.1.minus.sum$Group.1 <- c(151:1)
RNAseq.1.sum <- merge(RNAseq.1.plus.sum, RNAseq.1.minus.sum, by="Group.1")
RNAseq.1.sum$x <- RNAseq.1.sum$x.x + RNAseq.1.sum$x.y
RNAseq.1.sum$norm <- RNAseq.1.sum$x / length(unique(RNAseq.1$V5)) # normalise by the number of positions
RNAseq.1.sum$map <- c(-75:75)
remove(RNAseq.1)

# RNAseq cDNA-starts: 
######################################
RNAseq.3 <- fread(args[2], sep='\t')
RNAseq.3.plus <- RNAseq.3[which(RNAseq.3$V6 == "+"),]
RNAseq.3.minus <- RNAseq.3[which(RNAseq.3$V6 == "-"),]
RNAseq.3.plus.sum <- aggregate(RNAseq.3.plus$V8, list(RNAseq.3.plus$V7), FUN=sum)
RNAseq.3.minus.sum <- aggregate(RNAseq.3.minus$V8, list(RNAseq.3.minus$V7), FUN=sum)
RNAseq.3.minus.sum$Group.1 <- c(151:1)
RNAseq.3.sum <- merge(RNAseq.3.plus.sum, RNAseq.3.minus.sum, by="Group.1")
RNAseq.3.sum$x <- RNAseq.3.sum$x.x + RNAseq.3.sum$x.y
RNAseq.3.sum$norm <- RNAseq.3.sum$x / length(unique(RNAseq.3$V5))
RNAseq.3.sum$map <- c(-75:75)
remove(RNAseq.3)

#### plot ####

# save to pdf
pdf(args[3],         # file name
    bg = "white",          # background color
    colormodel = "cmyk",    # color model 
    paper = "A4")          # paper size

# set line paramters
tans <- 1.0
g_size <- 0.4
adj <- 0.4

# plot figure
ggRNAseq.3utr <- ggplot() + theme_bw() + 
  geom_line(aes(RNAseq.1.sum$map, as.vector(RNAseq.1.sum$norm), colour="cDNAends"),size=g_size, alpha=tans) + 
  geom_line(aes(RNAseq.2.sum$map, as.vector(RNAseq.2.sum$norm), colour="cDNAstarts"),size=g_size, alpha=tans) + 
  ggtitle("RNAseq cDNA start/end coverage around 3'UTR CAGE peaks") + 
  xlab("position relative to CAGE peak") + 
  ylab("coverage of read starts/ends") + 
  theme(text=element_text(size=8),axis.text=element_text(size=8), axis.title=element_text(size=8,face="plain")) + 
  scale_colour_manual(values=c("#08519c", "#B2463B")) +
  scale_x_continuous(limits = c(-75, 75)) 

# closing the graphical device
dev.off()
