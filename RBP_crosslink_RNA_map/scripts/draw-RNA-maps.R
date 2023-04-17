library("ggplot2")
library("smoother")
library("cowplot")

args<-commandArgs(TRUE)
smoothing_window <- 3   # size of the smoothing window

# RNA map of coverage in region -100 to +100 bps relative to the target size.
RBP.target <- read.table(args[1], sep='\t')

RBP.target.plus <- RBP.target[which(RBP.target$V6 == "+"),] # strand needs to be seperated for the correct counting
RBP.target.minus <- RBP.target[which(RBP.target$V6 == "-"),]

RBP.target.plus.sum <- aggregate(RBP.target.plus$V8, list(RBP.target.plus$V7), FUN=sum)
RBP.target.minus.sum <- aggregate(RBP.target.minus$V8, list(RBP.target.minus$V7), FUN=sum)
RBP.target.minus.sum$Group.1 <- c(201:1)

RBP.target.sum <- merge(RBP.target.plus.sum, RBP.target.minus.sum, by="Group.1")
RBP.target.sum$x <- RBP.target.sum$x.x + RBP.target.sum$x.y

RBP.target.sum$norm <- RBP.target.sum$x / nrow(unique(RBP.target[c("V5")])) # normalise by the number of 
RBP.target.sum$smooth <- smth(RBP.target.sum$norm, window = smoothing_window, method = "gaussian") # smoothing the curve
RBP.target.sum$map <- c(-100:100)
remove(RBP.target)



### plot ###
pdf(paste(args[2],".pdf", sep=""))

tans <- 1.0
g_size <- 0.4

ggRNAmap <- ggplot() + theme_bw() + 
  geom_line(aes(RBP.target.sum$map, as.vector(RBP.target.sum$smooth), colour="RBP"),size=g_size, alpha=tans) +
  ggtitle("RNA map") + 
  xlab("position relative to the peak") + 
  ylab("coverage of crosslink enrichment") + 
  theme(text=element_text(size=8),axis.text=element_text(size=8), axis.title=element_text(size=8,face="plain")) + 
  scale_x_continuous(limits = c(-90, 90))
ggRNAmap



