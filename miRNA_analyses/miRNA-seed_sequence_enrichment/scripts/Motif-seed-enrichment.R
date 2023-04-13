library("ggplot2")
library("smoother")
library("cowplot")

smoothing_window <- 1
args<-commandArgs(TRUE)

#coverage <- read.table("../ams-a2-hela_trimmed_single-sum-TARGET.bed.seed.coverage.tab", header=FALSE, sep="\t")
coverage <- read.table(args[1], header=FALSE, sep="\t")
colnames(coverage) <- c("coverage","motif","name")
coverage$position <- rep(c(-30:30),length(unique(coverage$motif)))

tans <- 1.0
g_size <- 0.4
adj <- 0.4

pdf(paste(args[1], ".pdf", sep = ""))

gg <- ggplot(coverage, aes(x=position, y=coverage)) + theme_cowplot() + 
  geom_line(aes(colour=motif),size=g_size, alpha=tans) + 
  ggtitle("6-mer enrichment around Ago2 (eiCLIP) binding") + 
  xlab("position relative to Ago2 crosslinking position") + 
  ylab("normalised motif coverage")  
gg


# sum of all motifs
coverage.sum <- aggregate(coverage$coverage, list(coverage$position), FUN = sum)
colnames(coverage.sum) <- c("position", "summed.coverage")

gg.all <- ggplot(coverage.sum, aes(x=position, y=summed.coverage)) + theme_cowplot() + 
  geom_line(size=g_size, alpha=tans) + 
  ggtitle("miR-seed enrichment around Ago2 (eiCLIP) binding") + 
  xlab("position relative to Ago2 crosslinking position") + 
  ylab("normalised motif coverage")  
gg.all


# heatmap
library(ggplot2)
library(RColorBrewer)
library(reshape2)

hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  
ggplot(coverage, aes(x = position, y = as.factor(name))) + theme_bw() + 
  geom_tile(aes(fill=coverage)) +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(100)) +
  ggtitle("miR-seed enrichment around Ago2 (eiCLIP) binding") + 
  xlab("position relative to Ago2 crosslinking position") + 
  ylab("normalised motif coverage")  



# order 
coverage.mean <- aggregate(coverage$coverage, list(coverage$motif, coverage$name), FUN = mean)
colnames(coverage.mean) <- c("motif","name","mean")
coverage.mean$name.motif <- paste(coverage.mean$name, coverage.mean$motif)
coverage.mean <- coverage.mean[order(-coverage.mean$mean),]
coverage.mean$index <- c(1:nrow(coverage.mean))
coverage.ordered <- merge(coverage,coverage.mean, by="motif")


hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  
ggplot(coverage.ordered, aes(x = position, y = as.factor(-index))) + theme_bw() + 
  geom_tile(aes(fill=coverage)) +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(100)) +
  scale_y_discrete(labels=as.factor(rev(coverage.mean$name.motif))) +
  ggtitle("miR-seed enrichment around Ago2 (eiCLIP) binding") + 
  xlab("position relative to Ago2 crosslinking position") + 
  ylab("normalised motif coverage")  


# threshold !!!!!!!! something is not right- TTTTT is not on the top anymore 1!!!!!
#coverage.mean <- coverage.mean[order(-coverage.mean$mean),]
coverage.mean.top100 <- head(coverage.mean,100)
coverage.mean.top100$index <- c(1:nrow(coverage.mean.top100))
coverage.ordered <- merge(coverage,coverage.mean.top100, by="motif")


m.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  
ggplot(coverage.ordered, aes(x = position, y = as.factor(-index))) + theme_bw() + 
  geom_tile(aes(fill=coverage)) +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(100)) +
  scale_y_discrete(labels=as.factor(rev(coverage.mean.top100$name.motif))) +
  ggtitle("miR-seed enrichment around Ago2 (eiCLIP) binding") + 
  xlab("position relative to Ago2 crosslinking position") + 
  ylab("normalised motif coverage")  


# # threshold - order by max peak
# coverage.mean <- aggregate(coverage$coverage, list(coverage$motif), FUN = max)
# colnames(coverage.mean) <- c("motif","max")
# coverage.max <- coverage.mean[order(-coverage.mean$max),]
# coverage.max.top100 <- head(coverage.max,100)
# coverage.max.top100$index <- c(1:nrow(coverage.max.top100))
# coverage.ordered <- merge(coverage,coverage.max.top100, by="motif")
# 
# m.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  
# ggplot(coverage.ordered, aes(x = position, y = as.factor(-index))) + theme_bw() + 
#   geom_tile(aes(fill=coverage)) +
#   coord_equal() +
#   scale_fill_gradientn(colours = hm.palette(100)) +
#   scale_y_discrete(labels=as.factor(coverage.mean$motif)) +
#   ggtitle("miR-seed enrichment around Ago2 (eiCLIP) binding") + 
#   xlab("position relative to Ago2 crosslinking position") + 
#   ylab("normalised motif coverage")  
# 


# threshold top 10
coverage.mean.top10 <- head(coverage.mean,10)
coverage.mean.top10$index <- c(1:nrow(coverage.mean.top10))
coverage.ordered <- merge(coverage,coverage.mean.top10, by="motif")

# RNA map
gg <- ggplot(coverage.ordered, aes(x=position, y=coverage)) + theme_cowplot() + 
  geom_line(aes(colour=name.motif),size=g_size, alpha=tans) + 
  ggtitle("miR-seed enrichment around Ago2 (eiCLIP) binding") + 
  xlab("position relative to Ago2 crosslinking position") + 
  ylab("normalised motif coverage")  
gg

