args<-commandArgs(TRUE)

peaks <- read.table(args[1], header=FALSE)

# name calumns
colnames(peaks)[1] <- "read.chr"
colnames(peaks)[2] <- "read.start"
colnames(peaks)[3] <- "read.end"
colnames(peaks)[4] <- "read.count"
colnames(peaks)[6] <- "read.strand"
colnames(peaks)[7] <- "cluster.chr"
colnames(peaks)[8] <- "cluster.start"
colnames(peaks)[9] <- "cluster.end"
colnames(peaks)[12] <- "cluster.strand"

# select max peak
max.peaks <-aggregate(peaks$read.count, list(peaks$cluster.chr, peaks$cluster.start, peaks$cluster.end, peaks$cluster.strand), FUN=max, na.rm=TRUE)
colnames(max.peaks)[1] <- "cluster.chr"
colnames(max.peaks)[2] <- "cluster.start"
colnames(max.peaks)[3] <- "cluster.end"
colnames(max.peaks)[4] <- "cluster.strand"
colnames(max.peaks)[5] <- "read.count"

max.peaks.positions <- merge(max.peaks, peaks, by=c("cluster.chr","cluster.start","cluster.end","cluster.strand","read.count"), all.x=TRUE)
max.peaks.positions <- max.peaks.positions[c("read.chr", "read.start","read.end","read.count","V10","read.strand")]

write.table(max.peaks.positions,  args[2], row.names=FALSE, sep="\t", quote=FALSE)

