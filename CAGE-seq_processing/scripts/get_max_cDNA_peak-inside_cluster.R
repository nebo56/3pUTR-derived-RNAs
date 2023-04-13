args<-commandArgs(TRUE)

#peaks <- read.table("CAGE-Hs-K562-ENCSR000CJN-rep1-rep2-merged-cDNAstart-sum-2counts_limit-UTR3-merged15nt-peaks.bed", header=FALSE)
peaks <- read.table(args[1], header=FALSE)

colnames(peaks)[1] <- "read.chr"
colnames(peaks)[2] <- "read.start"
colnames(peaks)[3] <- "read.end"
colnames(peaks)[4] <- "read.count"
colnames(peaks)[5] <- "read.strand"
colnames(peaks)[6] <- "cluster.chr"
colnames(peaks)[7] <- "cluster.start"
colnames(peaks)[8] <- "cluster.end"
colnames(peaks)[11] <- "cluster.strand"

max.peaks <-aggregate(peaks$read.count, list(peaks$cluster.chr, peaks$cluster.start, peaks$cluster.end, peaks$cluster.strand), FUN=max, na.rm=TRUE)
colnames(max.peaks)[1] <- "cluster.chr"
colnames(max.peaks)[2] <- "cluster.start"
colnames(max.peaks)[3] <- "cluster.end"
colnames(max.peaks)[4] <- "cluster.strand"
colnames(max.peaks)[5] <- "read.count"

max.peaks.positions <- merge(max.peaks, peaks, by=c("cluster.chr","cluster.start","cluster.end","cluster.strand","read.count"), all.x=TRUE)

max.peaks.positions <- max.peaks.positions[c("read.chr", "read.start","read.end","read.count","V9","read.strand")]

write.table(max.peaks.positions,  args[2], row.names=FALSE, sep="\t", quote=FALSE)

# median_max_pos <-aggregate(max.peaks.positions$read.start, list(max.peaks.positions$cluster.chr, max.peaks.positions$cluster.start, max.peaks.positions$cluster.end, max.peaks.positions$cluster.strand, max.peaks.positions$read.count), FUN=median, na.rm=TRUE)
# colnames(median_max_pos)[1] <- "cluster.chr"
# colnames(median_max_pos)[2] <- "cluster.start"
# colnames(median_max_pos)[3] <- "cluster.end"
# colnames(median_max_pos)[4] <- "cluster.strand"
# colnames(median_max_pos)[5] <- "read.peak.count"
# colnames(median_max_pos)[6] <- "read.peak.position"
# median_max_pos$read.peak.position <- as.integer(round(median_max_pos$read.peak.position, digits = 0))
# 
# write.table(median_max_pos, args[2], row.names=FALSE, sep="\t", quote=FALSE)


