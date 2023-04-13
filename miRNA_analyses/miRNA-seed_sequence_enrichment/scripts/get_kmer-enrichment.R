args<-commandArgs(TRUE)

kmers <- read.table(args[1], header = TRUE)
kmers.control <- read.table(args[2], header = TRUE)

merged <- merge(kmers, kmers.control, by="kmer", all.x = TRUE)
merged$enrichment <- merged$percentage_of_inclusion.x / merged$percentage_of_inclusion.y

#summary(merged$enrichment)
merged <- merged[order(-merged$enrichment),]

# filter 1.5 enriched
merged.filtered <- merged[which(merged$enrichment >= 1.5),]

# export as table
write.table(merged.filtered, paste(args[3], "_1.5_enriched.tab", sep = ""), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#export as fasta
write.table(merged.filtered$kmer, paste(args[3], "_1.5_enriched.fa", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
