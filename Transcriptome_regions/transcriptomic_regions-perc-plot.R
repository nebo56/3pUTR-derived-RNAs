library("ggplot2")
library("smoother")
library(RColorBrewer)
library(reshape2)

ratio <- read.table("CAGE-Hs-HepG2-HeLa-K562-samples-unique_positions-genomic_regions.txt", sep = "\t", header=TRUE)
colnames(ratio) <- c("sample", "read.number", "perc.3UTR", "perc.5UTR", "perc.CDS", "perc.Intronic")

# reshape the table using melt function
ratiom <- melt(ratio[,c("sample", "read.number", "perc.3UTR", "perc.5UTR", "perc.CDS", "perc.Intronic")],id.vars = 1)

# ratio of genomic regions
gg.ratio <- ggplot(ratiom[which(ratiom$variable != "read.number"),],aes(x = variable,y = value)) + theme_bw()
gg.ratio + geom_bar(stat="identity", width = 0.7, aes(fill = sample), position = "dodge") + 
  scale_y_continuous(name="%", labels = scales::comma) + 
  ggtitle("CAGE samples - uniq positions - genomic ratio") +
  scale_fill_brewer(palette="Spectral")

#read.count
gg.read.number <- ggplot(ratiom[which(ratiom$variable == "read.number"),],aes(x = variable,y = value)) + theme_bw()
gg.read.number + geom_bar(stat="identity", width = 0.7, aes(fill = sample), position = "dodge") + 
  scale_y_continuous(name="total number of unique positions", labels = scales::comma) + 
  ggtitle("CAGE samples - uniq positions - genomic ratio") +
  scale_fill_brewer(palette="Spectral")


# Stacked Percent
ggplot(ratiom[which(ratiom$variable != "read.number"),], aes(fill=variable, y=value, x=sample)) + theme_bw() +
  geom_bar( stat="identity", position="fill") + 
  coord_flip() +
  ggtitle("CAGE samples uniq positions  - genomic ratio") +
  scale_fill_brewer(palette="Spectral") +
  xlab("ratio") 

