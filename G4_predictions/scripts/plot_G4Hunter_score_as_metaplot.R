library("ggplot2")
library("smoother")
library("cowplot")

# in order to get the right format you need to remove specific lines
# $cat CAGE-Hs-K562-ENCSR000CJN-rep1-rep2-W25-S1.2.txt | grep -v "Start" | grep -v '>' > test.txt
# and add header "Start 	 End 	 Sequence	 Length 	 Score"

args<-commandArgs(TRUE)

# inport G4scores
G4score <- read.table(args[1], header = FALSE, sep = '\t')
colnames(G4score) <- c("Start", "End", "Sequence", "Length", "Score")
G4score$Centre <- G4score$Start + round((G4score$End - G4score$Start) / 2)
G4score.sum <- aggregate(G4score$Score, list(G4score$Centre), FUN=sum)
colnames(G4score.sum) <- c("position", "G4H.score.sum")

# plot scores as metaplot

pdf(args[2])

tans <- 1.0
g_size <- 0.4
adj <- 0.4

ggRNAmap <- ggplot(rep1.all, aes(x=position-100, y=G4H.score.sum, color=region)) + theme_cowplot() + theme(plot.title = element_text(size = 10, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(size=g_size, alpha=tans) +
  ggtitle("G4Hunter score plot: CAGE") + 
  xlab("position relative to CAGE peak") + 
  ylab("G4-Hunter score sum") + 
  theme(text=element_text(size=10),axis.text=element_text(size=10), axis.title=element_text(size=10,face="plain")) + 
  scale_colour_manual(values=c("#42476E", "#4F80E1", "#439915", "#FF8100", "#D94D3D", "#F2CC17", "#CFB7A2", "#3F8782")) +
  scale_x_continuous(limits = c(-100, 100)) 
ggRNAmap


