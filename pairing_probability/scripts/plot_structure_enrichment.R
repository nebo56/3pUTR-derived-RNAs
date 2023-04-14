library("ggplot2")
library("cowplot")

args<-commandArgs(TRUE)

# read normalised table
structure <- read.table(args[1], sep=",", header = FALSE)

# store figure as pdf
pdf(args[2])

# plot
library(ggplot2)
ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 10, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(c(-75:75), as.numeric(structure))) + 
  ggtitle("pairing probability around 3'UTR CAGE peaks") + 
  xlab("position relative to CAGE peak") + 
  ylab("normalised coverage of pairing probability") + 
  theme(text=element_text(size=10),axis.text=element_text(size=10), axis.title=element_text(size=10,face="plain")) + 
  scale_colour_manual(values=c("#D81C1E","#F7AD64","#ADD4A2","#2C84BC")) +
  scale_x_continuous(limits = c(-60, 60)) +
  scale_y_continuous(limits = c(0.4, 0.8)) 

dev.off()