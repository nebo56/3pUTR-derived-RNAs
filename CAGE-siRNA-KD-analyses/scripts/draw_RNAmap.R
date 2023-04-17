library("ggplot2")
library("smoother")
library("data.table")

smoothing_window <- 5
args<-commandArgs(TRUE)

# rep1
######################################
#CAGE.1 <- fread("../ISL1-Start.bed-flank100-ISL1_batch3_rep1.ctss.bed.coverage.bed", sep='\t')
CAGE.1 <- fread(args[1], sep='\t')
CAGE.1$norm <- CAGE.1$V8 #/ mean(CAGE.1$V8)

# rep2
######################################
#CAGE.2 <- fread("../ISL1-Start.bed-flank100-ISL1_batch3_rep2.ctss.bed.coverage.bed", sep='\t')
CAGE.2 <- fread(args[2], sep='\t')
CAGE.2$norm <- CAGE.2$V8 #/ mean(CAGE.2$V8)

# rep3
######################################
#CAGE.3 <- fread("../ISL1-Start.bed-flank100-ISL1_batch3_rep3.ctss.bed.coverage.bed", sep='\t')
CAGE.3 <- fread(args[3], sep='\t')
CAGE.3$norm <- CAGE.3$V8 #/ mean(CAGE.3$V8)

# rep4
######################################
#CAGE.4 <- fread("../TCF25-Start.bed-flank30-TCF25_batch1A_rep3.ctss.bed.coverage.bed", sep='\t')
#CAGE.4 <- fread(args[4], sep='\t')
#CAGE.4$norm <- CAGE.4$V8 #/ mean(CAGE.4$V8)

#in case of minus strand
if (CAGE.1$V6[1] == "-") {
  CAGE.1$V7 <- order(-CAGE.1$V7) -2
  CAGE.2$V7 <- order(-CAGE.2$V7) -2
  CAGE.3$V7 <- order(-CAGE.3$V7) -2
  #CAGE.4$V7 <- order(-CAGE.4$V7) -2
}

red <- c("#F7827A","#ED403E","#832120","#560F17")
blue <- c("#6DCCE8","#0D7AE8","#023B85","#02195E")
grey <- c("#CACFD9","#898E99","#5A6670","#2F373D","#000000")
tans <- 0.8
g_size <- 0.5
adj <- 0.4

pdf(paste(args[4],"-CAGE-KD_relative_to_siRNA_target.pdf", sep = ""))

# gg <- ggplot() + theme_bw() + 
#   geom_line(aes(c(0:31), as.vector(CAGE.1$norm), colour="rep1"),size=g_size, alpha=tans) + 
#   geom_line(aes(c(0:31), as.vector(CAGE.2$norm), colour="rep2"),size=g_size, alpha=tans) +
#   geom_line(aes(c(0:31), as.vector(CAGE.3$norm), colour="rep3"),size=g_size, alpha=tans) +
#   ggtitle(paste(args[4], "CAGE-KD signal relative to siRNA target", sep=": ")) + 
#   xlab("position relative to target start") + 
#   ylab("normalised coverage of CAGE ctss") + 
#   theme(text=element_text(size=8),axis.text=element_text(size=8), axis.title=element_text(size=8,face="plain")) #+ 
#   #scale_x_continuous(limits = c(-5, 35)) #+
#   #scale_y_continuous(limits = c(0, 8)) +
#   #scale_colour_manual(values=c(blue[1], blue[2], blue[3], blue[4], red[1], red[2], red[3], red[4]))
# gg


g1 <- ggplot(data=CAGE.1, aes(x=V7, y=V8, fill="rep1")) + theme_bw() + 
  geom_bar(stat="identity", position=position_dodge(), fill = "#6DCCE8") +
  ggtitle(paste(args[4], "CAGE-KD signal relative to siRNA target - rep1", sep=": ")) + 
  xlab("position relative to target start") + 
  ylab("coverage of CAGE ctss") 

g2 <- ggplot(data=CAGE.2, aes(x=V7, y=V8, fill="rep2")) + theme_bw() + 
  geom_bar(stat="identity", position=position_dodge(), fill = "#0D7AE8") +
  ggtitle(paste(args[4], "CAGE-KD signal relative to siRNA target - rep2", sep=": ")) + 
  xlab("position relative to target start") + 
  ylab("coverage of CAGE ctss") 

g3 <- ggplot(data=CAGE.3, aes(x=V7, y=V8, fill="rep3")) + theme_bw() + 
  geom_bar(stat="identity", position=position_dodge(), fill = "#023B85") +
  ggtitle(paste(args[4], "CAGE-KD signal relative to siRNA target - rep3", sep=": ")) + 
  xlab("position relative to target start") + 
  ylab("coverage of CAGE ctss") 

#g4 <- ggplot(data=CAGE.4, aes(x=V7, y=V8, fill="rep3")) + theme_bw() + 
#  geom_bar(stat="identity", position=position_dodge(), fill = "#023B85") +
#  ggtitle(paste(args[4], "CAGE-KD signal relative to siRNA target - rep4", sep=": ")) + 
#  xlab("position relative to target start") + 
#  ylab("coverage of CAGE ctss") 

library(gridExtra)
#grid.arrange(g1, g2, g3, g4, ncol = 1)
grid.arrange(g1, g2, g3, ncol = 1)

# merge all into one !!! needs to be separated by strand or to the position !!!
CAGE.all <- rbind(CAGE.1, CAGE.2, CAGE.3)
CAGE.sum <- aggregate(CAGE.all$V8, list(CAGE.all$V7), FUN=sum)
colnames(CAGE.sum) <- c("position", "read.count")

g.all <- ggplot(data=CAGE.sum, aes(x=position, y=read.count, fill="all rep sum")) + theme_bw() + 
  geom_bar(stat="identity", position=position_dodge(), fill = "#023B85") +
  ggtitle(paste(args[4], "CAGE-KD signal relative to siRNA target - all rep", sep=": ")) + 
  xlab("position relative to target start") + 
  ylab("coverage of CAGE ctss") 

grid.arrange(g.all, ncol = 1)
CAGE.sum$sum.all.norm <- CAGE.sum$read.count / sum(CAGE.sum$read.count)

write.table(cbind((args[4]), t(CAGE.sum[,3])), paste(args[4],"-sum.norm.counts.txt", sep = ""), sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(cbind((args[4]), t(CAGE.sum[,2])), paste(args[4],"-sum.counts.txt", sep = ""), sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
# CAGE.1$replicate <- "rep1"
# CAGE.2$replicate <- "rep2"
# CAGE.3$replicate <- "rep3"
# 
# CAGE.all.rep <- rbind(CAGE.1, CAGE.2, CAGE.3)
#   
# ggplot(data=CAGE.all.rep, aes(x=V7, y=V8, fill=replicate)) + theme_bw() + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   scale_x_continuous(limits = c(50, 80))

# old

# rep2
######################################
#CAGE.2 <- fread("../ISL1-Start.bed-flank50-ISL1_batch3_rep2.ctss.bed.coverage.bed", sep='\t')
#CAGE.2.sum <- aggregate(CAGE.2$V8, list(CAGE.2$V7), FUN=sum)
#CAGE.2.sum$norm <- CAGE.2.sum$x / mean(CAGE.2.sum$x)
#CAGE.2.sum$smooth <- smth(CAGE.2.sum$norm, window = smoothing_window, method = "gaussian") #SMOOTHING
#CAGE.2.sum$map <- c(-30:31)
#remove(CAGE.2)
