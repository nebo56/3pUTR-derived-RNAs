library(CAGEr)
library("BSgenome.Hsapiens.UCSC.hg38")
library(GenomicRanges)
seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "ENSEMBL"

# for more detailed information follow instructions from CAGEr bioconductor site: http://www.bioconductor.org/packages/devel/bioc/vignettes/CAGEr/inst/doc/CAGEexp.html

# read mapped by bowtie2 BAM CAGE samples. Folder path needs to be specified here
bam_file <- list.files("./BAMs-bowtie2/", full.names=T, pattern="*.bam")
bam_file

ce <- CAGEexp( genomeName     = "BSgenome.Hsapiens.UCSC.hg38"
               , inputFiles     = bam_file
               , inputFilesType = "bam"
               , sampleLabels   = c('rep1','rep2')
)

# create CAGE / TSS data 
ce <- getCTSS(ce, correctSystematicG=FALSE,removeFirstG=TRUE) # remove first G can be TRUE


CTSStagCountSE(ce)
CTSScoordinatesGR(ce)
CTSStagCountDF(ce)
CTSStagCountGR(ce, 1)  # GRanges for one sample with expression count.
sampleLabels(ce)
corr.m <- plotCorrelation2( ce, samples = "all"
                            , tagCountThreshold = 1, applyThresholdBoth = FALSE
                            , method = "pearson")

corr.m <- plotCorrelation2( ce, samples = "all"
                            , tagCountThreshold = 1, applyThresholdBoth = FALSE
                            , method = "pearson")

ce <- CAGEr::normalizeTagCount(ce, method = "powerLaw",
                               fitInRange = c(10, 10^4), alpha = 1.18, T = 10^6)

ce <- CAGEr::clusterCTSS(object = ce,
                         threshold = 1,
                         thresholdIsTpm = TRUE,
                         nrPassThreshold = 1,
                         method = "distclu",
                         maxDist = 20,
                         removeSingletons = TRUE,
                         keepSingletonsAbove = 5)

# create objects
CTSStagCountSE(ce)
CTSScoordinatesGR(ce)
CTSStagCountDF(ce)
CTSStagCountGR(ce, 1)  # GRanges for one sample with expression count.
sampleLabels(ce)

# Correlation between samples
corr.m <- plotCorrelation2( ce, samples = "all"
                            , tagCountThreshold = 1, applyThresholdBoth = FALSE
                            , method = "pearson")


# Library size
librarySizes(ce)

# Reverse cumulative distribution of CAGE tags
plotReverseCumulatives(ce, fitInRange = c(5, 1000), onePlot = TRUE)