#!/usr/bin/env -S Rscript --vanilla
library("sequenza")

args = commandArgs(trailingOnly=TRUE)

# test if there are exactly two or four arguments:
if (length(args)==2) {
    cellularity_input = NULL
    ploidy_input = NULL
} else if (length(args)==4) {
    cellularity_input = as.numeric(args[3])
    ploidy_input = as.numeric(args[4])
} else {
    stop("Exactly 2 or 4 arguments must be supplied: sample_id, sample.binned.seqz[, cellularity, ploidy].n", call.=FALSE)
}

#-----------------------------------------------------------------------------------------------
#RUN SEQUENZA, R
#-----------------------------------------------------------------------------------------------
data.file <- args[2]
seqz.data <- read.seqz(data.file)
gc.stats <- gc.sample.stats(data.file)
test <- sequenza.extract(data.file)
CP.example <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = args[1], out.dir="./")
cint <- get.ci(CP.example)

#Plot cellularity
jpeg(paste(args[1],".nitz.cellularity.jpg", sep="", collapse=NULL))
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE, likThresh=c(0.95))
dev.off()

#Call CNVs
if (is.null(cellularity_input)) {
    cellularity <- cint$max.cellularity
    ploidy <- cint$max.ploidy
} else {
    cellularity <- cellularity_input
    ploidy <- ploidy_input
}
seg_table <- read.table(paste(args[1],"_segments.txt", sep="", collapse=NULL), header = TRUE, sep = "\t", dec = ".")
avg.depth.ratio <- mean(seg_table$depth.ratio)

#Save parameters to file
cellularity
write(cellularity, file = paste(args[1],".nitz.cellularity.txt", sep="", collapse=NULL))
write(ploidy, file = paste(args[1],".nitz.ploidy.txt", sep="", collapse=NULL))
write(avg.depth.ratio, file = paste(args[1],".nitz.ave_depth.txt", sep="", collapse=NULL))

#Detect variant alleles
mut.tab <- na.exclude(do.call(rbind, test$mutations))
mut.alleles <- mufreq.bayes(mufreq = mut.tab$F, depth.ratio = mut.tab$adjusted.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)

#Detect CN variation
seg.tab <- na.exclude(do.call(rbind, test$segments))
cn.alleles <- baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)
seg.tab <- cbind(seg.tab, cn.alleles)
seg.tab

#write sequenza matrix to file, this will serve as input to loss score script's 2nd arg
write.table(seg.tab, file = paste(args[1],".nitz.copynumber_calls.txt", sep="", collapse=NULL), append = FALSE)

#exit
q()
n
