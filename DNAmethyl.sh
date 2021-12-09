args <- commandArgs(TRUE)

if (length(args) == 2) {
    for(i in 1:length(args)){
        eval(parse(text = args))
    }
} else {
    stop()
}

library(methylKit)
library(GenomicRanges)
library(openxlsx)

#--format adjust
data <- read.table("Aged_HSC.cpg.methexport.tab")
data <- read.table("Aged_test.txt")#data2 <- data
col1 <- paste(data2[,1],".",data2[,2])
strand <- data2[,3]
freqC <- data2[,4]/data2[,5]*100
freqT <- 100 - freqC
Aged <- data.frame(chrBase=col1,chr=data2[,1],base=data2[,2],strand=strand,coverage=data2[,5],freqC=freqC,freqT=freqT)

#--cont
data <- read.table("Young_Y_HSC.cpg.methexport.tab")
data2 <- data
col1 <- paste(data2[,1],".",data2[,2])
strand <- data2[,3]
freqC <- data2[,4]/data2[,5]*100
freqT <- 100 - freqC
Young <- data.frame(chrBase=col1,chr=data2[,1],base=data2[,2],strand=strand,coverage=data2[,5],freqC=freqC,freqT=freqT)

write.table(Aged,"Aged", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "0",  row.names = FALSE,col.names = TRUE)
write.table(Young,"Young", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "0",  row.names = FALSE,col.names = TRUE)

file.list <- list("Aged","Young")
myobj=methRead(file.list,
               sample.id=list("age","young"),
               assembly="mm10",
               treatment=c(1,0),
               context="CpG",
               mincov = 10)

meth=unite(myobj, destrand=FALSE)

bed <- read.table(bedfile)


region <- GRanges(seqnames=bed[,1], ranges=IRanges(start=bed[,2], end=bed[,3]),strand=bed[,4])
result  <- regionCounts(myobj,region)
write.table(result[[1]],file=paste0("result_Aged_",bedfile),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "0",  row.names = FALSE,col.names = FALSE)
write.table(result[[2]],file=paste0("result_Young_",bedfile),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "0",  row.names = FALSE,col.names = FALSE)
