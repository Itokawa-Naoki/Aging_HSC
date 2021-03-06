library("DESeq2")

# getting data generated by ATAC-seq read count in merged peak file
countData <- read.csv(url("https://github.com/Itokawa-Naoki/Aging_HSC/raw/main/data/ATAC_mapped_count.csv"),header=TRUE)
row.names(countData) <- countData[,1]
countData <- countData[,-1]
colData <- read.csv(url("https://github.com/Itokawa-Naoki/Aging_HSC/raw/main/data/ATAC_phenodata_DESeq2.csv"),row.names = 1)
colData[,1] <- as.factor(colData[,1])

# DESeq2 calculation
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, design = ~ covariate)
dds <- DESeq(dds)

# comparison of each fraction between young and aged
levels <- levels(factor(colData[,1]))
for(i in 1:8) {
  res1 <- results(dds, contrast = c("covariate", levels[i], levels[i+8]))  
  rres1 <- cbind(ID=rownames(res1),res1)
  write.csv(rres1,file = paste0(levels[i],"-vs-",levels[i+8],"_ATAC.csv"), row.names = F)
}


# normalized count by DESeq2
ntd <- normTransform(dds)
antd <- assay(ntd)
antd2 <- 2^antd-1
write.csv(antd,file="ATACseq_DESeq2_norm_count.csv",row.names = TRUE)
