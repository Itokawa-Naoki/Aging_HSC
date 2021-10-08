# Download from SRA
# Example file is young HSC rep1(RNA-seq)
example="SRR13192285"
fastq-dump $exapmle	

# get mm10 annotation file
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -zxvf Mus_musculus_UCSC_mm10.tar.gz

# get index for tophat2(bowtie)
wget "ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip"
mkdir indexes
unzip mm10.zip -d ./indexes/

# Mapping to coding region
tophat2 -r 100 -p 1 -G Mus_musculus_UCSC_mm10.gtf indexes "${example}.fastq"

# bam sorting and indexing
samtools
