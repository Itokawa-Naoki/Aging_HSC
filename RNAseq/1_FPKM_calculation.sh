# Download from SRA
# Example file is young HSC rep1(RNA-seq)
example="SRR13192285"
fastq-dump $exapmle	

# get mm10 annotation file
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -zxvf Mus_musculus_UCSC_mm10.tar.gz
bowtie2-build -f Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa mm10_index

# get index for tophat2(bowtie)
wget "ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip"
mkdir indexes
unzip mm10.zip -d ./indexes/

# Mapping to coding region
tophat2 -r 100 -p 1 -G Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf mm10_index "${example}.fastq"

# bam sorting and indexing
samtools sort -O bam -o sort_accepted_hits.bam accepted_hits.bam 
samtools index sort_accepted_hits.bam

# calculation of FPKM and read count for DESeq2 by stringtie
stringtie -e -B -p 16 -G Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf -o ${example}.gtf -A ${example}_abd.txt sort_accepted_hits.bam
