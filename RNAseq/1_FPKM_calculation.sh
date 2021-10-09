# Download from SRA
# Example file is young HSC rep1(RNA-seq)
example="SRR13192285"
fastq-dump $exapmle	

# get mm10 reference and annotation file
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -zxvf Mus_musculus_UCSC_mm10.tar.gz

# Mapping to coding region
tophat2 -r 100 -p 1 -G Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf \
Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
"${example}.fastq"

# bam sorting and indexing
cd tophat_out
samtools sort -O bam -o sort_accepted_hits.bam accepted_hits.bam 
samtools index sort_accepted_hits.bam
cd ..

# calculation of raw read count for DESeq2 and FPKM by stringtie
stringtie -e -B -p 16 -G Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf -o ${example}.gtf -A ${example}_abd.txt tophat_out/sort_accepted_hits.bam
