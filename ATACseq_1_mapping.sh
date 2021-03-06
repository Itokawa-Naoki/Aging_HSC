# Download from SRA
# Example file is young HSC rep1(ATAC-seq)
example="SRR13186609"
fastq-dump $example	

# get mm10 reference file
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -zxvf Mus_musculus_UCSC_mm10.tar.gz

#mapping
bowtie2 -x Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -q "${example}.fastq" -S ${example}.sam

samtools view -bS ${example}.sam > ${example}.bam

#elimination of mitochondria reads
samtools view -hF 4 ${example}.bam | grep -vF chrM | samtools view -bS > ${example}_uM.bam

samtools sort ${example}_uM.bam > sorted_${example}_uM.bam
samtools index sorted_${example}_uM.bam

# subsampling. The number of reads in all samples was subsampled to match the lowest in the sample.
ratio=".6378"
samtools view -h -b -s 660${ratio} sorted_${example}_uM.bam > ${example}_d${ratio}.bam
samtools index ${example}_d${ratio}.bam
