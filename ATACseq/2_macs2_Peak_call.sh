example="SRR13186609"
ratio=".6378"

bedtools bamtobed -i ${example}_d${ratio}.bam > ${example}_d${ratio}.bed
macs2 callpeak -t ${example}_d${ratio}.bed -n ${example} -q 0.001 -f BED -g mm --keep-dup all --SPMR --call-summits
