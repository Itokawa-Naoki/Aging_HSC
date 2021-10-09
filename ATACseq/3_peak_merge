wget 

zcat  *_q0.001_peaks.narrowPeak.gz | sort -k1,1 -k2,2n > Temp.bed 
bedtools merge -i Temp.bed > Marged_peak.bed 

bedtools map -a Marged_peak.bed -b ${example}.bed -o count -c 4 -null 0 > mapped_count_${example}.bed
