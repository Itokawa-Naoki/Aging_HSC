wget https://github.com/Itokawa-Naoki/Aging_HSC/raw/main/data/Y_peaks.tar.gz
wget https://github.com/Itokawa-Naoki/Aging_HSC/raw/main/data/A_peaks.tar.gz
tar -zxvf Y_peaks.tar.gz
tar -zxvf A_peaks.tar.gz

cat  *.narrowPeak | sort -k1,1 -k2,2n > Temp.bed 
bedtools merge -i Temp.bed > Marged_peak.bed 
