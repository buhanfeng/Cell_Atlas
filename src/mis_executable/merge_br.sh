

for line in `cat ./rename_cnc/rename_cnc.txt`
do
	line=$line |tr -d '\r\n'
	line=$line |tr -d '\n'
    echo 'precessing mis: '$line
    fastq='/RAID_32T/fbh/at_ze/x_ze_cnc_36t5040_wi_8969/'$line'/fastq/'
    cd $fastq
    cat SRR????_1.fastq.gz > ${line}_S1_L001_R1_001.fastq.gz
    cat SRR????_2.fastq.gz > ${line}_S1_L001_R2_001.fastq.gz
done