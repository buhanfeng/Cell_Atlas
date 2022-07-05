

project_root='/RAID_32T/fbh/at_ze/x_ze_cnc_36t5040_wi_8969/test_root'
config_file='./rename_cnc/rename_cnc_test.txt'
log_file=${project_root}/SSo_sh_log.txt

run_ds()
{
    /RAID_32T/fbh/tools/STAR-master/source/STAR \
     --genomeDir /RAID_32T/fbh/at_ze/ref_genome/STAR_Index_GRCz11_EN/star \
     --outSAMmultNmax -1 \
     --runThreadN 30 \
     --readNameSeparator space \
     --outSAMunmapped Within KeepPairs \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat \
     --readFilesIn ${2} ${3} \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist None \
     --soloBarcodeReadLength 0 \
     --soloCBstart 1 \
     --soloCBlen 12 \
     --soloUMIstart 13 \
     --soloUMIlen 8 \
     --outFileNamePrefix ${1}
}

run_v2()
{
    /RAID_32T/fbh/tools/STAR-master/source/STAR \
     --genomeDir /RAID_32T/fbh/at_ze/ref_genome/STAR_Index_GRCz11_EN/star \
     --outSAMmultNmax -1 \
     --runThreadN 30 \
     --readNameSeparator space \
     --outSAMunmapped Within KeepPairs \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat \
     --readFilesIn ${2} ${3} \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist None \
     --soloBarcodeReadLength 0 \
     --outFileNamePrefix ${1}
}

for line in `cat ${config_file}`
do
	line=$line |tr -d '\r\n'
	line=$line |tr -d '\n'
    echo 'precessing mis: '$line
    `echo precessing mis: ${line} >> ${log_file}`
    fastq=${project_root}/${line}/fastq
    R2=${fastq}/${line}_S1_L001_R2_001.fastq.gz
    R1=${fastq}/${line}_S1_L001_R1_001.fastq.gz
    out=${project_root}/${line}/SSo/
    if [ ! -d ${out} ];then
        echo ${out}' is not exist, make dir '${out}
        `echo ${out} is not exist, make dir ${out} >> ${log_file}`
        mkdir ${out}
    else
            echo 'removing and creating '${out}
            rm -R ${out}
            mkdir ${out}
    fi
    #run_ds $out $R2 $R1
    run_v2 $out $R2 $R1
    gzip ${out}Solo.out/Gene/filtered/*
    gzip ${out}Solo.out/Gene/raw/*
done



