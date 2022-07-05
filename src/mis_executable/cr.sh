



project_root='/RAID_32T/fbh/at_ze/x_ze_br_12t360_wi_8142'
config_file='./sp_mis/br.txt'
log_file=${project_root}/cr_sh_log.txt
ref_genome_root=/RAID_32T/fbh/at_ze/ref_genome
refv='10x_Index_GRCz11_EN'
ref_genome=${ref_genome_root}/${refv}


for line in `cat ${config_file}`
do
	line=$line |tr -d '\r\n'
	line=$line |tr -d '\n'
    echo 'precessing mis: '$line
    `echo precessing mis: ${line} >> ${log_file}`
    fastq=${project_root}/${line}/fastq
    out=${project_root}/${line}/CR/
    if [ ! -d ${out} ];then
        echo ${out}' is not exist, make dir '${out}
        `echo ${out} is not exist, make dir ${out} >> ${log_file}`
        mkdir ${out}
    fi
    id=${refv}
    target=${out}${id}
    if [ ! -d ${target} ];then
        echo ${target}' is not exist, start running CR ... '
        `echo ${target} is not exist, start running CR ...  >> ${log_file}`
    else
            echo 'removing old file under: '${target}
            rm -R ${target}
    fi
    echo ${id}
    echo ${fastq}
    echo ${line}
    echo ${ref_genome}
    cd ${out}
    cellranger count --id=${id} \
     --fastqs=${fastq} \
     --sample=${line} \
     --transcriptome=${ref_genome}
done



