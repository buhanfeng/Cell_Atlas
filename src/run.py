
import os
import constant
import enConvert
import subprocess
import shlex
import shutil

##
# x1_ze_cht_52_wi_120581
# x1_ze_en_12_wi_141664
# x1_ze_en_12_wi_157816
# x1_ze_en_24_wi_564810_12723307
# x1_ze_en_24_wi_564810_12723308
# x1_ze_en_30_wi_141664
# x1_ze_en_30_wi_152982
# x1_ze_en_48_wi_564810_12723309
# x1_ze_en_48_wi_564810_12723310
# x1_ze_en_8_wi_141664
# x1_ze_hindbrain_16_wi_141428
# x1_ze_hindbrain_24_wi_141428
# x1_ze_hindbrain_44_wi_141428
# x1_ze_Wretina_14_wi_122680_10440903
# x1_ze_Wretina_14_wi_122680_10440904
# x1_ze_Wretina_24_wi_122680
# x1_ze_Wretina_36_wi_122680
# x1_ze_Wretina_48_wi_122680_10440905
# x1_ze_Wretina_48_wi_122680_13683432
# x1_ze_Wretina_72_wi_122680_10440902
# x1_ze_Wretina_72_wi_122680_10440908

# C2 CR/Ut/SSo/zU
CELL_RANGER = 'CR'
UMI_TOOLS = 'Ut'
STAR_SOLO = 'SSo'
Z_UMI = 'zU'

# CPU
CPU = 30
# project root
# PROJECT_ROOT = '/RAID_32T/fbh/at_ze'
PROJECT_ROOT = r'C:\Users\fbh\Desktop\temp\jobs'
# pattern for Ut
BC_PATTERN_10X = 'CCCCCCCCCCCCCCCCNNNNNNNNNN'
# minimum_length for cutadapt
MINIMUM_LENGTH = 20
# ref genome root
REF_GENOME_ROOT = os.path.join(PROJECT_ROOT, 'ref_genome')
# STAR root
STAR_EXE = '/RAID_32T/fbh/tools/STAR-master/source/STAR'
STAR_10X_EXE = ''
# C1
C1_LIST = [
    {'name': 'x1_ze_en_8_wi_141664', 'C2_LIST': [
        {
            'name': 'CR', 'C3_LIST': ['genomeXX'],
            'name': 'SSo', 'C3_LIST': ['genomeXX'],
            'name': 'Ut', 'C3_LIST': ['genomeXX']
        }
    ]}
]
# UMI-tools step
RUN_UMI_TOOLS_PRE = 2


def run_cell_ranger(id, cpu, fastqs, sample, refg):
    command = 'cellranger count --id=$id \
                --localcores=$cpu \
                --fastqs=$fastqs \
                --sample=$sample \
                --transcriptome=$refg'
                
    command = command.replace('$id', id).replace('$cpu', cpu).replace('$fastqs', fastqs).replace('$sample', sample).replace('$refg', refg)
    log('start running CR ...\n', 'command line: ', command)
    # os.system(command)

def run_STARsolo(refg, cpu, readFilesIn, soloType):
    command = '$STAR_EXE \
                --genomeDir $refg \
                --outSAMmultNmax -1 \
                --runThreadN $cpu \
                --readNameSeparator space \
                --outSAMunmapped Within KeepPairs \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesCommand zcat \
                --readFilesIn $readFilesIn \
                --soloType $soloType \
                --soloCBwhitelist None \
                --soloBarcodeReadLength 0 \
                --outReadsUnmapped Fastx'
    command = command.replace('$refg', refg).replace('$cpu', cpu).replace('$readFilesIn', readFilesIn).replace('$soloType', soloType)
    log_command_line(command)
    # os.system(command)

def run_UMI_tools_after_trim(r1_name, r2_name, c2_path, c3_path, cpu, minimum_length, refg, gtf):
    r1_extract = r1_name.replace('.fastq.gz', '') + '_extracted' + '.fastq.gz'
    r2_extract = r2_name.replace('.fastq.gz', '') + '_extracted' + '.fastq.gz'
    r1_trim = r1_extract.replace('_extracted', '_trimmed')
    r2_trim = r2_extract.replace('_extracted', '_trimmed')
    r1_extract = os.path.join(c2_path, 'extract', r1_extract)
    r2_extract = os.path.join(c2_path, 'extract', r2_extract)
    r1_trim = os.path.join(c3_path, 'extract', r1_extract)
    r2_trim = os.path.join(c3_path, 'extract', r2_extract)
    command3 = 'cutadapt \
    --cores=$cpu \
    --minimum-length $minimum_length \
    -a "CCCATGTACTCTGCGTTGATACCACTGCTTX;min_overlap=5" \
    -A "A{100}X;min_overlap=5" \
    -g "XT{100};min_overlap=10" \
    -G "XAAGCAGTGGTCTCAACGCAGAGTACATGGG;min_overlap=15" \
    -o $r1_trim \
    -p $r2_trim $r1_extract $r2_extract'
    command3 = command3.replace('$cpu', cpu).replace('$minimum_length', minimum_length).replace('$r1_trim', r1_trim)\
        .replace('$r2_trim', r2_trim).replace('$r1_extract', r1_extract).replace('$r2_extract', r2_extract)
    log_command_line(command3)
    # os.system(command3)
    command4 = 'STAR \
                --genomeDir $refg \
                --outSAMmultNmax -1 \
                --runThreadN $cpu \
                --readNameSeparator space \
                --outSAMunmapped Within KeepPairs \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesCommand zcat \
                --readFilesIn $r2_trim'
    command4 = command4.replace('$refg', refg).replace('$cpu', cpu).replace('$r2_trim', r2_trim)
    log_command_line(command4)
    # os.system(command4)
    bam = os.path.join(c3_path, 'STAR', 'Aligned.sortedByCoord.out.bam')
    command5 = 'featureCounts \
                -a $gtf \
                -o gene_assigned \
                -R BAM $bam \
                -T $cpu'
    command5 = command5.replace('$gtf', gtf).replace('$bam', bam).replace('$cpu', cpu)
    log_command_line(command5)
    # os.system(command5)
    bam1 = os.path.join(c3_path, 'featureCounts', 'Aligned.sortedByCoord.out.bam.featureCounts.bam')
    bam2 = os.path.join(c3_path, 'featureCounts', 'assigned_sorted.bam')
    command6 = 'samtools sort $bam1 -o $bam2'
    command6 = command6.replace('$bam1', bam1).replace('$bam2', bam2)
    log_command_line(command6)
    # os.system(command6)
    bam = os.path.join(c3_path, 'featureCounts', 'Aligned.sortedByCoord.out.bam.featureCounts.bam')
    command7 = 'samtools index $bam2'
    command7 = command7.replace('$bam2', bam2)
    log_command_line(command7)
    # os.system(command7)
    command8 = 'umi_tools count \
    --per-gene \
    --gene-tag=XT \
    --wide-format-cell-counts \
    --assigned-status-tag=XS \
    --per-cell \
    -I $bam2 \
    -S $matrix'
    matrix = os.path.join(c3_path, 'count', 'counts.tsv.gz')
    command8 = command8.replace('$bam2', bam2).replace('$matrix', matrix)
    log_command_line(command8)
    # os.system(command8)
    command9 = 'gunzip $matrix'.replace('$matrix', matrix)
    log_command_line(command9)
    # os.system(command9)
    
def run_UMI_tools_beforeTrim(r1_name, r2_name, bc_pattern, c1_path, c2_path):
    r1_extract = r1_name.replace('.fastq.gz', '') + '_extracted' + '.fastq.gz'
    r2_extract = r2_name.replace('.fastq.gz', '') + '_extracted' + '.fastq.gz'
    r1_extract = os.path.join(c2_path, 'extract', r1_extract)
    r2_extract = os.path.join(c2_path, 'extract', r2_extract)
    r1_name = os.path.join(c1_path, 'fastq', r1_name)
    r2_name = os.path.join(c1_path, 'fastq', r2_name)
    wp = os.path.join(c2_path, 'whitelist', 'whitelist.txt')
    command1 = 'umi_tools whitelist \
                --stdin $r1_name \
                --bc-pattern=$bc_pattern \
                --log2stderr > $wp'
    command1 = command1.replace('$r1_name', r1_name).replace('$bc_pattern', bc_pattern).replace('$wp', wp)
    log_command_line(command1)
    # os.system(command1)
    command2 = 'umi_tools extract \
                --bc-pattern=$bc_pattern \
                --stdin $r1_name \
                --stdout $r1_extract \
                --read2-in $r2_name \
                --read2-out=$r2_extract \
                --whitelist=$wp'
    command2 = command2.replace('$bc_pattern', bc_pattern).replace('$r1_name', r1_name)\
        .replace('$r1_extract', r1_extract).replace('$r2_name', r2_name).replace('$r2_extract', r2_extract).replace('$wp', wp)
    log_command_line(command2)
    args = shlex.split(command2)
    print(args)
    ## os.system(command2)
    # p=subprocess.Popen(args)
    

def run_zUMI():
    pass
    
def log(*a):
    mes = ''
    for e in a:
        mes = mes + ' ' + e
    print(mes)
    
def log_command_line(cl):
    print('Start execusing command line: ' + cl)
    
def set_refg(name):
    refg = os.path.join(REF_GENOME_ROOT, name)
    if not os.path.exists(refg):
        log(refg, 'is not exist')
        return 'none'
    else:
        return refg
        
def set_gtf(name):
    name = name.replace('10x_Index_' ,'').replace('STAR_Index_', '')
    gtf = os.path.join(REF_GENOME_ROOT, name)
    gtf = os.path.join(gtf, 'genes.gtf')
    if not os.path.exists(gtf):
        log(gtf, 'is not exist')
        return 'none'
    else:
        return gtf
        
def set_fastq_filename(i, c1_name):
    return c1_name + '_S1_L001_R' + str(i) + '_001.fastq.gz'
    

def set_UMI_evn_after(c3_path):
    featureCounts = os.path.join(c3_path, 'featureCounts')
    STAR = os.path.join(c3_path, 'STAR')
    TPT = os.path.join(c3_path, 'TPT')
    count = os.path.join(c3_path, 'count')
    dl = [featureCounts, STAR, TPT, count]
    for d in dl:
        if os.path.exists(d):
            shutil.rmtree(d) 
        else:
            os.mkdir(d)
    
def set_UMI_evn_before(c2_path):
    whitelist = os.path.join(c2_path, 'whitelist')
    extract = os.path.join(c2_path, 'extract')
    dl = [whitelist, extract]
    for d in dl:
        if os.path.exists(d):
            shutil.rmtree(d) 
        else:
            os.mkdir(d)
    
def change_dir(d):
    log('current pwd: ')
    # os.system('pwd')
    log('try to change to dir:', d)
    if not os.path.exist(d):
        log('dir:', d, 'is not exist!')
    else:
        # os.system('cd ' + d)
        log('change to dir :' + d + 'success!')
    
    
def run():
    cpu = str(CPU)
    for c1 in C1_LIST:
        c1_name = c1['name']
        c1_path = os.path.join(PROJECT_ROOT, c1_name)
        fastq_path = os.path.join(c1_path, 'fastq')
        if not (os.path.exists(c1_path)):
            print(c1_path + ' is not exist, skipping this C1 mission!')
            continue
        if not (os.path.exists(fastq_path)):
            print(fastq_path + ' is not exist, skipping this C1 mission!')
            continue
        else:
            r1_name = set_fastq_filename(1, c1_name)
            r2_name = set_fastq_filename(2, c1_name)
        for c2 in c1['C2_LIST']:
            c2_name = c2['name']
            c2_path = os.path.join(c1_path, c2_name)
            if not (os.path.exists(c2_path)):
                print(c2_path + 'is not exist, creating path!')
                os.mkdir(c2_path)
            else:
                print(c2_path + ' is currently exist, we working under this path!')
            for c3 in c2['C3_LIST']:
                c3_name = c3
                c3_path = os.path.join(c2_path, c3_name)
                refg = set_refg(c3_name)
                gtf = set_gtf(c3_name)
                if refg == 'none':
                    log('no reference genome name: ', c3_name, 'under reference genome root: ', REF_GENOME_ROOT, '\nskip this C3 mission!' )
                    continue
                if os.path.exists(c3_path):
                    print(c3_path + ' is currently exist, remove old file!')
                    shutil.rmtree(c3_path)
                    os.mkdir(c3_path)
                else:
                    print(c3_path + ' is not exist, creating path!')
                    os.mkdir(c3_path)
                if c2_name == CELL_RANGER:
                    id = c3_name
                    shutil.rmtree(c3_path)
                    fastqs = fastq_path
                    sample = c1_name
                    change_dir(c2_path)
                    log('start running cell ranger...')
                    run_cell_ranger(id, cpu, fastqs, sample, refg)
                elif c2_name == STAR_SOLO:
                    readFilesIn = r2_name + ' ' + r1_name
                    soloType = 'CB_UMI_Simple'
                    log('start running STARsolo...')
                    run_STARsolo(refg, cpu, readFilesIn, soloType)
                elif c2_name == UMI_TOOLS:
                    set_UMI_evn_after(c3_path)
                    bc_pattern = BC_PATTERN_10X
                    minimum_length = str(MINIMUM_LENGTH)
                    log('start running UMI_tools...')
                    run_UMI_tools_after_trim(r1_name, r2_name, c2_path, c3_path, cpu, minimum_length, refg, gtf)
                    matrix = os.path.join(c3_path, 'count', 'counts.tsv.gz')
                    enConvert.en_geneId2geneName(matrix)
                elif c2_name == Z_UMI:
                    run_zUMI()


def ut_in_C2_LIST(C2_LIST):
    ll = []
    for e in C2_LIST:
        ll.append(e['name'])
    return (UMI_TOOLS in ll)

def run_UMI_tools_pre():
    for c1 in C1_LIST:
        c1_name = c1['name']
        b = ut_in_C2_LIST(c1['C2_LIST'])
        if not b:
            log(c1_name + UMI_TOOLS + ' is not assigned, skipping this C1 mission!')
            continue
        c1_path = os.path.join(PROJECT_ROOT, c1_name)
        fastq_path = os.path.join(c1_path, 'fastq')
        if not (os.path.exists(c1_path)):
            print(c1_path + ' is not exist, skipping this C1 mission!')
            continue
        if not (os.path.exists(fastq_path)):
            print(fastq_path + ' is not exist, skipping this C1 mission!')
            continue
        else:
            r1_name = set_fastq_filename(1, c1_name)
            r2_name = set_fastq_filename(2, c1_name)
        c2_name = UMI_TOOLS
        c2_path = os.path.join(c1_path, c2_name)
        if not (os.path.exists(c2_path)):
            print(c2_path + ' is not exist, creating path!')
            os.mkdir(c2_path)
        else:
            print(c2_path + ' is currently exist, we working under this path!')
        set_UMI_evn_before(c2_path)
        bc_pattern = BC_PATTERN_10X
        run_UMI_tools_beforeTrim(r1_name, r2_name, bc_pattern, c1_path, c2_path)


if RUN_UMI_TOOLS_PRE == 1:
    run_UMI_tools_pre()
else:
    run()
	
	
	
	
