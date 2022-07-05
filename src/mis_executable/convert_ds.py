import os
import subprocess
import shutil

project_root = '/RAID_32T/fbh/at_ze/ds_ze_en_4t12_wi_6474'
uniform_q = 'F'
log_fre = 25000000
log_path = os.path.join(project_root, 'log.txt')
config_path = './convert_ds_con/convert_ds_con.txt'

#fastq_path = os.join.path(project_root, 'fastq')

def bam2sam():
    mis_list = os.listdir(project_root)
    for mis in mis_list:
        mis_path = os.path.join(project_root, mis)
        sam_path = os.path.join(mis_path, 'sam')
        bam_path = os.path.join(mis_path, 'bam')
        if not os.path.exists(bam_path):
            print(bam_path + ' is not exists, skip mis: ', mis)
            continue
        # if os.path.exists(sam_path):
        #     print(sam_path + ' is currently exists, pls check, skip mis: ', mis')
        #     return
        # os.mkdir(sam_path)
        bam_list = os.listdir(bam_path)
        for bam in bam_list:
            bam_file = os.path.join(bam_path, bam)
            sam_name = bam.replace('.bam', '.sam')
            sam_file = os.path.join(sam_path, sam_name)
            cmd = ['samtools', 'view', '-o', sam_file, bam_file]
            print(cmd)
            p=subprocess.Popen(cmd)

# NB501164:249:HGJ55BGXY:4:22601:23347:10702	16	1	2420	1	50M	*	0	0	TCCTAGTAAATTATGAGGTAACAAATCTTCAGTAGATCTCATGTTTGAAT	EEEEEEEAE</EEEEEEEEEE/EE<EEE<EEEEEEEEEEE/EEEEAAAAA	XC:Z:TAGTCTTGAAGA	MD:Z:12C11G25	XF:Z:INTERGENIC	NM:i:2	XM:Z:AACGACTC	UQ:i:63	AS:i:-9	RG:Z:A	PG:Z:bowtie2
# @SRR6890845.1 AGGCCTAAG-CTATAGAG:CATCAT:NS500422_345_H3VKVBGXY_1_11101_17527_1084_ length=33
# GTCTAGCGATTTGGCACGCATGTACCGTGATGA
# +SRR6890845.1 AGGCCTAAG-CTATAGAG:CATCAT:NS500422_345_H3VKVBGXY_1_11101_17527_1084_ length=33
# AAAAAEEEEEEEEEAEEEEEEEEEAAEEEEEEE
def run():
    # mis_list = os.listdir(project_root)
    mis_list = get_mist_list(config_path)
    print(mis_list)
    for mis in mis_list:
        if mis == 'log.txt':
            continue
        mis_path = os.path.join(project_root, mis)
        if not os.path.exists(mis_path):
            print('can not find: ', mis_path, ', skipping this mission.')
            continue
        sam_path = os.path.join(mis_path, 'sam')
        if not os.path.exists(sam_path):
            print('can not find: ', sam_path, ', skipping this mission.')
            continue
        fastq_path = os.path.join(mis_path, 'fastq')
        if os.path.exists(fastq_path):
            shutil.rmtree(fastq_path)
        os.mkdir(fastq_path)
        sam_list = os.listdir(sam_path)
        for sam in sam_list:
            sam_file = os.path.join(sam_path, sam)
            r1_file = sam_file.replace('.sam', '_R1.fastq')
            r2_file = sam_file.replace('.sam', '_R2.fastq')
            r1_file = r1_file.replace(r'/sam/', r'/fastq/')
            r2_file = r2_file.replace(r'/sam/', r'/fastq/')
            to = countlines(sam_file)
            with open(sam_file, mode='r') as f1, open(r1_file, mode='a') as f2, open(r2_file, mode='a') as f3, open(log_path, mode='a') as logf:
                log(logf, 'processing sam file:', sam_file, 'total line:', str(to))
                lineNo = 0
                for l in f1:
                    lineNo = lineNo + 1
                    if lineNo%log_fre == 0:
                        log(logf, str(lineNo), 'finished, total:', str(to))
                    name, sequence, quality, barcode, umi = None, None, None, None, None
                    try:
                        ll = l.split('\t')
                        name = ll[0]
                        sequence = ll[9]
                        quality = ll[10]
                        sequenceL = len(sequence)
                        qualityL = len(quality)
                        if sequenceL != qualityL:
                            log(f, 'read:', name, 'has inconsist length, we justify it.', 'sequence:', sequence, 'length:', length)
                            if qualityL < sequenceL:
                                quality = quality + uniform_q * (sequenceL - qualityL)
                            else:
                                quality = quality[0:sequenceL]
                        if len(ll) > 16 and ll[11].startswith('XC:') and ll[16].startswith('XM:'):
                            barcode = ll[11].split(':')[2]
                            umi = ll[16].split(':')[2].strip()
                        else:
                            i = 11
                            for i in range(0, len(ll)):
                                tag = ll[i]
                                if tag.startswith('XC:'):
                                    barcode = tag.split(':')[2]
                                if tag.startswith('XM:'):
                                    umi = tag.split(':')[2].strip()
                                i = i + 1
                        if barcode == None or umi == None:
                            log(logf, 'read', name, 'has very wrong formation, skipping this read, pls check.')
                            continue
                        if len(barcode) != 12 or len(umi) != 8:
                            log(logf, 'read', name, 'has volitated length, skipping this read.', 'barcode:', barcode, 'umi:', umi)
                            continue
                    except Exception as err:
                        log(logf, 'read:', name, 'parsing wrong, pls check!')
                        print(err)
                        continue
                    else:
                        line1 = '@' + name + '\n'
                        line3 = '+\n'
                        f2.write(line1)
                        f2.write(barcode + umi + '\n')
                        f2.write(line3)
                        f3.write(line1)
                        f3.write(sequence + '\n')
                        f3.write(line3)
                        f2.write(uniform_q*len(barcode + umi) + '\n')
                        f3.write(quality + '\n')
                        #if lineNo != to:
                        #    f2.write(uniform_q*len(barcode + umi) + '\n')
                        #    f3.write(quality + '\n')
                        #else:
                        #    f2.write(uniform_q*len(barcode + umi))
                        #    f3.write(quality)
        
        
def log(f, *msg):
    o = ' '.join(msg) + '\n'
    print(o)
    f.write(o)
    
def countlines(fp):
    with open(fp, 'r') as fp:
        for count, line in enumerate(fp):
            pass
        return count + 1
        
def get_mist_list(p):
    mis_list = []
    with open(p, mode='r') as f:
        for l in f:
            mis_list.append(l.strip())
    return mis_list
        

def gzip_mis():
    mis_list = get_mist_list(config_path)
    for mis in mis_list:
        mis_root = os.path.join(project_root, mis)
        fastq_path = os.path.join(mis_root, 'fastq')
        cmd = ['gzip', '-r', fastq_path]
        print(cmd)
        p=subprocess.Popen(cmd)

gzip_mis()
#run()