
import os

# demo of indrop R1 and R2 file
# TGAGCCACATC GAGTGATTGCATGTGACGCCTT TTGCATAT ATTGGT TTTT
# @TGAGCCACATC-TTGCATAT_ATTGGT_NS500663_216_HYM7NBGXX_1_22307_16608_8582
# GCCTGTCGGAGGAGCGCATCAGCGCGCACCACGTG

# GCTTACCT GAGTGATTGCTTGTGACGCCTT ATGCATGG TCAACT TTTTTTT
# @GCTTACCT-ATGCATGG_TCAACT_NS500663_216_HYM7NBGXX_1_22307_18303_8570
# CCCCCTTTAGTCATGGATTTCTATTTGTTTTTTAA

project_root = '/RAID_32T/fbh/at_ze/id_ze_en_4t24_wi_2294'
uniform_q = 'F'
adapter = 'GAGTGATTGCTTGTGACGCCTT'
log_fre = 100000000


def run():
    mis_list = os.listdir(project_root)
    for mis in mis_list:
        mis_path = os.path.join(project_root, mis)
        exFastq_path = os.path.join(mis_path, 'exFastq')
        fastq_path = os.path.join(mis_path, 'fastq')
        exFastq_list = os.listdir(exFastq_path)
        for exFastq in exFastq_list:
            exFastq_file = os.path.join(exFastq_path, exFastq)
            handleExFastq_e(exFastq_file)
        


def handleExFastq(exFastq_file):
    r1 = []
    r2 = []
    r1_file = exFastq_file.replace('.fastq', '_R1.fastq')
    r2_file = exFastq_file.replace('.fastq', '_R2.fastq')
    r1_file = r1_file.replace(r'/exFastq/', r'/fastq/')
    r2_file = r2_file.replace(r'/exFastq/', r'/fastq/')
    exFastq_file_lines = countlines(exFastq_file)
    print('exFastq_file lines: ' + str(exFastq_file_lines))
    with open(exFastq_file, mode='r') as f:
        lineNo = 0
        ll = [None, None, None, None]
        for l in f:
            if lineNo%log_fre == 0:
                print(str(lineNo) + ' finished, total: ' + str(exFastq_file_lines))
            ind = lineNo%4
            ll[ind] = l
            if ind == 0 and not l.startswith('@'):
                print(' '.join(['fastq parse erorr in line:', lineNo]))
            elif ind == 2 and not l.startswith('+'):
                print(' '.join(['fastq parse erorr in line:', lineNo]))
            if ind == 3:
                r1Re = assembleR1(ll)
                r2Re = assembleR2(ll)
                if len(r1Re) == 0 or len(r2Re) == 0:
                    continue
                else:
                    r1 = r1 + r1Re
                    r2 = r2 + r2Re
            lineNo = lineNo + 1
    print('starting writing R1 ... ')
    with open(r1_file, mode='a') as f:
        lo = len(r1)
        ind = 0
        for l in r1:
            if ind%log_fre == 0:
                print(str(ind) + ' finished, total: ' + str(lo))
            f.write(l)
            ind = ind + 1
            if ind != lo:
                f.write('\n')
        # f.write('\n'.join(r1))
    print('starting writing R2 ... ')
    with open(r2_file, mode='a') as f:
        lo = len(r2)
        ind = 0
        for l in r2:
            f.write(l)
            ind = ind + 1
            if ind != lo:
                f.write('\n')
        # f.write('\n'.join(r2))
        
def handleExFastq_e(exFastq_file):
    r1 = []
    r2 = []
    r1_file = exFastq_file.replace('.fastq', '_R1.fastq')
    r2_file = exFastq_file.replace('.fastq', '_R2.fastq')
    r1_file = r1_file.replace(r'/exFastq/', r'/fastq/')
    r2_file = r2_file.replace(r'/exFastq/', r'/fastq/')
    exFastq_file_lines = countlines(exFastq_file)
    print('processing: ' + exFastq_file + ', lines: ' + str(exFastq_file_lines))
    with open(exFastq_file, mode='r') as f1, open(r1_file, mode='a') as f2, open(r2_file, mode='a') as f3:
        lineNo = 0
        for l in f1:
            if lineNo%log_fre == 0:
                print(str(lineNo) + ' finished, total: ' + str(exFastq_file_lines))
            ind = lineNo%4
            if ind == 0:
                try:
                    ll = l.split(' ')
                    c = ll[1]
                    cc = c.split(':')
                    barcode = cc[0]
                    barcode = barcode.replace('-', adapter)
                    line2 = barcode + cc[1]
                    name = '@' + cc[2]
                    f2c = '\n'.join([name, line2, '+', len(line2)*uniform_q])
                    f2.write(f2c)
                    if lineNo != exFastq_file_lines:
                        f2.write('\n')
                    f3.write(name + '\n')
                except:
                    continue
                    print('wrong at line: ' + lineNo)
                    lineNo = lineNo + 4
            elif ind == 2:
                f3.write('+\n')
            else:
                f3.write(l)
            lineNo = lineNo + 1

                #f3.write('\n')
    r1_lines = countlines(r1_file)
    r2_lines = countlines(r2_file)
    log = r1_file.replace('_R1.fastq', '.log')
    with open(log, mode='w') as f:
        if r1_lines == r2_lines and r1_lines == exFastq_file_lines:
            f.write('length consisted')
        else:
            f.write('length not consisted')

def countlines(fp):
    with open(fp, 'r') as fp:
        for count, line in enumerate(fp):
            pass
        return count + 1
        

# @SRR11699721.16 CCAGACAG-CACAAGGC:ATATCT:NS500422_445_H77LWBGX2_1_11101_2327_1607 length=61
# GGCATGGATGAGCTCTACAAATAAGATACTCCAGCAAAACGATCATCATACGTATCCGGAA
# +SRR11699721.16 CCAGACAG-CACAAGGC:ATATCT:NS500422_445_H77LWBGX2_1_11101_2327_1607 length=61
# AAAAAEEEEEEEEEEEEEEEEEEEEEAEEEEAEAAAAEEEAEAEEEEE<EEAAEAAAAEAE
def assembleR1(ll):
    if len(ll) != 4 or len(ll[0].split(' ')) != 3 or len(ll[0].split(' ')[1].split(':')) != 3:
        print(ll)
        return []
    con = ll[0].split(' ')[1].split(':')
    line1 = '@' + con[2]
    barcode = con[0].replace('-', adapter)
    umi = con[1]
    line2 = barcode + umi
    line3 = '+' + con[2]
    line4 = uniform_q*len(line2)
    return [line1, line2, line3, line4]


def assembleR2(ll):
    con = ll[0].split(' ')[1].split(':')
    name = con[2]
    line1 = '@' + name
    line2 = ll[1]
    line3 = '+' + name
    line4 = ll[3]
    return [line1, line2, line3, line4]

run()
