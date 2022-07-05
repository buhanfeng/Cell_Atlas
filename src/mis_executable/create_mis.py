import subprocess
import os

project_root = '/RAID_32T/fbh/at_ze/x_ze_cnc_36t5040_wi_8969'
config_file = './rename_cnc/rename_cnc.txt'
buffer_path = ''


#def set_mis_root(mis_name):
#    print(['setting mis root:', mis_name].join(' '))
#    mis_root = os.path.join(project_root, mis_name)
#    if os.path.exists(mis_root):
#        print([mis_name, 'already exists, skip this mis'].join(' '))
#    else:
#        os.mkdir(mis_root)
#    
#def set_srr(mis_name, sl):
#    srr_path = os.path.join(project_root, mis_name, 'srr')
#    os.mkdir()
#    sl = sl.split(' ')
#    for s in sl:
#        cmd = ['mv', os.path.join(buffer_path, s), srr_path].join(' ')
#        os.system(cmd)


# maximum of 40 threads
def set_fastq():
    with open(config_file, mode='r') as f:
        for l in f:
            ll = l.split('\t')
            print(' '.join(['set_fastq for:', l]))
            mis_name = ll[0]
            sl = ll[1].split('/')
            mis_root = os.path.join(project_root, mis_name)
            fastq_path = os.path.join(mis_root, 'fastq')
            srr_path = os.path.join(mis_root, 'srr')
            # os.mkdir(fastq_path)
            for s in sl:
                cmd = []
                print(ll[2])
                s_path = os.path.join(srr_path, s)
                if ll[2] == 'single':
                    cmd = ['fastq-dump', '--gzip', '--origfmt', '--outdir', fastq_path, s_path]
                elif ll[2].strip() == 'pair':
                    print('jere')
                    cmd = ['fastq-dump', '--split-files', '--gzip', '--origfmt', '--outdir', fastq_path, s_path]
                print(cmd)
                p=subprocess.Popen(cmd)


def rename_ds():
    mis_list = get_mist_list(config_file)
    for mis in mis_list:
        mis_path = os.path.join(project_root, mis)
        fastq_path = os.path.join(mis_path, 'fastq')
        fastq_list = os.listdir(fastq_path)
        for fastq in fastq_list:
            fastq_file = os.path.join(fastq_path, fastq)
            new = ''
            if '_R1' in fastq:
                new = mis + '_S1_L001_R1_001.fastq.gz'
            else:
                new = mis + '_S1_L001_R2_001.fastq.gz'
            new = os.path.join(fastq_path, new)
            cmd = 'mv ' + fastq_file + ' ' + new
            print(cmd)
            os.system(cmd)
            
def rename_br():
    mis_list = get_mist_list(config_file)
    for mis in mis_list:
        mis_path = os.path.join(project_root, mis)
        fastq_path = os.path.join(mis_path, 'fastq')
        fastq_list = os.listdir(fastq_path)
        for fastq in fastq_list:
            ll = fastq.split('_')
            prefix = ll[0]
            suffix = ll[1]
            length = len(prefix)
            new = 'SRR' + prefix[length-4: length] + '_' + suffix
            fastq_file = os.path.join(fastq_path, fastq)
            new_file = os.path.join(fastq_path, new)
            cmd = 'mv ' + fastq_file + ' ' + new_file
            print(cmd)
            os.system(cmd)
        srr_path = os.path.join(mis_path, 'srr')
        srr_list = os.listdir(srr_path)
        for srr in srr_list:
            print(srr[length-4: length])
            new = 'SRR' + srr[length-4: length]
            srr_file = os.path.join(srr_path, srr)
            new_file = os.path.join(srr_path, new)
            cmd = 'mv ' + srr_file + ' ' + new_file
            print(cmd)
            os.system(cmd)
            

def get_mist_list(p):
    mis_list = []
    with open(p, mode='r') as f:
        for l in f:
            mis_list.append(l.strip())
    return mis_list

rename_br()
#rename_ds()
#set_fastq()