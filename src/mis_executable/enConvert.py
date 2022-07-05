
import os
import pandas as pd

def parse_gtf(gtf):
    d = {}
    with open(gtf, mode='r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            ll = line.split('\t')
            ann = ll[8]
            annl = ann.split('; ')
            geneId = annl[0]
            geneName = annl[2]
            geneId = geneId.replace('gene_id "', '').replace('"', '')
            geneName = geneName.replace('gene_name "', '').replace('"', '')
            d[geneId] = geneName
    id2name = os.path.join('/RAID_32T/fbh/at_ze/jobs', 'id2name.tsv')
    with open(id2name, mode='w') as f:
        f.write('gene_id' + '\t' 'gene_name' + '\n')
        for e in d:
            f.write(e + '\t' + d[e] + '\n')
            
def en_geneId2geneName(matrix):
    id2name = os.path.join('/RAID_32T/fbh/at_ze/jobs', 'id2name.tsv')
    d = {}
    if not os.path.exists(id2name):
        print(id2name + ' is not exist, pls run parse_gtf, halt!')
        return False
    else:
        with open(id2name, mode='r') as f:
            for line in f:
                geneId = line.split('\t')[0]
                geneId = geneId.strip()
                geneName = line.split('\t')[1]
                geneName = geneName.strip()
                d[geneId] = geneName
    df = pd.read_csv(matrix)
    df = df.fillna(0)
    df = df.rename(index = d)
    p = os.path.realpath(matrix)
    pp = os.path.join(p, 'counts_name.tsv')
    df.to_csv(pp, sep = '\t')