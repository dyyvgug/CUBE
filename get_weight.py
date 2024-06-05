#!/usr/bin/python
# coding:utf-8
# Yingying Dong. To calculate the mCAI weight.

import re
import os
try:
    import rpy2.robjects as robjects
except:
    os.system('pip install rpy2')
    import rpy2.robjects as robjects


pat1 = re.compile(r'\s+')

aa_codon = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC'],
    'K': ['AAA', 'AAG'], 'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M': ['ATG'],
    'N': ['AAT', 'AAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'Y': ['TAT', 'TAC'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'W': ['TGG'],
    'STOP': ['TAG', 'TAA', 'TGA']
}
codon_count = {
    'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0,
    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0,
    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
    'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
    'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0, 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,
    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0,
    'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0,
    'TGT': 0, 'TGC': 0, 'TGG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 'TGA': 0
}
aa_num = {
    'F': 2, 'Y': 2, 'C': 2, 'H': 2, 'Q': 2, 'N': 2, 'K': 2, 'D': 2, 'E': 2,
    'P': 4, 'T': 4, 'V': 4, 'A': 4, 'G': 4,
    'L': 6, 'R': 6, 'S': 6,
    'W': 1, 'M': 1, 'STOP': 1,
    'I': 3,
}


def rev_com(seq):
    intab = "ATCGatcg"
    outab = "TAGCTAGC"
    trantab = str.maketrans(intab, outab)
    result = seq.translate(trantab)
    return result[::-1]


def CDS_info(line):
    line = pat1.split(line)
    info = [line[0], line[3], line[4], line[6], line[7]]
    return info


def CDS_info2(line2):
    line2 = pat1.split(line2)
    info2 = [line2[0], line2[4], line2[5], line2[7], line2[8]]
    return info2


def ext_CDS(fna, gff):
    CDS = ''
    CDS_dict = {}
    sequ = ''
    for line in fna:
        if line.startswith('>') and sequ == '':
            keys = pat1.split(line)[0].replace('>', '')
            CDS_dict[keys] = []
        elif not line.startswith('>'):
            sequ = sequ + line.strip()
        elif line.startswith('>') and sequ != '':
            CDS_dict[keys] = sequ.upper()
            sequ = ''
            keys = pat1.split(line)[0].replace('>', '')
    fna.close()

    for line in gff:
        line = line.strip()
        if 'ribosomal' in line or 'Ribosomal' in line:
            if 'Mitochondrial' not in line and 'mitochondrial' not in line:
                if 'kinase' not in line and 'ubiquitin' not in line:
                    if 'apicoplast' not in line:
                        if 'CDS' in line:
                            try:
                                i = CDS_info(line)
                                sequence = CDS_dict[i[0]][(int(i[1]) - 1):int(i[2])]
                                sequence = sequence[int(i[4]):]
                                if i[3] == '-':
                                    sequence = rev_com(sequence)
                                CDS += '>' + i[0] + '\n' + sequence + '\n'
                            except:
                                try:
                                    i = CDS_info2(line)
                                    sequence = CDS_dict[i[0]][(int(i[1]) - 1):int(i[2])]
                                    sequence = sequence[int(i[4]):]
                                    if i[3] == '-':
                                        sequence = rev_com(sequence)
                                    CDS += '>' + i[0] + '\n' + sequence + '\n'
                                except:
                                    CDS += ''
                            else:
                                with open('err_file', 'w') as e:
                                    e.write(line + ' is wrong format')

    gff.close()
    return CDS


def calc_freq(codon_count, rscu):
    count_tot = {}
    for aa in aa_codon.keys():
        n = 0
        for codon in aa_codon[aa]:
            n = n + codon_count[codon]
        count_tot[aa] = float(n)
    for aa in aa_codon.keys():
        for codon in aa_codon[aa]:
            if count_tot[aa] != 0.0:
                freq = codon_count[codon] / count_tot[aa]
                RSCU = codon_count[codon] / count_tot[aa] * aa_num[aa]
            else:
                freq = 0.0
                RSCU = 0.0
            rscu += '{}\t{}\t{}\t{:.3f}\t{:.3f}\n'.format(aa, codon, codon_count[codon], freq, RSCU)
    return rscu


def cal_weight(cds):
    in_file = cds.split('\n')
    rs = ''
    dna = ''
    for j in in_file:
        if not j.startswith('>'):
            dna = dna + j.strip()

    rs += 'AA\tcodon\thits\tfrequency\tRSCU\n'

    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        if codon in codon_count:
            codon_count[codon] = codon_count[codon] + 1
    result = calc_freq(codon_count=codon_count, rscu=rs)
    return result


def read_file(fna, gff):
    fa = open('{}'.format(fna), 'r')
    gf = open('{}'.format(gff), 'r')
    cds = ext_CDS(fna=fa, gff=gf)
    rsc = cal_weight(cds)
    with open('rscu.txt', 'w') as f:
        f.write(rsc)

    r_script = '''
    gene_fre = read.table("rscu.txt",header = T,sep = '\t',quote = "")
    df <- gene_fre
    df$Weights <- ave(df$RSCU,df$AA,FUN=function(x) x/max(x))
    df = df[,-c(1,3,4,5)]
    df = df[-c(30,61,62,63,64),]
    write.table(df,file = "weight.txt",sep = '\t',quote = F,row.names = F,col.names = F) # for calculate CAI
    '''
    robjects.r(r_script)

    with open('weight.txt', 'r') as f2:
        lines = f2.readlines()

    f.close()
    f2.close()

    os.remove("rscu.txt")
    os.remove("weight.txt")

    return lines
