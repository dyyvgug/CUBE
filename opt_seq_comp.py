#!/usr/bin/python
# coding:utf-8

# ===========================================================================================================
# Author: Yingying Dong.Email: dyyvgug@163.com .
#   Change the sequence codons to optimal, and avoid polyN signal and restriction sites.
# ===========================================================================================================

import os
import sys
import re
import decimal
import argparse
import platform
import get_CDS_rscu as wei

parser = argparse.ArgumentParser(description='Change DNA sequence to improve expression.v1.0', prog='change_opt',
                                 usage='%(prog)s [options]')
parser.add_argument('-inp', nargs='?', type=argparse.FileType('r'), required=True,
                    help='(Required Parameters) The file name of the original sequence.\n\tThe sequence default type is DNA sequence, if that is protein sequence, please add \'-Pro\' parameter')
parser.add_argument('-DNA', action='store_true',
                    help='FASTA file for DNA sequences of genes that wish to increase expression')
parser.add_argument('-Pro', action='store_true',
                    help='FASTA file for protein sequences of genes that wish to increase expression')
parser.add_argument('-spe', nargs='?', type=str, help='Latin name of host species, separated by an underscore, for example: Caenorhabditis_elegans')
parser.add_argument('-genome', nargs='?', type=str, help='The FASTA file of the species genome')
parser.add_argument('-gff', nargs='?', type=str, help='The annotation file GFF3 format of the species')
parser.add_argument('-o', nargs='?', type=str, default='optimized_seq.fa',
                    help='The file name of output optimized sequence.The default file name is \'optimized_seq.fa\'')
parser.add_argument('-poly', action='store_true',
                    help='Need to remove polyN sequence')
parser.add_argument('-res', action='store_true',
                    help='Need to remove restriction enzyme sites')
parser.add_argument('-res_sites1', type=str, default=None,
                    help='Restriction enzyme recognition sequence')
parser.add_argument('-res_sites2', type=str, default=None,
                    help='Restriction enzyme recognition sequence')
parser.add_argument('-res_sites3', type=str, default=None,
                    help='Restriction enzyme recognition sequence')

args = parser.parse_args()


def hete(dataSource, species, sourceType, poly, res, res_sites1, res_sites2, res_sites3):

    if sourceType == 0:

        codon_table = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                       'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                       'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                       'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                       'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                       'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                       'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                       'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                       'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
                       'ATG': 'M', 'GAT': 'D', 'GAC': 'D',
                       'GAA': 'E', 'GAG': 'E', 'CAT': 'H', 'CAC': 'H',
                       'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
                       'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S',
                       'AGA': 'R', 'AGG': 'R', 'TTT': 'F', 'TTC': 'F',
                       'TTA': 'L', 'TTG': 'L', 'TGT': 'C', 'TGC': 'C',
                       'TGG': 'W', 'TAT': 'Y', 'TAC': 'Y',
                       'TAA': 'STOP', 'TAG': 'STOP', 'TGA': 'STOP'}

        dataSource += '\n>'
        f = dataSource.split('\n')

        aa_format = ''
        dna = ''
        for line in f:
            if line.startswith('>') and dna == '':
                header = line.strip()
                aa_format += header + '\n'
            elif not line.startswith('>'):
                dna = str.upper(dna) + line.strip()
            elif line.startswith('>') and dna != '':
                prot = ''
                for i in range(0, len(dna), 3):
                    codon = dna[i:i + 3]
                    if codon in codon_table:
                        if codon_table[codon] == 'STOP':
                            prot = prot + '*'
                        else:
                            prot = prot + codon_table[codon]
                    else:
                        prot = prot + '-'
                j = 0
                while j < len(prot):
                    aa_format += prot[j:j + 48] + '\n'
                    j = j + 48
                dna = ''
                header = line.strip()
                aa_format += header + '\n'
        aa_format += '>'

    if sourceType == 1:

        dataSource += '\n>'
        f2 = dataSource.split('\n')

        aa_format = ''
        AA_seq = ''

        for line in f2:
            if line.startswith('>') and AA_seq == '':
                header = line.strip()
                aa_format += header + '\n'
            elif not line.startswith('>'):
                AA_seq = AA_seq + line.strip()
            elif line.startswith('>') and AA_seq != '':
                AA_seq = AA_seq.upper()
                j = 0
                while j < len(AA_seq):
                    aa_format += AA_seq[j:j + 48] + '\n'
                    j = j + 48
                AA_seq = ''
                header = line.strip()
                aa_format += header + '\n'
        aa_format += '>'

    opt_seq_ori = ''

    aa_codon_fre = {}
    RSCU_table = []
    for i in species:
        RSCU_table.append(i.strip().split('\t'))
    for j in RSCU_table:
        aa_codon_fre[j[1]] = {j[0]: j[4]}
    optimalA = max(decimal.Decimal(x['A']) for x in aa_codon_fre.values() if
                   'A' in x)  # First determine if it exists, otherwise a keyerror will occur.
    keyA = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'A': str(optimalA)})]
    opt_dic = {
        'A': keyA}  # This dictionary stores the optimal codon and corresponding the amino acid.In order to design the amino acid sequence set as the optimal codon DNA sequence.
    optimale = max(decimal.Decimal(x['STOP']) for x in aa_codon_fre.values() if 'STOP' in x)
    keye = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'STOP': str(optimale)})]

    rareA = min(decimal.Decimal(x['A']) for x in aa_codon_fre.values() if 'A' in x)
    raKeyA = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({'A': str(rareA)})]
    sub_dic = {keyA: raKeyA}  # This dictionary stores the optimal codon and rarest codon for an amino acid.

    def findopt(aa, aa_codon_fre, opt_dic, sub_dic):
        optimal = max(decimal.Decimal(x[aa]) for x in aa_codon_fre.values() if aa in x)
        optKey = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({aa: str(optimal)})]
        opt_dic[aa] = optKey
        rare = min(decimal.Decimal(x[aa]) for x in aa_codon_fre.values() if aa in x)
        raKey = list(aa_codon_fre.keys())[list(aa_codon_fre.values()).index({aa: str(rare)})]
        sub_dic[optKey] = raKey
        return opt_dic, sub_dic

    list_aa = ['C', 'D', 'E', 'F', 'G', 'H', 'K', 'I', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'Y', 'T', 'V', 'W']
    for a in list_aa:
        findopt(a, aa_codon_fre, opt_dic, sub_dic)
    opt_dic['*'] = keye

    prot = ''
    aa_format2 = aa_format.split('\n')

    for line in aa_format2:
        if line.startswith('>') and prot == '':
            header = line
            opt_seq_ori += header + '\n'
        elif not line.startswith('>'):
            prot = prot + line.strip()
        elif line.startswith('>') and prot != '':
            opt_dna = ''
            for i in range(0, len(prot), 1):
                j = prot[i]
                if j in opt_dic.keys():
                    opt_dna = opt_dna + opt_dic[j]
                else:
                    print(j + ' It is not amino acid')
            j = 0
            while j < len(opt_dna):
                opt_seq_ori += opt_dna[j:j + 48] + '\n'
                j = j + 48
            prot = ''
            header = line
            opt_seq_ori += header + '\n'
    opt_seq_ori += '>'

    if poly is False and res is False:
        return opt_seq_ori
    else:
        re_rep = ''
        opt_seq = ''
        opt_seq_ori2 = opt_seq_ori.split('\n')
        for line in opt_seq_ori2:
            if line.startswith('>') and opt_seq == '':
                header = line
                re_rep += header + '\n'
            elif not line.startswith('>'):
                opt_seq = opt_seq + line.strip()
            elif line.startswith('>') and opt_seq != '':
                if poly is True and res is True:
                    opt_seq = findrep('C{5,}', opt_seq, sub_dic)  # rep_posC
                    opt_seq = findrep('G{5,}', opt_seq, sub_dic)  # rep_posG
                    opt_seq = findrep('T{5,}', opt_seq, sub_dic)  # rep_posT
                    opt_seq = findrep('A{5,}', opt_seq, sub_dic)  # rep_posA
                    opt_seq = findrep(res_sites1, opt_seq, sub_dic)  # res_site1
                    opt_seq = findrep(res_sites2, opt_seq, sub_dic)  # res_site2
                    opt_seq = findrep(res_sites3, opt_seq, sub_dic)  # res_site3
                elif poly is True and res is False:
                    opt_seq = findrep('C{5,}', opt_seq, sub_dic)  # rep_posC
                    opt_seq = findrep('G{5,}', opt_seq, sub_dic)  # rep_posG
                    opt_seq = findrep('T{5,}', opt_seq, sub_dic)  # rep_posT
                    opt_seq = findrep('A{5,}', opt_seq, sub_dic)  # rep_posA
                elif poly is False and res is True:
                    opt_seq = findrep(res_sites1, opt_seq, sub_dic)  # res_site1
                    opt_seq = findrep(res_sites2, opt_seq, sub_dic)  # res_site2
                    opt_seq = findrep(res_sites3, opt_seq, sub_dic)  # res_site3
                j = 0
                while j < len(opt_seq):
                    re_rep += opt_seq[j:j + 48] + '\n'
                    j = j + 48
                opt_seq = ''
                header = line
                re_rep += header + '\n'
        return re_rep


def findrep(polyN, seq, sub_dic):
    pN = re.compile(str(polyN))
    pos_list = []
    while len(pN.findall(seq)):
        for m in pN.finditer(seq):
            pos = m.start()
            pos_list.append(pos)
        if len(pos_list):
            for i in pos_list:
                if i % 3 == 0:
                    if seq[i:i + 3] in sub_dic.keys():
                        seq = seq[:i] + sub_dic[seq[i:i + 3]] + seq[i + 3:]
                elif i % 3 == 1:
                    if seq[i - 1:i + 2] in sub_dic.keys():
                        seq = seq[:i - 1] + sub_dic[seq[i - 1:i + 2]] + seq[i + 2:]
                elif i % 3 == 2:
                    if seq[i + 1:i + 4] in sub_dic.keys():
                        seq = seq[:i + 1] + sub_dic[seq[i + 1:i + 4]] + seq[i + 4:]
        else:
            print('These sequences no longer have restriction sites or consecutive same nucleotides')
    if len(pN.findall(seq)) == 0:
        print('These sequences no longer have restriction sites or consecutive same nucleotides')
    return seq


if __name__ == '__main__':
    if args.DNA is True and args.Pro is False:
        typ = 0
    elif args.DNA is False and args.Pro is False:
        typ = 0
        print('If no \'-DNA\' or \'-Pro\', the default is \'-DNA\'')
    elif args.DNA is True and args.Pro is True:
        typ = 0
        print('If \'-DNA\' and \'-Pro\'both parameters are added, the default is \'-DNA\'')
    else:
        typ = 1

    if args.spe is None and args.genome is not None and args.gff is not None:
        weight = wei.read_file(args.genome, args.gff)
        weight2 = weight[:-1]
        opt = hete(dataSource=args.inp.read(), species=weight2, sourceType=typ, poly=args.poly, res=args.res,
                   res_sites1=args.res_sites1, res_sites2=args.res_sites2, res_sites3=args.res_sites3)
        with open(args.o, 'w') as fi:
            fi.write(opt)
        fi.close()
    elif args.spe is None and args.genome is None and args.gff is None:
        print('\tThere are no parameters for \'-spe\',\'-genome\', and \'-gff\'.\n\tIf the species you need to calculate is in \'sup_spe.txt\', please follow the name after the \'-spe\' parameter.\n\tIf the calculation of the species is not supported, use \'-genome\' and \'-gff\' followed by genome and annotation files')
    elif args.spe is not None and args.genome is None and args.gff is None:
        syst = platform.system()
        if syst == "Windows":
            os.chdir('.\\')
            r_path = '.\\resource\\RSCU\\'
        elif syst == "Linux":
            os.chdir('./')
            r_path = './resource/RSCU/'

        if os.path.exists('{}{}'.format(r_path, args.spe)):
            RSCU = open('{}{}'.format(r_path, args.spe), 'r')
            opt = hete(dataSource=args.inp.read(), species=RSCU, sourceType=typ, poly=args.poly, res=args.res,
                       res_sites1=args.res_sites1, res_sites2=args.res_sites2, res_sites3=args.res_sites3)
            with open(args.o, 'w') as fi:
                fi.write(opt)
            fi.close()
        else:
            print('\tThe calculation of this species is not supported, and the species that supports calculation are mentioned in the \'supported_species.txt\'.\n\tIf you have the genome and GFF annotation files of the species, you can use \'-gff\' and \'-genome\'.\n\t\'-gff\' and \'-genome\' are followed by the annotation file GFF format of the species and the fasta format file of the genome respectively.')
    else:
        print('\tMissing parameters.\n\tIf you do not have the species name of the \'-spe\' parameter, you must provide the corresponding \'-genome\' and \'-gff\'.')
