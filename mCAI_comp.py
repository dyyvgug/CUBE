#!/usr/bin/python
# -*- coding:utf-8 -*-

# =================================================================================
# Author: Yingying Dong. Email: dyyvgug@163.com .This script is used to calculate
#  modified CAI(mCAI) value.
# =================================================================================

import os
import sys
import platform
import argparse
try:
    from scipy import stats
except:
    os.system('pip install scipy')
    from scipy import stats
import get_weight as rs


parser = argparse.ArgumentParser(description='Calculate mCAI.', prog='mCAI', usage='%(prog)s [options]')
parser.add_argument('-spe', nargs='?', type=str, help='The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans')
parser.add_argument('-inp', nargs='?', required=True, type=str, help='The FASTA file of genes sequence that you want to calculate the mCAI value')
parser.add_argument('-genome', nargs='?', type=str, help='The FASTA file of the species genome')
parser.add_argument('-gff', nargs='?', type=str, help='The annotation file GFF3 format of the species')
args = parser.parse_args()


def cal_mcai(file, species):
    CAI_file = open('mCAI.txt', 'w')
    weight_table = []
    for line in species:
        weight_table.append(line.strip().split('\t'))
    codon_weight = {}
    for i in weight_table:
        codon_weight[i[0]] = float(i[1])

    dna = ''
    weight_list = []
    CAI_file.write('gene_id\tmCAI_value\n')

    f = open(file, 'r')
    for line in f:
        if line.startswith('>') and dna == '':
            header = line.strip().replace('>', '')
        elif not line.startswith('>'):
            dna = str.upper(dna) + line.strip()
        elif line.startswith('>') and dna != '':
            for j in range(0, len(dna), 3):
                codon = dna[j:j + 3]
                if codon in codon_weight:
                    weight_list.append(codon_weight[codon])
            CAI = stats.gmean(weight_list)
            CAI_file.write('{}\t{}\n'.format(header, CAI))
            header = line.strip().replace('>', '')
            dna = ''
            weight_list = []

    f.close()
    CAI_file.close()


if __name__ == '__main__':
    #print("Current path is %s" % (os.path.abspath(sys.argv[0])))
    if args.spe is None and args.genome is not None and args.gff is not None:
        rrs = rs.read_file(args.genome, args.gff)
        cal_mcai(args.inp, rrs)
    elif args.spe is None and args.genome is None and args.gff is None:
        print('\tThere are no parameters for \'-spe\',\'-genome\', and \'-gff\'.\n\tIf the species you need to calculate is in \'sup_spe.txt\', please follow the name after the \'-spe\' parameter.\n\tIf the calculation of the species is not supported, use \'-genome\' and \'-gff\' followed by genome and annotation files')
    elif args.spe is not None and args.genome is None and args.gff is None:
        syst = platform.system()
        if syst == "Windows":
            os.chdir('.\\')
            we_path = '.\\resource\\weight\\'
        elif syst == "Linux":
            os.chdir('./')
            we_path = './resource/weight/'

        if os.path.exists('{}{}_weight.txt'.format(we_path, args.spe)):
            wei_file = open('{}{}_weight.txt'.format(we_path, args.spe), 'r')
            cal_mcai(args.inp, wei_file)
        else:
            print('\tThe calculation of this species is not supported, and the species that supports calculation are mentioned in the \'sup_spe.txt\'.\n\t If you have the genome and GFF annotation files of the species, you can use \'-gff\' and \'-genome\'.\n\t\'-gff\' and \'-genome\' are followed by the annotation file GFF format of the species and the fasta format file of the genome respectively.')
    else:
        print('\tMissing parameters.\nIf you do not have the species name of the \'-spe\' parameter, you must provide the corresponding \'-genome\' and \'-gff\'.')





