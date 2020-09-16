#!/usr/bin/python
# -*- coding:utf-8 -*-
# Author: Yingying Dong.

import os
import sys
import platform
import argparse
from scipy import stats


parser = argparse.ArgumentParser(description='Calculate mCAI.', prog='mCAI', usage='%(prog)s [options]')
parser.add_argument('-spe', nargs='?', required=True, type=str, help='The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans')
parser.add_argument('-inp', nargs='?', required=True, type=str, help='The FASTA file of the gene sequences that wants to calculate the mCAI value')
parser.add_argument('-o', nargs='?', type=str, default='mCAI.txt',
                    help='The file name of output mCAI value.The default file name is \'mCAI.txt\'')
args = parser.parse_args()


def cal_mcai(file, species, out):
    syst = platform.system()
    if syst == "Windows":
        os.chdir('.\\')
        we_path = '.\\resource\\weight\\'
    elif syst == "Linux":
        os.chdir('./')
        we_path = './resource/weight/'

    if os.path.exists('{}{}'.format(we_path, species)):
        weight_file = open('{}{}'.format(we_path, species), 'r')
        CAI_file = open(out, 'w')

        weight_table = []
        for line in weight_file:
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
        weight_file.close()
        CAI_file.close()
    else:
        print('\tThe calculation of this species is not supported, and the species that supports calculation are mentioned in the \'supported_species.txt\'.\n\t If you have the genome and GFF annotation files of the species, you can generate weight from the cal_RSCU.py and cal_weight.R file, and then use the script to calculate the mCAI value')


if __name__ == '__main__':
    cal_mcai(args.inp, args.spe, args.o)
