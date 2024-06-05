#!/usr/bin/python
# -*- coding:utf-8 -*-
# Author: Yingying Dong.
# =================================================================================
# Author: Yingying Dong. Email: dyyvgug@163.com .This script is used to calculate
#  CUB indices.
# =================================================================================

import os
import sys
import platform
import argparse
from scipy import stats
import codonw


parser = argparse.ArgumentParser(description='Calculate CUB indices.', prog='CUB', usage='%(prog)s [options]')
parser.add_argument('-spe', nargs='?', required=True, type=str, default='Caenorhabditis_elegans', help='The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans')
parser.add_argument('-i', nargs='?', required=True, type=argparse.FileType('r'), help='The FASTA file of the gene sequences that you want to calculate CUB values')
parser.add_argument('-o', nargs='?', type=str, default='cub.txt',
                    help='The file name of output CUB value.The default file name is \'cub.txt\'')
parser.add_argument('-cub',nargs='?', type=list or str, default=["CAI","ENC"],
                    help='The CUB indices you want to calculate, you can input one or more indices, such as ["CAI",""]')
args = parser.parse_args()


def cal_cub(dataSource, species, output,indices):
    syst = platform.system()
    if syst == "Windows":
        os.chdir('.\\')
        we_path = '.\\resource\\weight\\'
    elif syst == "Linux":
        os.chdir('./')
        we_path = './resource/weight/'

    if os.path.exists('{}{}'.format(we_path, species)):
        weight_file = open('{}{}'.format(we_path, species), 'r')
        result = open(output, 'w+')

        weight_table = []
        for line in weight_file:
            weight_table.append(line.strip().split('\t'))
        codon_weight = {}
        for i in weight_table:
            codon_weight[i[0]] = float(i[1])

        dna = ''
        header = ''
        weight_list = []
        result.write("gene_id\t")
        result.write("\t".join(indices))
        result.write("\n")
        indices = [i.lower() for i in indices]

        dataSource += '\n>'
        f = dataSource.split('\n')

        for line in f:
            if line.startswith(">") and dna == "":
                header = line.strip().replace(">", "")
            elif not line.startswith(">"):
                dna = str.upper(dna) + line.strip()
            elif line.startswith(">") and dna != "":
                for j in range(0, len(dna), 3):
                    codon = dna[j: j + 3]
                    if codon in codon_weight:
                        weight_list.append(codon_weight[codon])
                # print(type(dna))
                CAI = stats.gmean(weight_list)
                index_list = []
                cseq = codonw.CodonSeq(dna)
                for i in indices:
                    if i == "cai":
                        index_list.append(CAI)
                    elif i == "gc3s":
                        index_list.append(cseq.bases2()['GC3s'])
                    elif i == "gc":
                        index_list.append(cseq.bases2()['GC'])
                    elif i == "cbi":
                        if species == "Escherichia_coli":
                            index_list.append(cseq.cbi())
                        elif species == "Bacillus_subtilis":
                            index_list.append(cseq.cbi(1))
                        elif species == "Saccharomyces_cerevisiae_S288C":
                            index_list.append(cseq.cbi(2))
                        else:
                            index_list.append("NA")
                    elif i == "fop":
                        if species == "Escherichia_coli":
                            index_list.append(cseq.fop())
                        elif species == "Bacillus_subtilis":
                            index_list.append(cseq.fop(1))
                        elif species == "Dictyostelium_discoideum":
                            index_list.append(cseq.fop(2))
                        elif species == "Aspergillus_nidulans":
                            index_list.append(cseq.fop(3))
                        elif species == "Saccharomyces_cerevisiae_S288C":
                            index_list.append(cseq.fop(4))
                        elif species == "Drosophila_melanogaster":
                            index_list.append(cseq.fop(5))
                        elif species == "Caenorhabditis_elegans":
                            index_list.append(cseq.fop(6))
                        elif species == "Neurospora_crassa":
                            index_list.append(cseq.fop(7))
                        else:
                            index_list.append("NA")
                    else:
                        index_list.append(getattr(cseq, i)())
                index_list = [str(num) for num in index_list]

                result.write('{}\t'.format(header.upper()))
                result.write("\t".join([index.upper() for index in index_list]))
                result.write('\n')
                header = line.strip().replace('>', '')
                dna = ""
                weight_list = []

        weight_file.close()

    else:
        print('\tThe calculation of this species is not supported, and the species that supports calculation are mentioned in the \'supported_species.txt\'.\n\t If you have the genome and GFF annotation files of the species, you can generate weight from the cal_RSCU.py and cal_weight.R file, and then use the script to calculate the mCAI value')


if __name__ == '__main__':
    cal_cub(dataSource=args.i.read(), species=args.spe, output=args.o,indices=args.cub)
