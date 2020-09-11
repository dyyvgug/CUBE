# CAFE：Codon Adaptation Facile Estimation
### CAFE can calculate the mCAI(modified Codon Adaptation Index) value, and optimize gene sequences to increase expression.
##### Created By: Yingying Dong
##### Email: dyyvgug@163.com

&#8195;&#8195;Before using, please make sure that python 3.X has been installed on your computer. When using, download the repository to the local.

&#8195;&#8195;CAFE is mainly divided into two parts.The first part is used to calculate the mCAI value, and the second part is used to optimize the gene sequence to increase the amount of gene expression.

**PART1: calculate mCAI value**

&#8195;&#8195;If the species name you plan to calculate is in the **supported_species.txt**, you can use the mCAI.py script to calculate.If the species is not in supported_species.txt, you can use the mCAI_comp.py script to calculate, but the **rpy2** dependent library is difficult to install on Windows.

mCAI.py optional arguments:

  -h, --help        show this help message and exit
  -spe [SPE]        The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans
  -inp [INP]        The FASTA file of the gene sequence that wants to calculate the mCAI value
  
mCAI_comp.py optional arguments:

  -h, --help        show this help message and exit
  -spe [SPE]        The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans
  -inp [INP]        The FASTA file of genes sequence that you want to calculate the mCAI value
  -genome [GENOME]  The FASTA file of the species genome
  -gff [GFF]        The annotation file GFF3 format of the species













