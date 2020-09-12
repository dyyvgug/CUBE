# CAFE：Codon Adaptation Facile Estimation
### CAFE can calculate the mCAI(modified Codon Adaptation Index) value, and optimize gene sequences to increase expression.
##### Created By: Yingying Dong
##### Email: dyyvgug@163.com

&#8195;&#8195;Before using, please make sure that python 3.X has been installed on your computer. When using, download the repository to the local.

&#8195;&#8195;CAFE is mainly divided into two parts.The first part is used to calculate the mCAI value, and the second part is used to optimize the gene sequence to increase the amount of gene expression.

**PART1: calculate mCAI value**

&#8195;&#8195;If the species name you plan to calculate is in the **supported_species.txt**, you can use the mCAI.py script to calculate.If the species is not in supported_species.txt, you can use the mCAI_comp.py script to calculate, but the **rpy2** dependent library is difficult to install on Windows.

mCAI.py optional arguments:
```
  -h, --help        show this help message and exit
  -spe [SPE]        The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans
  -inp [INP]        The FASTA file of the gene sequences that wants to calculate the mCAI value
```
mCAI_comp.py optional arguments:
```
  -h, --help        show this help message and exit
  -spe [SPE]        The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans
  -inp [INP]        The FASTA file of gene sequences that you want to calculate the mCAI value
  -genome [GENOME]  The FASTA file of the species genome
  -gff [GFF]        The annotation file GFF3 format of the species
```
&#8195;&#8195;Among them, the ```-spe``` parameter followed by the species that support calculation is in **supported_species.txt**. The ```-inp``` parameter is followed by the gene sequences for which the mCAI value needs to be calculated. The example sequences is in the **example_files folder**.```-genome``` and ```-gff``` parameters are followed by genome sequence and GFF3 annotation file respectively.The file format can refer to the example in the example_files folder.Before using the GFF3 file, it is recommended to use the sed_gff.sh script format.The specific command is ```$ bash sed_gff.sh```.

Usage example:

(1). If the species name I plan to calculate is in the **supported_species.txt**:
On Linux bash or Windows cmd： ```python mCAI.py -spe Caenorhabditis_elegans -inp input_Ce.fa```

If the input sequence is not in the current script folder, remember to add the path.E.g:
On Linux bash:```dyy@Workstation:~$ python mCAI.py -spe Zoogloea_oleivorans -inp /home/disk1/input_Zo.fa```
On Windows cmd:``` C:\Users\dyy> python mCAI.py -spe Zoogloea_oleivorans -inp G:\github\CAFE\example_files\PART1\input_Zo.fa```

(2).If the species is not in supported_species.txt:
On Linux bash or Windows cmd： ```python mCAI_comp.py -inp input_Zo.fa -gff Zoogloea_oleivorans.gff -genome Zoogloea_oleivorans.fna```














