# CUBE：Codon Usage Bias Ensemble
### CUBE local toolkit can calculate some popular CUB indices such as the modified codon adaptation index (mCAI) value, codon bias index (CBI), effective number of codons (ENC), frequency of optimal codons (Fop), GC content at the third position of synonymous codon (GC3s).
### Besides, CUBE toolkit also can optimize gene sequences to increase expression.

##### Created By: Yingying Dong
##### Email: dyyvgug@hotmail.com

&#8195;&#8195;Before using, please make sure that python 3.X has been installed on your computer. When using, download the repository to the local.

&#8195;&#8195;CUBE toolkit is mainly divided into two parts.The first part is used to calculate the CUB indices, and the second part is used to optimize the gene sequence to increase the amount of gene expression.

-----
**PART1: calculate CUB indices**

&#8195;&#8195;If the species name you plan to calculate is in the **supported_species.txt**, you can use the cube.py script to calculate.If the species is not in supported_species.txt, you can use the cube_comp.py script to calculate, but the **rpy2** dependent library is difficult to install on Windows.

cube.py optional arguments:
```
  -h, --help        show this help message and exit
  -cub [CUB]        The list of CUB indices that you want to calculate, for example: [CAI,CBI,gc3s]
  -spe [SPE]        The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans
  -inp [INP]        The FASTA file of the gene sequences that wants to calculate the cube value
  -o [O]            The file name of output CUB value.The default file name is 'cube.txt'
```
cube_comp.py optional arguments:
```
  -h, --help        show this help message and exit
  -spe [SPE]        The Latin name of the species, separated by an underscore, for example: Caenorhabditis_elegans
  -inp [INP]        The FASTA file of gene sequences that you want to calculate the cube value
  -genome [GENOME]  The FASTA file of the species genome
  -gff [GFF]        The annotation file GFF3 format of the species
  -cub [CUB]        The list of CUB indices that you want to calculate, for example: [CAI,CBI,gc3s]
  -o [O]            The file name of output CUB value.The default file name is 'cube.txt'
```
&#8195;&#8195;Among them, the ```-spe``` parameter followed by the species that support calculation is in **supported_species.txt**. The ```-inp``` parameter is followed by the gene sequences for which the CUB indices need to be calculated. The example sequences is in the **example_files folder**.```-genome``` and ```-gff``` parameters are followed by genome sequence and GFF3 annotation file respectively.The file format can refer to the example in the example_files folder.Before using the GFF3 file, it is recommended to use the sed_gff.sh script format.The specific command is ```$ bash sed_gff.sh```.

##############

Usage example:

##############

(1). If the species name I plan to calculate is in the **supported_species.txt**:

On Linux bash or Windows cmd： ```python cube.py -spe Caenorhabditis_elegans -inp input_Ce.fa```

&#8195;&#8195;If the input sequence is not in the current script folder, remember to add the path.E.g:

&#8195;&#8195;On Linux bash:

&#8195;&#8195;```dyy@Workstation:~$ python cube.py -spe Zoogloea_oleivorans -inp /home/disk1/input_Zo.fa```

&#8195;&#8195;On Windows cmd:

&#8195;&#8195;``` C:\Users\dyy> python cube.py -spe Zoogloea_oleivorans -inp G:\github\CAFE\example_files\PART1\input_Zo.fa```

(2).If the species is **not** in supported_species.txt:

On Linux bash or Windows cmd： ```python cube_comp.py -inp input_Zo.fa -gff Zoogloea_oleivorans.gff -genome Zoogloea_oleivorans.fna```

-----
**PART2: optimize gene sequences**

&#8195;&#8195;This optSeq_comp.py script is used to optimize gene sequences in order to increase the expression of heterologous proteins and endogenous genes.

optSeq_comp.py optional arguments:
```
  -h, --help            show this help message and exit
  -inp [INP]            (Required Parameters) The file name of the original sequence. The
                        sequence default type is DNA sequence, if that is protein sequence,
                        please add '-Pro' parameter
  -DNA                  FASTA file for DNA sequences of genes that wish to increase expression
  -Pro                  FASTA file for protein sequences of genes that wish to increase
                        expression
  -spe [SPE]            Latin name of host species, separated by an underscore, for example:
                        Caenorhabditis_elegans
  -genome [GENOME]      The FASTA file of the species genome
  -gff [GFF]            The annotation file GFF3 format of the species
  -o [O]                The file name of output optimized sequence.The default file name is
                        'optimized_seq.fa'
  -poly                 Need to remove polyN sequence
  -res                  Need to remove restriction enzyme sites
  -res_sites1 [RES_SITES1]
                        Restriction enzyme recognition sequence
  -res_sites2 [RES_SITES2]
                        Restriction enzyme recognition sequence
  -res_sites3 [RES_SITES3]
                        Restriction enzyme recognition sequence
 ```
&#8195;&#8195;Among them,if the ```-DNA``` and ```-Pro``` parameters are not added, the default is ```-DNA```.That is, the input original gene sequences are DNA sequences.Here are examples of usage：

##############

Usage example:

##############

(1).If I want the host species of gene expression in **supported_species.txt**:

```python optSeq_comp.py -spe Komagataella_pastoris -DNA -inp input_original_DNA.fa```

(2).If the input sequence file contains protein sequences：

```python optSeq_comp.py -spe Komagataella_pastoris -Pro -inp input_original_protein.fa```

(3).If you want optimized sequences to avoid the polyN signal (polyA is known to cause premature termination of transcription and translation):

```python optSeq_comp.py -spe Komagataella_pastoris -DNA -inp input_original_DNA.fa -poly```

(4).If you want optimized sequences without corresponding restriction enzyme sites：

```python optSeq_comp.py -spe Komagataella_pastoris -DNA -inp input_original_DNA.fa -res -res_sites1 "GAATTC"```

(5).If the species is **not** in supported_species.txt:

```python optSeq_comp.py -inp input_original_DNA.fa -gff Zoogloea_oleivorans.gff -genome Zoogloea_oleivorans.fna```

------
&#8195;&#8195;Thank you and best wishes！
