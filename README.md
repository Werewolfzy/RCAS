# RCAS (Recombinant Analysis)
`RCAS` is a command line tool used for the fine-mapped analyses. 

A large number of SNPS were reduced to obtain new recombination loci based on the standard of allele frequency difference greater than 0.3 (default) between offspring and parents.

## Getting Started
0.Git clone RCAS

In order to download `RCAS`, you should clone this repository via the commands:
```
git clone https://github.com/Werewolfzy/RCAS.git
cd RCAS
```

1.Prepare the environment：

This script depends on the Python3 software and the Plink software and requires the following dependency packages:
```
pip install -r requirements.txt
```

2.Prepare the genotype file, the phenotype file and the significant snp file.

The genotype file is BED format of PLink software. The phenotype file and the significant snp file refer to the test_phenotype.txt and the test_significant.txt. In the phenotype file, 1 and 2 are treated as two different progenitor populations, and 3 is its progenitor population. By comparing 1 and 2, the genotypes with significant loci in the 3 populations are extracted.

3.The command line
```
python chongzu.py -g [genotype_file] -p [phenotype_file] -s [significant_snp] -f [freq_cutoff  (The default is 0.3)] -o [out_file  (The default is result.txt)]
```









