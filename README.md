# cleanSeq
* cleanSeq is a bioinformatic tool for contamination detection, contamination removal, mutation detection, and mutation verification.

# prework
## step1: install dependence tool
```shell
BWA (v0.7.17-r1188)
Samtools (v1.10)
GATK (v4.1.8.1)
BLAST+ (v1.12.0+)
Clustal Omega (v1.2.4)
FLASH (v1.2.11)
trimmomatic (v0.38)
```

## step2: bulid ntdatabase
```shell
download blast database: https://ftp.ncbi.nih.gov/blast/db/FASTA/
bulid blastdb: "makeblastdb -dbtype nucl -in ntdatabase -input_type fasta -out nt.blastdb -parse_seqids"

note: raw nt database contain 'N' characters, please remove it before using 'makeblastdb' command bulid blast database
```
## step3: bulid python environment (if you choose run cleanSeq.py)
```shell
conda create -n cleanSeq python==3.8
conda activate cleanSeq
conda install -c conda-forge biopython (v:1.79)
conda install -c anaconda pandas (v:1.1.3)
conda install -c anaconda reportlab (v:3.5.51)
conda install -c conda-forge bs4 (v:4.10.0)
conda install -c anaconda requests (v:2.24.0)
conda install -c anaconda openpyxl (v:3.0.5)
conda install -c conda-forge matplotlib (v:3.2.2)
conda install -c anaconda xlrd (v:1.79)
```

# Get cleanSeq
## download binary 
* This binary is only for Linux systems: 
  * download: https://pan.baidu.com/s/1hn277ozXEBKcvTlhGBBDbw  
  * extract code: cxgf

```shell
command: chmod a+x ./cleanSeq
         ./cleanSeq reference.fa raw1.fq.gz raw2.fq.gz ntPath
```

## or download cleanSeq.py (pre-bulid python environment)
* only for Linux systems:
    * download: https://github.com/bingxiao-wcy/cleanSeq/blob/main/cleanSeq.py

```shell
command: cleanSeq command: python cleanSeq.py reference.fa raw1.fq.gz raw2.fq.gz ntPath
```

# test cleanSeq
## example
* get TestData 
  * download: https://pan.baidu.com/s/1NNxzAXQZNM5ZEP_NomNw4w 
  * extract code: 8vsk
* run command: ./CleanSeq Escherichia_coli.fasta EcoliPsu.read1.fastq EcoliPsu.read2.fastq nt.blastdb -identity 90
* outptut result: https://github.com/bingxiao-wcy/cleanSeq/blob/main/output.rar



