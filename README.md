# Get cleanSeq
## dependence tool
```shell
BWA v0.7.17-r1188
Samtools v1.10
GATK v4.1.8.1
BLAST+ v1.12.0+
Clustal Omega v1.2.4
FLASH v1.2.11
trimmomatic v0.38
```
## Execute from source code
```shell
# bulid python environment.
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
* Command to test
```shell
python cleanSeq.py reference.fa raw1.fq.gz raw2.fq.gz ntPath
```

## or download binary 
This binary is only for Linux systems: 
link: https://pan.baidu.com/s/1hn277ozXEBKcvTlhGBBDbw  extract code: cxgf
```shell
chmod a+x ./cleanSeq
```

* Command to test
```shell
./cleanSeq reference.fa raw1.fq.gz raw2.fq.gz ntPath
```

# Take a quick glance
* Sample PDF report: https://github.com/bingxiao-wcy/cleanSeq/blob/main/Report.pdf
* Dataset for testing: http:



