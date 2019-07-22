Code_challenge
----
Code was written by R

Please install package "stringr" before running
```
install.packages('stringr')
```
gtf file was downloaded from ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz

When you use it in terminal of linux/mac, the first two arguments will be VCF file name and gtf/gff file name.
e.g.
```
Rscript code_challenge_templus.R /Users/mg/Downloads/Challenge_data.vcf /Users/mg/Downloads/Homo_sapiens.GRCh37.87.chr.gtf 
```
Input: VCF file and gtf/gff file

Output: annotated VCF file with "TYPE","DP/AO/r_AO" and "info_ExAC"

"TYPE" is type of variants

"DP/AO/r_AO" is deep of reads, variants reads and percentage of variants reads V.S. reference reads 

"info_ExAC" is allele frequence, Gene ID and Existing variation ID from ExAC
