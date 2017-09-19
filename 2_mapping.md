# Toy data for testing

```
zcat PF505_all_R2scythe_and_trimm_paired.cor.fastq.gz | head -n 1000 > toyR2.fastq 
```
problems with misnamed reads fixed with bbmap like this (from within each directory):
```
bbmap/bbmap/repair.sh in1=toyR1.fastq in2=toyR2.fastq out1=toy1fixed1.fq out2=toy2fixed2.fq outsingle=toysingle.fq
```
