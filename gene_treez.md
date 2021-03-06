# Gene treez
directory with full genotype files (not just SNPS):
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio
```
Get the coordinates of some interesting genes:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'MRPL1'
```
Here they are:
```
NDUFAF3  chr03	104846075	104847927
NDUFA2 chr05	138178553	138180781
TACO1 chr17	56920802	56928897
HARS2 chr05	138224999	138233197
COX16 chr14	132330617	132354632
C11orf83   chr11	11484179	11486231
MRPS21  chr01	124258900	124272013
UQCRC1   chr03	104359835	104371673
ATP5B    chr12	54910245	54920096
MRPL48  chr11	65083808	65165565
NDUFB11   chrX	47194137	47196604
NDUFB5  chr03	84594254	84613076
MRPL1  chr04	56478980	56572561
MRPL47 chr03	84577965	84594079
```
Extract section of gene with vcftools:
```
module load StdEnv/2020 vcftools/0.1.16
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered1.vcf.gz --chr chr01 --from-bp 124258900 --to-bp 124272013 --out MRPS21.vcf.gz --recode
```






vcftools --gzvcf FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered1.vcf.gz --chr chr05 --from-bp 138224999 --to-bp 138233197 --out HARS2.vcf.gz --recode

