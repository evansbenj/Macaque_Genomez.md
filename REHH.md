# Rehh  

I am going to try to use the R package rehh to detect long runs of homozygosity.

First step is to find some SNPs in each N_interact gene and name them in a vcf file.

I am working in this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/ROH
```
First I made a bed file (`Ninteract_chr03.bed`)that has the ranges of all N_interact genes on chr03. I got this information from the gff file (more conveniently summarized in my table). Paste it in this file:
```
chr03	84594254	84613076	NDUFB5
chr03	104846075	104847927	NDUFAF3
chr03	148084830	148114936	ACAD9
chr03	156309713	156315708	NDUFB4
chr03	157411293	157436295	TIMMDC1
```
Then extract the SNPs like this:
```
module load tabix
tabix FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz -B Ninteract_chr03.bed > Ninteract_chr03_SNPs.vcf
```
