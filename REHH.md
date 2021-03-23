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
Here are the first SNPs in each gene
```
chr03	84594254	84613076	NDUFB5 84594255
chr03	104846075	104847927	NDUFAF3  104846083
chr03	148084830	148114936	ACAD9  148085019
chr03	156309713	156315708	NDUFB4 156309746
chr03	157411293	157436295	TIMMDC1 157411323
```
Some (all?) of these are not variable in each species and the first one isn't even variable in the taxa at all (it just differs from the ref).  But this should be ok.

Now I am going to edit a version of the phased VCF file to add names to these SNPs (I'll call them each by the gene acronym).

Then I compressed it:
```
bgzip -c chr03_phased_named_SNPs.vcf > chr03_phased_named_SNPs.vcf.gz
```


