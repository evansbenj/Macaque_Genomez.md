# LD pruning on vcf SNP-only genotypes

SNP files are here:
```
/mnt/scratch/ben_evans/F_and_M/FandM_chr*_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz*
```

Convert to plink format like this:
```
vcftools --gzvcf /mnt/scratch/ben_evans/F_and_M/FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --plink --out /mnt/scratch/ben_evans/F_and_M/FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz_in_plink
```
This preserves biallelic SNPs only. THen use plink to identify SNPs that are not in LD, like this:
```
plink --file data --indep 50 5 2 --out prefix --double-id
```
for a zipped vcf file, this is the command I used:
```
plink --file FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz_in_plink --indep 50 5 2 --out /mnt/scratch/ben_evans/F_and_M/FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprune --double-id
```

The `--double-id` was needed to make plink ignore underscores in the sample IDs

according to : http://zzz.bwh.harvard.edu/plink/summary.shtml#prune

"The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously. That is, this considers the correlations between SNPs but also between linear combinations of SNPs. A VIF of 10 is often taken to represent near collinearity problems in standard multiple regression analyses (i.e. implies R^2 of 0.9). A VIF of 1 would imply that the SNP is completely independent of all other SNPs. Practically, values between 1.5 and 2 should probably be used; particularly in small samples, if this threshold is too low and/or the window size is too large, too many SNPs may be removed."

this generated two files for each chr: 
      FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprune.prune.in
      FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprune.prune.out
      
Then extract the SNPs listed the FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprune.prune.in file like this:      
```
plink --file data --extract plink.prune.in --make-bed --out pruneddata
```
or for these data like this:
```
plink --file FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz_in_plink --extract FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprune.prune.in --make-bed --out FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprunneddata.out
```
now make a new vcf file that contains only SNPs not in LD:
```
plink --bfile your_ped_map_input --recode vcf
```
```
plink --bfile FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprunneddata.out --recode vcf -out FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprunneddata
```
