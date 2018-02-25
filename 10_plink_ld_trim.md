# Need to do LD trimming on SNP file

SNP files are here:
```
/mnt/scratch/ben_evans/F_and_M/FandM_chr*_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz*
```

May need to convert to plink format like this:
```
vcftools --gzvcf /mnt/scratch/ben_evans/F_and_M/FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --plink --out /mnt/scratch/ben_evans/F_and_M/FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz_in_plink
```

Use plink like this:
```
plink --file data --indep 50 5 2 --out prefix --double-id
```

for a zipped vcf file, use this:
```
plink --vcf /mnt/scratch/ben_evans/F_and_M/FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --indep 50 5 2 --out /mnt/scratch/ben_evans/F_and_M/FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly_LDprune --double-id
```

The `--double-id` is needed to make plink ignore underscores in the sample IDs

according to : http://zzz.bwh.harvard.edu/plink/summary.shtml#prune

"The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously. That is, this considers the correlations between SNPs but also between linear combinations of SNPs. A VIF of 10 is often taken to represent near collinearity problems in standard multiple regression analyses (i.e. implies R^2 of 0.9). A VIF of 1 would imply that the SNP is completely independent of all other SNPs. Practically, values between 1.5 and 2 should probably be used; particularly in small samples, if this threshold is too low and/or the window size is too large, too many SNPs may be removed."

this will generate two files: 
      plink.prune.in
      plink.prune.out
      
Then save the plink.prune.in SNPs:      
```
plink --file data --extract plink.prune.in --make-bed --out pruneddata
```
