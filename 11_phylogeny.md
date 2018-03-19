# Phylogeny

First I am making a new vcf file that has only variable sites with no missing data:

```
zcat ../F_and_M/FandM_chrX_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz | /mnt/expressions/ben_evans/bin/vcftools/bin/vcftools --vcf - --max-missing-count 0 --out ../F_and_M/FandM_chrX_BSQR_jointgeno_allsites_filtered_SNPsonly_nomissing --recode 
```

Then I can convert this to a tab file and use my script to make this into a phylip file.
