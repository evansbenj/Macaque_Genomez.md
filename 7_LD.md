# LDHat

I would like to calculate recombination rates for each species.  LDHat can do this using unphased poolymorphism data.  

# Generating input files

Vcftools can generate input files for LDHat.  I am going to try this out using this file:
```
../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz
```
```
vcftools --vcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --out --ldhat-geno 
```

This should generate two output files with suffixes  ".ldhat.sites" and ".ldhat.locs",
