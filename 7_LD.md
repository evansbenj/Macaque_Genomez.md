# LDHat

I would like to calculate recombination rates for each species.  LDHat can do this using unphased poolymorphism data.  

# Generating input files

Vcftools can generate input files for LDHat.  I am going to try this out using this file:
```
../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz
```
on iqaluk
```
../bin/vcftools/bin/vcftools --gzvcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --keep tonk_individuals.txt --chr chr12 --ldhat-geno --out ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk 
```

This should generate two output files with suffixes  `.ldhat.sites` and `.ldhat.locs`. The `--keep` command refers to a file with a list of the individuals to include. The `--chr` option is needed because each chr is in a different linkage group and make sure to use the entire name of the chr or the output file will have no sites. The `--out` option provides a prefix for the output files.
