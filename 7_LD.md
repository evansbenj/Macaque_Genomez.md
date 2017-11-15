# LDHat

I would like to calculate recombination rates for each species.  LDHat can do this using unphased poolymorphism data.  

# Generating input files

Vcftools can generate input files for LDHat.  I am going to try this out using this file:
```
../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz
```
on iqaluk

* first filter out the species that will be analyzed (will have to do this for each chr
```
../bin/vcftools/bin/vcftools --gzvcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --keep tonk_individuals.txt --chr chr12 --out ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz
```

Need to check if this is really zipped.

* then thin the vcf file to include only one variable position within a 15 bp window. I modified a script that Laurie wrote ('thinVCF.pl')  take gz files as input.

```
./thinVCF.pl ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.gz ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk_thinned.vcf
```

Then this needs to be compressed and indexed.

```
../bin/htslib-1.6/bin/bgzip ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk_thinned.vcf
```
```
../bin/htslib-1.6/bin/tabix -p vcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk_thinned.vcf.gz
```

* Now make the LDHat input files

```
../bin/vcftools/bin/vcftools --gzvcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk_thinned.vcf.gz --chr chr12 --ldhat-geno --out ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk_thinned 
```
(no longer need --keep tonk_individuals.txt  because these were selected earlier)

This should generate two output files with suffixes  `.ldhat.sites` and `.ldhat.locs`. The `--keep` command refers to a file with a list of the individuals to include. The `--chr` option is needed because each chr is in a different linkage group and make sure to use the entire name of the chr or the output file will have no sites. The `--out` option provides a prefix for the output files.

  I'm trying now to pipe the output to bgzip to keep it small prior to inputting into vcftools.  One concern is that I probably should divide up the vcf files by taxon before thinning; that way I don't remove thin based on sites that are divergent between species in the multisample vcf.
