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
../bin/vcftools/bin/vcftools --gzvcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --keep tonk_individuals.txt --chr chr12 --recode --out ../SEAsian_macaques_bam/females_and_males/Fa
```

The `--recode ` option is needed to make it output the file. The `--out` option specifies a directory and the output file is written to this directory with the prefix of the input file followed by 'recode.vcf'

* then thin the vcf file to include only one variable position within a 15 bp window. I modified a script that Laurie wrote ('thinVCF.pl')  take gz files as input. But I don't need this because vcftools writes a vcf file (which I should later delete).

```
./thinVCF.pl ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz.recode.vcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz_thinned.vcf
```
The intermediate file should be deleted.
```
rm -f ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz.recode.vcf 
```

The thinned file should be compressed:

```
../bin/htslib-1.6/bin/bgzip ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz_thinned.vcf
```
No indexing needed, but if it was I would do it like this:

```
../bin/htslib-1.6/bin/tabix -p vcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz_thinned.vcf.gz
```

* Now make the LDHat input files

```
../bin/vcftools/bin/vcftools --gzvcf ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz_thinned.vcf.gz --chr chr12 --ldhat-geno --out ../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz_thinned
```
(no longer need --keep tonk_individuals.txt  because these were selected earlier)

This should generate two output files with suffixes  `.ldhat.sites` and `.ldhat.locs`. The `--keep` command refers to a file with a list of the individuals to include but this is not needed because we already selected them in an earlier step. The `--chr` option is required, I think, even though we are only inputting one chr anyhow. The `--out` option provides a prefix for the output files.
  

  
  
# Lookup files

LDhat needs lookupfiles that specify the number of sequences in the sample, theta per site, and a 'grid size' which is the maximum value of 4Ner and has a suggested value of 100. In the manual it says "It is worth noting that minor changes in the value of theta per site do not seem to have a large influence on the estimated recombination rate".  So I plan to use lkgen to make these files using the existing file (lk_n50_t0.001) for theta per site equal to 0.001 (which is pretty close to the estimate for M. tonkeana of 0.002 and probably simialr to the esimates for the other species.

This file is here:
```
/work/ben/2017_SEAsian_macaques/bin/LDhat/LDhat_lookup/lk_n50_t0.001
```

```
/work/ben/2017_SEAsian_macaques/bin/LDhat/lkgen -lk /work/ben/2017_SEAsian_macaques/bin/LDhat/LDhat_lookup/lk_n50_t0.001 -nseq 12
```
This generates a lookup file `new_lk.txt` that should be renamed to `/work/ben/2017_SEAsian_macaques/bin/LDhat/LDhat_lookup/lk_n12_t0.001`

Now I have the input files:
```
../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz_thinned.ldhat.locs
../SEAsian_macaques_bam/females_and_males/FandM_chr12_BSQR_jointgeno_allsites_filtered_SNPsonly_tonk.vcf.gz_thinned.ldhat.sites
/work/ben/2017_SEAsian_macaques/bin/LDhat/LDhat_lookup/lk_n12_t0.001
```
So I can run LDhat interval
```
/work/ben/2017_SEAsian_macaques/bin/LDhat/interval -seq <file_name> -loc <file_name> -lk <file_name>
```

