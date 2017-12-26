# ChrX can have females diploid and males haploid by depth or both haploid by depth

* Genotyping by depth
** males only by depth
```
sqsub -r 2d --mpp 16G -o chrX_geno_depth_onlymaleshaploid.log bash -c "perl 10_Genotypes_only_male_chrX_based_on_allelic_depth.pl /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered.vcf.gz 11111110001100001111110 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_females_and_males_haploid.vcf.gz.tab"
```
** both by depth
```
sqsub -r 2d --mpp 16G -o chrX_geno_depth_allhaploid.log bash -c "perl 10_Genotypes_only_male_chrX_based_on_allelic_depth.pl /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered.vcf.gz 00000000000000000000000 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_allhaploid.vcf.gz.tab"
```

# Abbababa for chrX

(I still need to figure out what is the difference between these Performs scripts)

```
sqsub -r 2d --mpp 16G -o abba_H3_hecki_H1maura_H2tonk_allhaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_females_and_males_haploid.vcf.gz.tab 00000000000000000000000 3_4_1-2-3-4-5_6-7-8-9-10_18-19-20-21-22-23 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_allhaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_allhaploid_chrX.stats"
```
```
sqsub -r 2d --mpp 16G -o abba_H3_hecki_H1maura_H2tonk_females_diploid_maleshaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX.pl /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_females_and_males_haploid.vcf.gz.tab 11111110001100001111110 3_4_1-2-3-4-5_6-7-8-9-10_18-19-20-21-22-23 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_females_diploid_maleshaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_females_diploid_maleshaploid_chrX.stats"
```
