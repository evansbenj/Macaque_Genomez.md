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

There are two "Performs' scripts for chrX.  The analysis is best done on the haploid calls using "Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl" because this one checks for the length of the call in a haploid.  For "Performs_ABBA_BABA_on_populations_onlychrX.pl" it is possible for a haploid call to be a microsat and still be counted. This detail is probably irrelevant but perhaps best to fix or not use diploid calls at all.

```
sqsub -r 2d --mpp 16G -o abba_H3_hecki_H1maura_H2tonk_allhaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl  /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_allhaploid.vcf.gz.tab 00000000000000000000000 3_4_1-2-3-4-5_6-7-8-9-10_18-19-20-21-22-23 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_allhaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_allhaploid_chrX.stats"
```
Or to include diploid genotypes from females (see note about about microsats:

```
sqsub -r 2d --mpp 16G -o abba_H3_hecki_H1maura_H2tonk_females_diploid_maleshaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX.pl /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_females_and_males_haploid.vcf.gz.tab 11111110001100001111110 3_4_1-2-3-4-5_6-7-8-9-10_18-19-20-21-22-23 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_females_diploid_maleshaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_hecki_H1_maura_H2_tonk_females_diploid_maleshaploid_chrX.stats"
```
# Other Preliminary chrX analyses:
H3tonk_H1nigra_H2hecki
```
sqsub -r 2d --mpp 16G -o abba_H3_tonk_H1_nigra_H2_hecki_allhaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl  /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_allhaploid.vcf.gz.tab 00000000000000000000000 3_4_18-19-20-21-22-23_17_1-2-3-4-5 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_tonk_H1_nigra_H2_hecki_allhaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_tonk_H1_nigra_H2_hecki_allhaploid_chrX.stats"
```
H3_borneo_no665_no1206_H1_nigra_H2_heck
```
sqsub -r 2d --mpp 16G -o abba_H3_borneo_no665_no1206_H1_nigra_H2_heck_allhaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl  /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_allhaploid.vcf.gz.tab 00000000000000000000000 3_4_11-14-16_17_1-2-3-4-5 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_heck_allhaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_heck_allhaploid_chrX.stats"
```

H3_borneo_no665_no1206_H1_nigra_H2_heck
```
sqsub -r 2d --mpp 16G -o abba_H3_borneo_no665_no1206_H1_nigra_H2_heck_allhaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl  /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_allhaploid.vcf.gz.tab 00000000000000000000000 3_4_11-14-16_17_1-2-3-4-5 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_heck_allhaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_heck_allhaploid_chrX.stats"
```
H3_borneo_no665_no1206_H1_nigra_H2_tonk
```
sqsub -r 2d --mpp 16G -o abba_H3_borneo_no665_no1206_H1_nigra_H2_tonk_allhaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl  /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_allhaploid.vcf.gz.tab 00000000000000000000000 3_4_11-14-16_17_18-19-20-21-22-23 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_tonk_allhaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_tonk_allhaploid_chrX.stats"
```
H3_borneo_no665_no1206_H1_nigra_H2_maura
```
sqsub -r 2d --mpp 16G -o abba_H3_borneo_no665_no1206_H1_nigra_H2_maura_allhaploid_chrX.log bash -c "perl Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup.pl  /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_allhaploid.vcf.gz.tab 00000000000000000000000 3_4_11-14-16_17_6-7-8-9-10 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_maura_allhaploid_chrX.abbababa /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/H3_borneo_no665_no1206_H1_nigra_H2_maura_allhaploid_chrX.stats"
```

