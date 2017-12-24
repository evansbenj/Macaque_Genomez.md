# Abbababa tests  

I have some files on goblin here:
```
/work/ben/2017_SEAsian_macaques/
```

* Genotype chrX by depth

from here `/work/ben/2017_SEAsian_macaques/ben_scripts`

```
sqsub -r 2d --mpp 16G -o chrX_geno_dept.log bash -c "perl 10_Genotypes_only_male_chrX_based_on_allelic_depth.pl /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered.vcf.gz 00000000000000000000000 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_females_and_males_haploid.vcf.gz.tab"
```

* Run abbababa on autosomes

```
./Wrapper_for_Performs_ABBA_BABA_populations_H3_hecki_H1_maura_H2_tonk.pl
```
