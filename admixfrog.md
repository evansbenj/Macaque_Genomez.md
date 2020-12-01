# Admixfrog

Final Vcf files (before thinning for admixfrog) are here on graham:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data
```
These files should be used for ABBA_BABA, Fst, pi, etc for the genome.

For admixfrog, I need to filter these files so they can be processed in a manageable way.

On graham I installed admixfrog here `/home/ben/.local/bin/admixfrog` as follows:
```
module load scipy-stack/2019b
pip install cython scipy --upgrade --user
pip install git+https://github.com/benjaminpeter/admixfrog --user
```

Now I think it works based on this:
```
module load scipy-stack/2019b
/home/ben/.local/bin/admixfrog --help
```
# Problems
I got errors when I tried to run commands.  So I modified this file:
/home/ben/.local/lib/python3.7/site-packages/admixfrog/interface.py
to have this command `import numpy as np` instead of this one : `import numpy`
but I am still getting another error.  Ugh.

# Filtering

I initially filtered by missingness per species and then thinned based on distance between SNPs as detailed below. Instead, I decided to filter only on the number of missing genotypes being less than or equal to 2 and setting the minimum genotype quality at 30 as follows:
```
vcftools --gzvcf FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz --max-missing-count 2 --minQ 30 --recode --recode-INFO-all --out ./FandM_chr03BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf
```
The files that this generated will be used as input for admix frog; they have names like this: 
```
all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz
```
and they are in this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/maxmissing_2_for_admixfrog
```


This is what I did previously:

on graham, load vcftools:
```
module load nixpkgs/16.09
module load intel/2018.3
module load vcftools/0.1.14
```
Keep variants that have been successfully genotyped in 50% of individuals and a minimum quality score of 30
The --recode flag tells the program to write a new vcf file with the filters, 
The --recode-INFO-all keeps all the INFO flags from the old vcf file in the new one. 
```
vcftools --gzvcf ../FandM_chr01_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --max-missing 0.5 --minQ 30 --recode --recode-INFO-all --out ./FandM_chr01_mm_0.5_minQ_30
```

The next step is to assess individual levels of missing data.
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30_minDP_3.vcf.gz --missing-indv
```

Examine the output file

```
cat out.imiss
```

This shows that, for chr01, papio, nigrescens_PM654, hecki_PF644, hecki_PF647, tonk_PF597 all have a proportion of missing sites respectively of 20.8%, 37.5%, 40.6%, 45.4%, 56.6%. All others have 10.6% (tonk_PF559) or <6.1%.

Now restrict the data to variants called in a high percentage of non-missing genotypes across individuals from every species.
Make a tab delimited file calleed 'species.txt' with the species assignment of each sample.
```
bru_PF707	brunnescens
hecki_PF505	hecki
hecki_PF643	hecki
hecki_PF644	hecki
hecki_PF647	hecki
hecki_PF648	hecki
maura_PF615	maura
maura_PF713	maura
maura_PM613	maura
maura_PM614	maura
maura_PM616	maura
nem_GumGum_female	nemestrina
nem_Ngsang_sumatra_female	nemestrina
nem_PM1206	nemestrina
nem_PM664	nemestrina
nem_PM665	nemestrina
nem_Sukai_male	nemestrina
nigra_PF1001	nigra
nigra_PF660	nigra
nigra_PM1003	nigra
nigrescens_PM1011	nigrescens
nigrescens_PM654	nigrescens
tog_PF549	togeanus
tonk_PF511	tonkeana
tonk_PF559	tonkeana
tonk_PF563	tonkeana
tonk_PF597	tonkeana
tonk_PF626	tonkeana
tonk_PM592	tonkeana
download  papio
```
```
mawk '$2 == "brunnescens"' species.txt > 1.keep && mawk '$2 == "hecki"' species.txt > 2.keep && mawk '$2 == "maura"' species.txt > 3.keep && mawk '$2 == "nemestrina"' species.txt > 4.keep && mawk '$2 == "nigra"' species.txt > 5.keep && mawk '$2 == "nigrescens"' species.txt > 6.keep && mawk '$2 == "togeanus"' species.txt > 7.keep && mawk '$2 == "tonkeana"' species.txt > 8.keep && mawk '$2 == "papio"' species.txt > 9.keep 
```
```
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 1.keep --missing-site --out 1_chr01
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 2.keep --missing-site --out 2_chr01 
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 3.keep --missing-site --out 3_chr01
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 4.keep --missing-site --out 4_chr01 
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 5.keep --missing-site --out 5_chr01
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 6.keep --missing-site --out 6_chr01 
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 7.keep --missing-site --out 7_chr01
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 8.keep --missing-site --out 8_chr01
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --keep 9.keep --missing-site --out 9_chr01
```
Now make a list of the bad loci we want to filter:
```
cat 1_chr01.lmiss 2_chr01.lmiss 3_chr01.lmiss 4_chr01.lmiss 5_chr01.lmiss 6_chr01.lmiss 7_chr01.lmiss 8_chr01.lmiss 9_chr01.lmiss | mawk '!/CHR/' | mawk '$6 > 0.2' | cut -f1,2 >> bad_chr01_loci
```
Now filter these loci:
```
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz --exclude-positions bad_chr01_loci --recode --recode-INFO-all --out FandM_chr01_exclude_missingness
```

After filtering, this reeduced SNPs by about half for chr01:
```
After filtering, kept 3178562 out of a possible 7229006 Sites
```
Now filter these loci (this level reduces the number of SNPs to 10% of the original):
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness.recode.vcf --thin 500 --recode --recode-INFO-all --out FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned
```
now, for chr01: 
```
After filtering, kept 322613 out of a possible 3178562 Sites
```
I now have these in sbatch scrips: `vcftools_1.sh` and `vcftools_2.sh`

# Admixfrog

OK first load baftools and index the filtered file
```
module load bcftools/1.9
bgzip -c FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.recode.vcf > FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
bcftools index FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
tabix -p vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
```


# Admixfrog seems to be working
I was in this directory on graham:
`/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/filtered_for_admixfrog`
But now I am in this directory on graham:
`/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/maxmissing_2_for_admixfrog`

Now create a reference file for each individual.  Before running admixfrog, do this:
```
module load scipy-stack/2019b
```
then for each sample do this:
```
admixfrog-ref [-h] --outfile OUTFILE [--states [STATES [STATES ...]]]
                     [--state-file STATE_FILE] [--cont-id CONT_ID]
                     [--ancestral ANCESTRAL]
                     [--random-read-samples [RANDOM_READ_SAMPLES [RANDOM_READ_SAMPLES ...]]]
                     [--vcf-ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz

```

This works
make the ref file
```
/home/ben/.local/bin/admixfrog-ref --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --out FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --states TON HEC NEM --pop-file pops.yaml 
```
or just do this with `sbatch admixfrog_make_refs.sh all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode_maxmissingcount_2.vcf.recode.vcf.gz` (and modify the vcf filename)

make the target (input) file
```
admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF511sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF511.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF559sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF559.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF563sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF563.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF597sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF597.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF626sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF626.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PM592.in.xz
```
or just do this with `sbatch admixfrog_make_target.sh chrX` (and modify the chr)

run the analysis (for Wallace's Line):
```
sbatch admixfrog_do_analysis_allchrs.sh nem_GumGum_female NEM SUM HEC
sbatch admixfrog_do_analysis_allchrs.sh nem_Ngsang_sumatra_female NEM SUM HEC
sbatch admixfrog_do_analysis_allchrs.sh nem_PM1206 NEM SUM HEC
sbatch admixfrog_do_analysis_allchrs.sh nem_PM664 NEM SUM HEC
sbatch admixfrog_do_analysis_allchrs.sh nem_PM665 NEM SUM HEC
sbatch admixfrog_do_analysis_allchrs.sh nem_Sukai_male NEM SUM HEC


sbatch admixfrog_do_analysis_allchrs.sh nem_GumGum_female NEM SUM TON
sbatch admixfrog_do_analysis_allchrs.sh nem_Ngsang_sumatra_female NEM SUM TON
sbatch admixfrog_do_analysis_allchrs.sh nem_PM1206 NEM SUM TON
sbatch admixfrog_do_analysis_allchrs.sh nem_PM664 NEM SUM TON
sbatch admixfrog_do_analysis_allchrs.sh nem_PM665 NEM SUM TON
sbatch admixfrog_do_analysis_allchrs.sh nem_Sukai_male NEM SUM TON


sbatch admixfrog_do_analysis_allchrs.sh nem_GumGum_female NEM SUM MAU
sbatch admixfrog_do_analysis_allchrs.sh nem_Ngsang_sumatra_female NEM SUM MAU
sbatch admixfrog_do_analysis_allchrs.sh nem_PM1206 NEM SUM MAU
sbatch admixfrog_do_analysis_allchrs.sh nem_PM664 NEM SUM MAU
sbatch admixfrog_do_analysis_allchrs.sh nem_PM665 NEM SUM MAU
sbatch admixfrog_do_analysis_allchrs.sh nem_Sukai_male NEM SUM MAU

```

where the sbatch file `admixfrog_do_analysis_allchrs.sh` is this:
```
#!/bin/sh
#SBATCH --job-name=AF_ref
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=2gb
#SBATCH --output=AF_ref.%J.out
#SBATCH --error=AF_ref.%J.err
#SBATCH --account=def-ben

# BRU HEC MAU NEM SUM NGA NGE TOG TON

bru_PF707   hecki_PF505 hecki_PF643 hecki_PF644 hecki_PF647 hecki_PF648 maura_PF615 maura_PF713 
maura_PM613 maura_PM614 maura_PM616 nem_GumGum   nem_Ngsang   nem_PM1206  nem_PM664   nem_PM665   
nem_Sukai  nigra_PF1001    nigra_PF660 nigra_PM1003    nigrescens_PM1011   nigrescens_PM654    
tog_PF549   tonk_PF511  tonk_PF559  tonk_PF563  tonk_PF597  tonk_PF626  tonk_PM592

# execute like this:
# sbatch admixfrog_do_analysis_allchrs.sh nigra_PF1001 HEC NGA NGE 

module load scipy-stack/2019b
# run the analyses for each chr
admixfrog --infile ${1}.chr01.in.xz --ref FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr01_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr02a.in.xz --ref FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.re
code_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr02a_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont
-est-contamination

admixfrog --infile ${1}.chr02b.in.xz --ref FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.re
code_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr02b_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont
-est-contamination

admixfrog --infile ${1}.chr03.in.xz --ref FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr03_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr04.in.xz --ref FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr04_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr05.in.xz --ref FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr05_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr06.in.xz --ref FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr06_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr07.in.xz --ref FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr07_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr08.in.xz --ref FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr08_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr09.in.xz --ref FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr09_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr10.in.xz --ref FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr10_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr11.in.xz --ref FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr11_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr12.in.xz --ref FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr12_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr13.in.xz --ref FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr13_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr14.in.xz --ref FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr14_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr15.in.xz --ref FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr15_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr16.in.xz --ref FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr16_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr17.in.xz --ref FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr17_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr18.in.xz --ref FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr18_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chr19.in.xz --ref FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.reco
de_maxmissingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chr19_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-es
t-contamination

admixfrog --infile ${1}.chrX.in.xz --ref all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode_max
missingcount_2.vcf.recode.vcf.gz.xz --out ${1}_chrX_${2}_${3}_${4}.out -b 10000 --states ${2} ${3} ${4} --c0 0 --dont-est-conta
mination
```

# Strategy for admixfrog
** within Sulawesi:
I think the best strategy is to focus on adjacent species when possible.  So for nigrescens - nigra and hecki; for hecki - tonk and nigresc; for tonk - hecki, maura; for maura - tonk and tog.

** betweeh Sulawesi and Borneo
I think the best strategy is to compare Borneo to Sumatra plus either hec, tonk, or mau.

# Problems overcome
Something was wrong with nem_PM1206 chr17 but I fixed it by changing the - 10000 flag to -b 20000 as follows:
```
admixfrog --infile nem_PM1206.chr17.in.xz --ref FandM_chr17_mm_0.5_minQ_30_exclude_missingness_thinned.recode.xz --out nem_PM1206_chr17_NEM_SUM_HEC.out -b 20000 --states NEM SUM HEC --c0 0 --dont-est-contamination
```


# Plotting circular plots

```R
## Working directory
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/admixfrog/HEC_NGE_TON")
library(tidyverse)
library(ggplot2)
library(circlize)
library(dplyr)
library(stringr)
# https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
# https://cran.r-project.org/web/packages/circlize/circlize.pdf
# initilize
rm(list=ls()) # removes all variables

# bru_PF707   
# hecki_PF505 hecki_PF643 hecki_PF644 hecki_PF647 hecki_PF648 
# maura_PF615 maura_PF713 maura_PM613 maura_PM614 maura_PM616 
# nem_Ngsang_sumatra_female   
# nem_GumGum   nem_PM1206  nem_PM664   nem_PM665   nem_Sukai   #GumGum is FEMALE; Sukai is a MALE  
# nigra_PF1001    nigra_PF660 nigra_PM1003    
# nigrescens_PM1011   nigrescens_PM654    
# tog_PF549   
# tonk_PF511  tonk_PF559  tonk_PF563  tonk_PF597  tonk_PF626  tonk_PM592

sample_vector <- c("hecki_PF505", "hecki_PF647", "hecki_PF648", "hecki_PF643", "hecki_PF644")
analysis <-"_HEC_NGE_TON"
chrs <- factor(c("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09",
         "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"))

png(paste("hecki",eval(analysis),".png",sep=""),
    width = 500, height = 500, units='mm', res = 300)
# https://jokergoo.github.io/circlize_book/book/introduction.html
# Initialize library




# trying ideogram (this works nicely!)
# https://mran.microsoft.com/snapshot/2014-12-11/web/packages/circlize/vignettes/genomic_plot.pdf
circos.clear()
# start at the top
# insert a space after chrX

circos.par("gap.degree" = c(rep(3, 20),30), 
           "cell.padding" = c(0, 0, 0, 0), 
           "start.degree" = 75) 


# this is a list of chr lengths
chr_length_list = list("1" = 225002135,
                       "2a" = 108967917,
                       "2b" = 131519175,
                       "3" = 198060209,
                       "4" = 190369981,
                       "5" = 179725205,
                       "6" = 171866349,
                       "7" = 185267708,
                       "8" = 144034664,
                       "9" = 111027318,
                       "10" = 129655328,
                       "11" = 127102482,
                       "12" = 133317794,
                       "13" = 95368959,
                       "14" = 169736342,
                       "15" = 92674614,
                       "16" = 74750809,
                       "17" = 76917410,
                       "18" = 70128972,
                       "19" = 53113586,
                       "X" = 148935249
)

# this is a vector of chr lenths
begins <- rep(1,21)
ends = c(225002135,
         108967917,
         131519175,
         198060209,
         190369981,
         179725205,
         171866349,
         185267708,
         144034664,
         111027318,
         129655328,
         127102482,
         133317794,
         95368959,
         169736342,
         92674614,
         74750809,
         76917410,
         70128972,
         53113586,
         148935249)

# this is needed to initalize the graph below
chr_lengths <- cbind(begins,ends)

# chr names
chr_names <- c("chr1","chr2a","chr2b","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
               "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
               "chrX")
chr_names_ordered <- factor(chr_names, ordered = TRUE, 
                            levels = c("chr1","chr2a","chr2b","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                       "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                       "chrX"))

chr_names_ordered_simple <- factor(chr_names, ordered = TRUE, 
                                   levels = c("1","2a","2b","3","4","5","6","7","8","9","10",
                                              "11","12","13","14","15","16","17","18","19",
                                              "X"))


# loop through each sample
for (sample in sample_vector){
  # loop through chrs
  for(i in levels(chrs)){
    print(eval(i))
    a <- read_csv(paste(eval(sample), "_",i,eval(analysis),".out.bin.xz", sep=""))
    assign(i,a)
  }  

  # rename chromosome column
  chr01$chrom <- "chr01"
  chr02a$chrom <- "chr02a"
  chr02b$chrom <- "chr02b"
  chr03$chrom <- "chr03"
  chr04$chrom <- "chr04"
  chr05$chrom <- "chr05"
  chr06$chrom <- "chr06"
  chr07$chrom <- "chr07"
  chr08$chrom <- "chr08"
  chr09$chrom <- "chr09"
  chr10$chrom <- "chr10"
  chr11$chrom <- "chr11"
  chr12$chrom <- "chr12"
  chr13$chrom <- "chr13"
  chr14$chrom <- "chr14"
  chr15$chrom <- "chr15"
  chr16$chrom <- "chr16"
  chr17$chrom <- "chr17"
  chr18$chrom <- "chr18"
  chr19$chrom <- "chr19"
  chrX$chrom <- "chrX"
  
  Allchrs <- rbind(chr01,chr02a,chr02b,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10,
                   chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX)
  
  # make the chr a factor so we can use it for faceting
  Allchrs$chrom <- as.factor(Allchrs$chrom)
  
  # make a dataframe for circular plotting
  Allchr_circular <- as.data.frame(Allchrs[,c(1,3,8,9,10,11,12,13)])
  # create end coordinates based on next start site and have NA for last start site
  Allchr_circular$end <- c(Allchr_circular$pos[-1]-1, NA)
  #View(Allchr_circular)
  # check for changes in chr at the last window
  temp <- ifelse(Allchr_circular$end != Allchr_circular$pos + 999999,
                                Allchr_circular$pos+999999,
                                Allchr_circular$end)
  # add this to the dataframe
  Allchr_circular$end <- temp
  # make an entry for last end positioin 
  Allchr_circular$end[nrow(Allchr_circular)]<-Allchr_circular$pos[nrow(Allchr_circular)]+999999
  #View(Allchr_circular)
  # Now fix the first entry of each chr
  Allchr_circular$end <- ifelse(Allchr_circular$pos < 999999,
                                999999,
                                Allchr_circular$end)
  
  # reorder the columns
  Allchr_circular <- Allchr_circular[, c(1,2,9,3,4,5,6,7,8)]
  names(Allchr_circular)[names(Allchr_circular) == "pos"] <- "start"
  # ok looks good
  
  
if(sample == sample_vector[1]) { 
  circos.initializeWithIdeogram(Allchr_circular,sort.chr = FALSE,
                                plotType = c("axis", "labels"),
                                tickLabelsStartFromZero = FALSE,
                                major.by = 50000000,
                                axis.labels.cex = 1,
                                labels.cex = 2)
}  
  
  circos.genomicTrackPlotRegion(data=Allchr_circular,panel.fun=function(region,value,...) {
    circos.genomicLines(region,value,type="l",
                        col=c("gray","blue","light blue","red","yellow","purple"),
                        area = TRUE,
                        border = c("gray","blue","light blue","red","yellow","purple"),
                        lwd=2)
                        }, track.height = 0.1)
  
  # circos.info(plot = TRUE) # this prints sector and track names
 # set.current.cell(sector.index = "chr01", track.index = 1)
  # label the track
  circos.text(0, 0.5, # label the middle of the track
              sub(".*hecki_", "", sample), # extract only sample ID  
              facing = "bending.inside",
              cex = 2,
              sector.index = "chr01", #next to chr01
              track.index = get.current.track.index(), 
              pos=2, offset = 0.05, # this makes the text offset to the left
              adj=1) # this makes the text right justified
  
  # label the species in the center
  text(0, 0, expression(italic("M. hecki")), cex = 4)
  
}
dev.off()
```
        

