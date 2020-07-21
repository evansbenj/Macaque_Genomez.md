# Admixfrog

Vcf files are here on graham:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM*vcf.gz
```

On graham I installed admixfrog here `/home/ben/.local/bin/admixfrog` as follows:
```
module load scipy-stack/2019b
pip install cython scipy --upgrade --user
pip install git+https://github.com/benjaminpeter/admixfrog --user
```

Now I think it works based on this:
```
/home/ben/.local/bin/admixfrog --help
```


# Filtering

I'm going to agressively filter the data before analysis using admix frog as follows:

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

This shows that, for chr01, nigrescens_PM654, hecki_PF644, hecki_PF647, tonk_PF597 have a proportion of missing sites respectively of 0.369508, 0.389771, 0.435798, 0.548834). All others have 10% (tonk_PF559) or <6%.

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
```
```
mawk '$2 == "brunnescens"' species.txt > 1.keep && mawk '$2 == "hecki"' species.txt > 2.keep && mawk '$2 == "maura"' species.txt > 3.keep && mawk '$2 == "nemestrina"' species.txt > 4.keep && mawk '$2 == "nigra"' species.txt > 5.keep && mawk '$2 == "nigrescens"' species.txt > 6.keep && mawk '$2 == "togeanus"' species.txt > 7.keep && mawk '$2 == "tonkeana"' species.txt > 8.keep 
```
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 1.keep --missing-site --out 1_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 2.keep --missing-site --out 2_chr01 
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 3.keep --missing-site --out 3_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 4.keep --missing-site --out 4_chr01 
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 5.keep --missing-site --out 5_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 6.keep --missing-site --out 6_chr01 
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 7.keep --missing-site --out 7_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 8.keep --missing-site --out 8_chr01 
```
Now make a list of the bad loci we want to filter:
```
cat 1_chr01.lmiss 2_chr01.lmiss 3_chr01.lmiss 4_chr01.lmiss 5_chr01.lmiss 6_chr01.lmiss 7_chr01.lmiss 8_chr01.lmiss | mawk '!/CHR/' | mawk '$6 > 0.2' | cut -f1,2 >> bad_chr01_loci
```
Now filter these loci:
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --exclude-positions bad_chr01_loci --recode --recode-INFO-all --out FandM_chr01_mm_0.5_minQ_30_exclude_missingness
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

# Admixfrog

OK first load baftools and index the filtered file
```
module load bcftools/1.9
bgzip -c FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.recode.vcf > FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
bcftools index FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
tabix -p vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
```


# Admixfrog seems to be working
I am in this directory on graham:
`/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/filtered_for_admixfrog`

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

BRU
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_BRU --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states BRU=bru_PF707 --state-file pops.yaml     
```
HEC
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_HEC1 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states HEC1=hecki_PF505 --state-file pops.yaml
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_HEC2 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states HEC2=hecki_PF643 --state-file pops.yaml
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_HEC3 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states HEC3=hecki_PF644 --state-file pops.yaml
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_HEC4 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states HEC4=hecki_PF647 --state-file pops.yaml   
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_HEC5 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states HEC5=hecki_PF648 --state-file pops.yaml  
```
NEM
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_BNEM1 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states NEM1=nem_GumGum_female --state-file pops.yaml 
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_BNEM2 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states NEM2=nem_PM1206 --state-file pops.yaml  
```
``` 
/home/ben/.local/bin/admixfrog-ref --outfile ref_BNEM3 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states NEM3=nem_PM664 --state-file pops.yaml  
```
```
 /home/ben/.local/bin/admixfrog-ref --outfile ref_BNEM4 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states NEM4=nem_PM665 --state-file pops.yaml  
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_BNEM5 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states NEM5=nem_Sukai_male --state-file pops.yaml 
```

TON:

```
/home/ben/.local/bin/admixfrog-ref --outfile ref_TON1 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states TON1=tonk_PF511 --state-file pops.yaml 
```  
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_TON2 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states TON2=tonk_PF559 --state-file pops.yaml 
``` 
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_TON3 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states TON3=tonk_PF563 --state-file pops.yaml 
``` 
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_TON4 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states TON4=tonk_PF597 --state-file pops.yaml 
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_TON5 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states TON5=tonk_PF626 --state-file pops.yaml 
```
```
/home/ben/.local/bin/admixfrog-ref --outfile ref_TON6 --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --states TON6=tonk_PM592 --state-file pops.yaml 
```

make the target (input) file
```
admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF511sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF511.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF559sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF559.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF563sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF563.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF597sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF597.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF626sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF626.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PM592.in.xz
```
run the analysis
```
admixfrog --infile tonk_PF511.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF511_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF559.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF559_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF563.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF563_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF597.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF597_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF626.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF626_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PM592.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PM592_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination

```


randomstuff
```
        --states AFR BNEM=nem_GumGum_female,nem_Ngsang_sumatra_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male DEN=Denisova \
        --pop-file data.yaml \
        (no need for this until we have separate files by chr):
        --rec-file rec.{CHROM}
        
--states MA=MA1+MA2+MA3+MA4+MA5 MB=MB1+MB2
        

