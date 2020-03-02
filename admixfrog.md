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
vcftools --gzvcf ../FandM_chr01_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --max-missing 0.5 --minQ 30 --recode --recode-INFO-all --out ./FandM_chr01_mm_0.5_minQ_30.vcf.gz
```

The next filter we will apply is a minimum depth for a genotype call 
This command will recode genotypes that have less than 3 reads.
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.vcf.gz --minDP 3 --recode --recode-INFO-all --out FandM_chr01_mm_0.5_minQ_30_minDP_3.vcf.gz 
```

The next step is to get rid of individuals that did not sequence well. We can do this by assessing individual levels of missing data.
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30_minDP_3.vcf.gz --missing-indv
```

Examine the output file

```
cat out.imiss
```

Make a list of individuals with <50% missing data
```
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
```

If needed, we can remove individuals with lots of missing data like this:
```
vcftools --vcf raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out raw.g5mac3dplm
```

Now restrict the data to variants called in a high percentage of individuals from each species.
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
vcftools --vcf DP3g95maf05.recode.vcf --keep 1.keep --missing-site --out 1
vcftools --vcf DP3g95maf05.recode.vcf --keep 2.keep --missing-site --out 2 
vcftools --vcf DP3g95maf05.recode.vcf --keep 3.keep --missing-site --out 3
vcftools --vcf DP3g95maf05.recode.vcf --keep 4.keep --missing-site --out 4 
vcftools --vcf DP3g95maf05.recode.vcf --keep 5.keep --missing-site --out 5
vcftools --vcf DP3g95maf05.recode.vcf --keep 6.keep --missing-site --out 6 
vcftools --vcf DP3g95maf05.recode.vcf --keep 7.keep --missing-site --out 7
vcftools --vcf DP3g95maf05.recode.vcf --keep 8.keep --missing-site --out 8 
```
Now make a list of the bad loci we want to filter:
```
cat 1.lmiss 2.lmiss 3.lmiss 4.lmiss 5.lmiss 6.lmiss 7.lmiss 8.lmiss | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci
```
Now filter these loci:
```
vcftools --vcf DP3g95maf05.recode.vcf --exclude-positions badloci --recode --recode-INFO-all --out DP3g95p5maf05
```


