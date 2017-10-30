# Plots

# Make some graphs of variants
```
bcftools stats -F <ref.fa> -s - <study.vcf.gz> > <study.vcf.gz.stats>
mkdir plots
plot-vcfstats -p plots/ <study.vcf.gz.stats>
```

example:

```
../bin/bcftools-1.6/bin/bcftools stats -F ../MacaM/MacaM_mt_female.fa -s - ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf.gz > ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf.gz.stats
```
```
module load intel/12.1.3
module load python/intel/2.7.8

```
```
../bin/bcftools-1.6/bin/plot-vcfstats -p plots/ ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf.gz.stats
```

# Density plots

Based on this site:
http://mbontrager.org/blog/2016/08/17/Variant-Exploration

but modified slightly because of problems with sed command to make this corrected pipe:

```bash
zcat  ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf | egrep -v "^#" | cut -f 8 | sed 's/^.*;DP=\([0-9]*\)*$/\1/' > ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf_depth.txt
```
