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

Begin by making a file with information from the vcf that can be fed into R for plotting:

```
java -jar GenomeAnalysisTK.jar \
     -R reference.fasta
     -T VariantsToTable \
     -V ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf.gz \
     -F CHROM -F POS -F ID -F QUAL -F MQ -F MQRankSum -F DP -F ReadPosRankSum -F SOR \
     -o results.table
```
example:

```bash
sqsub -r 1d --mpp 6G -o temp.log /usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java  -Xmx2g -jar ../bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/GenomeAnalysisTK.jar -R ../MacaM/MacaM_mt_y.fa -T VariantsToTable -V ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf.gz -F CHROM -F POS -F ID -F QUAL -F MQ -F MQRankSum -F DP -F ReadPosRankSum -F SOR -o results.table
```
