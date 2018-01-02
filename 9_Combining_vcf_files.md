# Combining VCF files

Sometimes the genotyper does not complete an entire chr in the allotted time.  So we can do the rest in an interval and then combine the sections we need.

For the truncated file, we first need to remove the last line, which is prematurely terminated.
```
gzip -cd "$files" | sed -e '$d' | gzip > "$files".tmp
```
repeat this to ensure last line removal with no error, rename and uncompress to avoid issues:
```
gunzip myvcf.vcf.gz
```
then compress again
```
/work/ben/2017_SEAsian_macaques/bin/htslib-1.6/bin/bgzip myvcf.vcf
```
then use tabix to create a tbi index
```
/work/ben/2017_SEAsian_macaques/bin/htslib-1.6/bin/tabix -p vcf myvcf.vcf.gz
```

For small bit at the end, extract only the section needed using gatk. This section should begin at the line after the last line from the previous file.
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -jar /work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/GenomeAnalysisTK.jar -T SelectVariants -R /work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_y.fa --variant nem_PM665sorted_rg_realign_dedup.bam_chr14_155000000_169736342_noBSQR.g.vcf.gz -L chr14:159035771-169736342 -o nem_PM665sorted_rg_realign_dedup.bam_chr14_159035771_169736342_noBSQR.g.vcf.gz
```
then use cat variants to combine
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -cp /work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R /work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_y.fa -V nem_PM665sorted_rg_realign_dedup.bam_1_159035770_chr14_noBSQR.g.vcf.gz -V nem_PM665sorted_rg_realign_dedup.bam_chr14_159035771_169736342_noBSQR.g.vcf.gz -out nem_PM665sorted_rg_realign_dedup.bam_chr14_noBSQR.g.vcf.gz -assumeSorted
```
Done!
