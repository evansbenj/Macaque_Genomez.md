# Workflow

I plan to use bwa to map the reads to the MacaM genome sequence.  This genome assembly supposedly is much better than rhemac2 and it has a good annotation file.  The drawback is that it does not have the pairwise alignments to the other genomes that the UCSC genome browser provides. And I won't have the repeat masker annotation file that rhemac2 has. I can look into trying to generate this.

After mapping and sorting with bwa version 0.7.16a-r1185-dirty, I plan to dedup with samtools 1.5 (using htslib 1.5), dedup with picardtools version 2.12.1, and possibly to filter using bcftools.  In this way I could avoid GATK completely, which is an annoying package to work with.  The only drawback is that I'd not use base recalibration, but for high converage genomes, this shouldn't be too big of a deal. Added to this, I do not have a know set of SNPs for the SE Asian macaques, although one does exist for rhesus.

Filtering should be agressive and potentially include (1) repetitive regions, (2) regions that depart from HWE, (3) regions near indels, (4) regions with very high or very low coverage.

Additionally I would like to use a HMM to identify and screen out heterozygous positions on the male chrX.


# Toy data for testing

```
zcat PF505_all_R2scythe_and_trimm_paired.cor.fastq.gz | head -n 1000 > toyR2.fastq 
```
problems with misnamed reads fixed with bbmap like this (from within each directory):
```
bbmap/bbmap/repair.sh in1=toyR1.fastq in2=toyR2.fastq out1=toy1fixed1.fq out2=toy2fixed2.fq outsingle=toysingle.fq
```

# Running Java 8 on sharcnet (no problem on bionc)

```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -version
```

# Dedup with picard

on iqaluk
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -jar /project0/ben/bin/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=maura_PM613sorted.bam O=maura_PM613sorted_dedup.bam  M=PM613marked_dup_metrics.txt TMP_DIR=`pwd`/TMP
```
on bionc03, from within '/mnt/scratch/ben_evans/SEAsian_bam_files'
```
qsub -l h_vmem=120g -cwd -b y -N PM613sorted bash -c "java -jar /mnt/expressions/ben_evans/bin/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=maura_PM613sorted.bam O=maura_PM613sorted_dedup.bam M=PM613marked_dup_metrics.txt"
```


```
