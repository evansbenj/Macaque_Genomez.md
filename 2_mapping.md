# Workflow

I plan to use bwa to map the reads to the MacaM genome sequence.  This genome assembly supposedly is much better than rhemac2 and it has a good annotation file.  The drawback is that it does not have the pairwise alignments to the other genomes that the UCSC genome browser provides. And I won't have the repeat masker annotation file that rhemac2 has. I can look into trying to generate this.

After mapping and sorting with bwa version 0.7.16a-r1185-dirty, I plan to dedup with samtools 1.5 (using htslib 1.5), dedup with picardtools version 2.12.1, and possibly to filter using bcftools.  In this way I could avoid GATK completely, which is an annoying package to work with.  The only drawback is that I'd not use base recalibration, but for high converage genomes, this shouldn't be too big of a deal. Added to this, I do not have a know set of SNPs for the SE Asian macaques, although one does exist for rhesus.

Filtering should be agressive and potentially include (1) repetitive regions, (2) regions that depart from HWE, (3) regions near indels, (4) regions with very high or very low coverage.

Additionally I would like to use a HMM to identify and screen out heterozygous positions on the male chrX.

# If problems, validate fastq files with fastQValidator:

on Bionc

```
/mnt/expressions/ben_evans/bin/fastQValidator/bin/fastQValidator --file SRR1951019_1scythe_and_trimm_paired.fastq.gz
```

# Repair with bbmap
On iqaluk, needed to repair unpaired reads from quake because I did not include the paired fastq files on the same line with a whitespace separator.

```
/project0/ben/bin/bbmap/bbmap/repair.sh -Xmx30g in1=/project0/ben/SEAsian_scy_trim_quake/nem_PM1206/nem_PM1206_all_R1scythe_and_trimm_paired.cor.fastq.gz in2=/project0/ben/SEAsian_scy_trim_quake/nem_PM1206/nem_PM1206_all_R2scythe_and_trimm_paired.cor.fastq.gz out1=/project0/ben/SEAsian_scy_trim_quake/nem_PM1206/nem_PM1206_all_R1scythe_and_trimm_paired.corfixed.fastq.gz out2=/project0/ben/SEAsian_scy_trim_quake/nem_PM1206/nem_PM1206_all_R2scythe_and_trimm_paired.corfixed.fastq.gz outsingle=/project0/ben/SEAsian_scy_trim_quake/nem_PM1206/nem_PM1206_all_R1scythe_and_trimm_pairsingletons.fq.gz 
```

# Mapping with bwa

For male genomes the reference will be MacaM_mt_y.fa, which includes the chrY.  For females the reference will be MacaM_mt_female.fa, which does not include chrY. Both have mtDNA.  Here is an example with a female, with a pipe to samtools for sorting the bam.

```
/project0/ben/bin/bwa/bwa mem /project0/ben/MacaM/MacaM_mt_female.fa /project0/ben/SEAsian_scy_trim_quake/hecki_PF644/PF644_all_R1scythe_and_trimm_paired.corfixed.fq.gz /project0/ben/SEAsian_scy_trim_quake/hecki_PF644/PF644_all_R2scythe_and_trimm_paired.corfixed.fq.gz -t 4 | /project0/ben/bin/samtools-1.5/samtools view -Shu - | /project0/ben/bin/samtools-1.5/samtools sort - -o /project0/ben/SEAsian_macaques_bam/hecki_PF644sorted.bam
```


# Running Java 8 on sharcnet (no problem on bionc)

```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -version
```


# Add readgroups

This should have been done with bwa mem, but I can do it now with Picard AddOrReplaceReadGroups (http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) 

'@RG\tID:foo\tSM:bar\tLB:library1'

```
java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
      O=output.bam \
      RGID=1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=SAMPLE_NAME
```

on iqaluk
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -jar /project0/ben/bin/picard/picard.jar AddOrReplaceReadGroups I=hecki_PF505sorted_dedup.bam O=hecki_PF505sorted_dedup_rg.bam RGID=hecki_PF505 RGLB=hecki_PF505 RGPL=illumina RGPU=hecki_PF505 RGSM=hecki_PF505
```

# Checking bam files

on bionc
```
/mnt/expressions/ben_evans/bin/samtools/bin/samtools quickcheck tonk_PF511.bam  && echo 'all ok' || echo 'fail!'
```

# Index bamfiles
on iqaluk
```
/project0/ben/bin/samtools-1.5/samtools index hecki_PF505sorted_dedup_rg.bam
```

# Indel realigner with GATK

(this pipeline is adapted from here http://www.htslib.org/workflow/)

```
java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R <ref.fa> -I <lane.bam> -o <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf>
```
```
java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R <ref.fa> -I <lane.bam> -targetIntervals <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf> -o <lane_realigned.bam>
```

on iqaluk (need to change to different ref for males)
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8g -jar /project0/ben/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /project0/ben/MacaM/MacaM_mt_female.fa -I hecki_PF648sorted_dedup_rg.bam -o forIndelRealigner_hecki_PF648.intervals
```
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8g -jar /project0/ben/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/GenomeAnalysisTK.jar -T IndelRealigner -R /project0/ben/MacaM/MacaM_mt_female.fa -I hecki_PF505sorted_dedup_rg.bam -targetIntervals forIndelRealigner_hecki_PF505.intervals -o hecki_PF505sorted_dedup_rg_realigned.bam
```

# index the realigned bams

```
samtools index <sample.bam>
```
Actually this is not needed because GATK does this on the fly

# Dedup with picard 

on iqaluk
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -jar /project0/ben/bin/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=hecki_PF505sorted_dedup_rg_realigned.bam O=hecki_PF505sorted_ddedup_rg_realigned.bam  M=hecki_PF505marked_ddup_metrics.txt TMP_DIR=`pwd`/TMP
```
on bionc03, from within '/mnt/scratch/ben_evans/SEAsian_bam_files'
```
qsub -l h_vmem=120g -cwd -b y -N PM613sorted bash -c "java -jar /mnt/expressions/ben_evans/bin/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=hecki_PF505sorted_dedup_rg_realigned.bam O=hecki_PF505sorted_ddedup_rg_realigned.bam M=PM613marked_ddup_metrics.txt"
```
on orca in directory `/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males`:
```
sqsub -r 1d --mpp 6G -o maura_PM613_dedup /usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx2g -jar /work/ben/2017_SEAsian_macaques/bin/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=maura_PM613sorted_dedup_rg_realigned.bam O=maura_PM613sorted_ddedup_rg_realigned.bam M=maura_PM613marked_ddup_metrics.txt

```


# BSQR 
## Call Variants
First I generated a list of known sites first for each sample for each chromosome. (4_gatk_haplotypecaller_sqsub_females.pl).  I did one for males using a reference genome with the Y chromosome.

```
#!/usr/bin/perl
# This script will run quake on trimmed fq files

my $gatkpath = "/work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/";
#my $referencegenome="/scratch/ben/MacaM/MacaM_mt_y.fa";
my $referencegenome="/work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_female.fa";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/";
my @chromosomes =("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrM");

@files = glob($majorpath."*ddedup_rg_realigned.bam");

foreach my $file (@files){
    if($file =~ /hecki_PF505/){
    foreach my $chromosome (@chromosomes){
	my $commandline = "sqsub -r 4d --mpp 16G -o ".$file."\.log ";
	$commandline = $commandline."/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -jar ".$gatkpath."GenomeAnalysisTK.jar -T HaplotypeCaller -R ".$referencegenome." -I ".$file." -L ".$chromosome." --emitRefConfidence GVCF -o ".$file."_".$chromosome."_noBSQR.g.vcf.gz";
# -out_mode EMIT_ALL_CONFIDENT_SITES
	print $commandline,"\n";
    $status = system($commandline);
    }
    }
}
```
## Cat Variants

Then I concatenated these variants for each individual to generate one SNP file for use in BQSR (4.5_gatk_catvariants_sqsub_females.pl):
```
#!/usr/bin/perl
# This script will run quake on trimmed fq files

my $gatkpath = "/work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/";
#my $referencegenome="/scratch/ben/MacaM/MacaM_mt_y.fa";
my $referencegenome="/work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_female.fa";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/";
my @chromosomes =("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrM");

@files = glob($majorpath."*ddedup_rg_realigned.bam");

foreach my $file (@files){
    my $commandline = "sqsub -r 2d --mpp 16G -o ".$file."\.log ";
    $commandline = $commandline."/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -cp ".$gatkpath."GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ".$referencegenome;
    foreach my $chromosome (@chromosomes){
	$commandline = $commandline." -V ".$file."_".$chromosome."_noBSQR.g.vcf.gz";
    }
    $commandline = $commandline." -out ".$file."_allchrs_noBSQR.g.vcf.gz";
    print $commandline,"\n";
#    $status = system($commandline);
}
```

## BQSR

Then I generated the BQSR table (5_gatk_BSQR_sqsub_females.pl) and used this table to print a new BQSR bam file (6_gatk_BSQR_printreads_sqsub_females.pl):
```
#!/usr/bin/perl
# This script will run quake on trimmed fq files

my $gatkpath = "/work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/";
#my $referencegenome="/scratch/ben/MacaM/MacaM_mt_y.fa";
my $referencegenome="/work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_female.fa";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/";
#my @chromosomes =("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrM");

@files = glob($majorpath."*ddedup_rg_realigned.bam");

foreach my $file (@files){
    if($file =~ /643/){
			my $commandline = "sqsub -r 4d --mpp 32G -o ".$file."_baserecalibrator\.log ";
			$commandline = $commandline."/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -jar ".$gatkpath."GenomeAnalysisTK.jar -T BaseRecalibrator -R ".$referencegenome." -knownSites ".$file."_allchrs_noBSQR.g.vcf.gz"." -I ".$file." -o ".$file."_lane_recal.table";
		    print $commandline,"\n";
#		    $status = system($commandline);
    }
}

```
```
#!/usr/bin/perl
# This script will print a new bam file with BSQR for females

my $gatkpath = "/work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/";
#my $referencegenome="/scratch/ben/MacaM/MacaM_mt_y.fa";
my $referencegenome="/work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_female.fa";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/";
my @chromosomes =("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrM");

@files = glob($majorpath."*ddedup_rg_realigned.bam");

foreach my $file (@files){
			my $commandline = "sqsub -r 5d --mpp 32G -o ".$file."_".$chromosome."\.log ";
			$commandline = $commandline."/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -jar ".$gatkpath."GenomeAnalysisTK.jar -T PrintReads -R ".$referencegenome." -I ".$file." -BQSR ".$file."_lane_recal.table -o ".$file."BSQR.bam";
		    print $commandline,"\n";
		    $status = system($commandline);
}

```



or use GATK:
```
java -Xmx32G -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R XXX.fa --variant XXX.g.vcf --variant YYY.vcf --variant ZZZ.g.vcf --includeNonVariantSites -o all_combined.vcf
```

Probably should break this up into chromosomes.


# Make some graphs of variants
```
bcftools stats -F <ref.fa> -s - <study.vcf.gz> > <study.vcf.gz.stats>
mkdir plots
plot-vcfstats -p plots/ <study.vcf.gz.stats>
```

# Filter

```
bcftools filter -O z -o <study_filtered..vcf.gz> -s LOWQUAL -i'%QUAL>10' <study.vcf.gz>
```


