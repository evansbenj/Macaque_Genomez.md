# Sanity check and coverage

I'd like to make a PCA of all the data and also check the coverage of each sample. I plan to do a PCA using SNPRelate, maybe on individual chromosomes?

# Depth by bam files (less accurate because includes mtDNA)
```
#!/usr/bin/perl                                                                                                                             # This script will check depth of multiple files                                                                                           

my $samtoolspath="/work/ben/2017_SEAsian_macaques/bin/samtools-1.6/bin/";
my $gatkpath = "/work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/";
#my $referencegenome="/scratch/ben/MacaM/MacaM_mt_y.fa";                                                                                                                  
my $referencegenome="/work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_female.fa";
my $majorpath = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/*/";
my @chromosomes =("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr1\
8","chr19","chrX","chrM");

@files = glob($majorpath."*ddedup_rg_realigned.bamBSQR.bam");

foreach my $file (@files){
        my $commandline = "sqsub -r 4d --mpp 16G -o ".$file."\.log bash -c \"";
        $commandline = $commandline.$samtoolspath."samtools depth ".$file." |  awk \'{sum+=$3} END { print \"Depth of ".$file." = \",sum/NR}\' >> depth.txt\"";
        print $commandline,"\n";
    $status = system($commandline);
}

```

On iqaluk using screen:
```
../bin/samtools-1.6/samtools depth ../SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR.bam |  awk '{sum+=$3} END { print "Depth = ",sum/NR}' >> ../SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR.bam.depth.txt
```


# Depth by chr (better)

vcftools:
```
../bin/vcftools/bin/vcftools --gzvcf ../SEAsian_macaques_bam/females_and_males/FandM_chr19_BSQR_jointgeno_allsites_filtered.vcf.gz --depth --out ../SEAsian_macaques_bam/females_and_males/FandM_chr19_BSQR_jointgeno_allsites_filtered_depth
```
