# Workflow

I am going to trim with scythe, then trimmomatic, and then quake, assuming a cutoff kmer of 1 for all genomes.  This should generate some very nice high quality data.  Only paired reads will be used for mapping because I read somewhere that using single end reads can alter the map qualities in an undesirable way.

# Initial QC

FastQC gives us an initial glimpse at data quality:
```
#!/usr/bin/perl
# This script will use trimmomatic to trim all of the fastq reads

my $status;
my @files;

@files = glob("../SEAsian_macaques_rawdata/hecki*/*.fastq.gz");

# Take off the ".fastq.gz" suffix
foreach(@files){
    $_=substr($_, 0, -9);
}

foreach(@files){
    my $commandline = "../bin/scythe/scythe -a ../bin/scythe/adap_TruSeqonly.fa -p 0.1 -q sanger ".$_.".fastq.gz | gzip > ".$_."scythe.fastq.gz";
    print $commandline,"\n\n";
    $status = system($commandline);
}
ben_evans@bionc03:/mnt/expressions/ben_evans/ben_scripts$ emacs 0.7_scythe_nigra.pl
ben_evans@bionc03:/mnt/expressions/ben_evans/ben_scripts$ emacs 0.7_scythe_nem.pl
ben_evans@bionc03:/mnt/expressions/ben_evans/ben_scripts$ emacs 0.7_scythe_tonk.pl
ben_evans@bionc03:/mnt/expressions/ben_evans/ben_scripts$ emacs 0.7_scythe_maura.pl
ben_evans@bionc03:/mnt/expressions/ben_evans/ben_scripts$ emacs 0.7_scythe_hecki.pl
ben_evans@bionc03:/mnt/expressions/ben_evans/ben_scripts$ ls
0.7_scythe_hecki.pl   0.7_scythe_nem.pl~    0_fastqc.pl			      1_trimmomatic_hecki.pl~  1_trimmomatic_nigra.pl	1_trimmomatic_tonk.pl~	temp.fasta
0.7_scythe_hecki.pl~  0.7_scythe_nigra.pl   0_fastqc.pl~		      1_trimmomatic_maura.pl   1_trimmomatic_nigra.pl~	2_quake_new.pl~
0.7_scythe_maura.pl   0.7_scythe_nigra.pl~  1.5_fastqc_after_trimmomatic.pl   1_trimmomatic_maura.pl~  1_trimmomatic.pl		2_quake_OLD.pl
0.7_scythe_maura.pl~  0.7_scythe_tonk.pl    1.5_fastqc_after_trimmomatic.pl~  1_trimmomatic_nem.pl     1_trimmomatic.pl~	2_quake.pl
0.7_scythe_nem.pl     0.7_scythe_tonk.pl~   1_trimmomatic_hecki.pl	      1_trimmomatic_nem.pl~    1_trimmomatic_tonk.pl	2_quake.pl~
ben_evans@bionc03:/mnt/expressions/ben_evans/ben_scripts$ more 0_fastqc.pl
#!/usr/bin/perl
# This script will use fastqc to do qc on raw reads

my $status;
my @files;

@files = glob("../SEAsian_macaques_rawdata/*/*.fastq.gz");

foreach(@files){
    my $commandline = "/mnt/expressions/ben_evans/bin/FastQC/fastqc ".$_;
    print $commandline,"\n";
    $status = system($commandline);
}
```

One of the findings from this analysis is that M. nemestrina individual PM1206 has far fewer reads than the other individuals.  This is a concern that I may address with more sequencing.

# Trimming 3' ends with Scythe

I plan to use three appreachws to trim the data: scythe, which trims of 3' adapter sequences, trimmomatic, which trims 5' and 3' sequences, and quake, which is a kmer based approach designed to get rid of rare reads that are mostly sequencing errors.

First scythe with this script: 

0.7_scythe_hecki.pl:

```
#!/usr/bin/perl
# This script will use scythe to trim 3' adapter seqs

my $status;
my @files;

@files = glob("../SEAsian_macaques_rawdata/hecki*/*.fastq.gz");

# Take off the ".fastq.gz" suffix
foreach(@files){
    $_=substr($_, 0, -9);
}

foreach(@files){
    my $commandline = "../bin/scythe/scythe -a ../bin/scythe/adap_TruSeqonly.fa -p 0.1 -q sanger ".$_.".fastq.gz | gzip > ".$_."scythe.fastq.gz";
    print $commandline,"\n\n";
    $status = system($commandline);
}
```

After scythe, I trimmed with trimmomatic and then with quake eliminating for all reads that had a kmer with a frequency of 1.

Trimmomatic like this:
```
java -jar ../bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog _log.txt /mnt/scratch/ben_evans/rhesus_macaque_genomez_original_raw_data/scythe_only/SRR1952145_1scythe.fastq.gz /mnt/scratch/ben_evans/rhesus_macaque_genomez_original_raw_data/scythe_only/SRR1952145_2scythe.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_1scythe_and_trimm_paired.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_1scythe_and_trimm_single.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_2scythe_and_trimm_paired.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_2scythe_and_trimm_single.fastq.gz ILLUMINACLIP:../bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```
and this
```
#!/usr/bin/perl                                                                                                                                                   
# This script will use trimmomatic to trim all of the fastq reads                                                                                                 

my $status;
my @files;

#@files = glob("../SEAsian_macaques_rawdata/*/*.fastq.gz");                                                                                                       

@files = glob("/mnt/scratch/ben_evans/rhesus_macaque_genomez_original_raw_data/scythe_only/SRR1952145*scythe.fastq.gz");

foreach(@files){
    $_=substr($_, 0, -9);
}


my $commandline = "java -jar ../bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog ".$_."_log.txt ";
foreach(@files){
    $commandline = $commandline.$_.".fastq.gz ";
}
foreach(@files){
    $commandline = $commandline.$_."_and_trimm_paired.fq.gz ".$_."_and_trimm_single.fq.gz ";
}
$commandline = $commandline."ILLUMINACLIP:../bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36";

print $commandline,"\n";
#$status = system($commandline);                                                                                                         ```
Two were truncated so I redid them like this:
```java -jar ../bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog _log.txt /mnt/scratch/ben_evans/rhesus_macaque_genomez_original_raw_data/scythe_only/SRR1952145_1scythe.fastq.gz /mnt/scratch/ben_evans/rhesus_macaque_genomez_original_raw_data/scythe_only/SRR1952145_2scythe.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_1scythe_and_trimm_paired.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_1scythe_and_trimm_single.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_2scythe_and_trimm_paired.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SRR1952145_2scythe_and_trimm_single.fastq.gz ILLUMINACLIP:../bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36


java -jar ../bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog _log.txt /mnt/scratch/ben_evans/rhesus_macaque_genomez_original_raw_data/scythe_only/SRR1929379_1scythe.fastq.gz /mnt/scratch/ben_evans/rhesus_macaque_genomez_original_raw_data/scythe_only/SRR1929379_2scythe.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SAMN03264624/SRR1929379_1scythe_and_trimm_paired.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SAMN03264624/SRR1929379_1scythe_and_trimm_single.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SAMN03264624/SRR1929379_1scythe_and_trimm_paired.fastq.gz /mnt/expressions/ben_evans/rhesus_macaque_genomez/SAMN03264624/SRR1929379_1scythe_and_trimm_single.fastq.gz ILLUMINACLIP:../bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```




Then I used paired reads from the scythe/trimmomatic/quake corrections only.  The argument for discarding single end reads is here (https://gatkforums.broadinstitute.org/gatk/discussion/2493/mixing-paired-end-and-single-end-reads); basically single end reads mess up the map qualities.

# Quake

I'm using quake to reduce the number of sequencing errors, fixing for all genomes a kmer of 1.

Here's an example of the commandline on info:
```
/usr/local/quake/src/correct -f filenamez.txt -z -k 19 -c 1 -m jelly_dump_all_19mers -p 4 -q 33
```
Important to remember is that the file 'filename' needs to include the paired reads on the same line separated by whitespace.  This is important in order to preserve the paired reads.  If this is forgotten, they can be repaired with bbmap.



