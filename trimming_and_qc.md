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

Then I used paired reads from the scythe/trimmomatic/quake corrections only.  The argument for discarding single end reads is here (https://gatkforums.broadinstitute.org/gatk/discussion/2493/mixing-paired-end-and-single-end-reads); basically single end reads mess up the map qualities.
