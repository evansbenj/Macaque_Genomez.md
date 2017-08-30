# Trimming

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
