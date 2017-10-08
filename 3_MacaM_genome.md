# MacaM genome

I decided to go with the MacaM genome because it fixes lots of errors with the rhemac2 genome and also seems better than the other two macaque genomes.  The problem is that, apart from a gene location gff file, there is not much else.  So minimally I will need to generate repeat masker annotations and also do whole genome alignments to outgroups such as baboons and humans.  Or maybe it is better to just align to rhemac2 and then use those files?  I'll try both.

# Index

```
/project0/ben/bin/samtools-1.5/samtools faidx /project0/ben/MacaM/MacaM_mt_female.fa
```
```
/project0/ben/bin/samtools-1.5/samtools faidx /project0/ben/MacaM/MacaM_mt_y.fa
```

# Create .dict file

```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -jar /project0/ben/bin/picard/picard.jar CreateSequenceDictionary REFERENCE=/project0/ben/MacaM/MacaM_mt_female.fa OUTPUT=/project0/ben/MacaM/MacaM_mt_female.dict
```
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -jar /project0/ben/bin/picard/picard.jar CreateSequenceDictionary REFERENCE=/project0/ben/MacaM/MacaM_mt_y.fa OUTPUT=/project0/ben/MacaM/MacaM_mt_y.dict
```

# Repeat masker

I installed RepeatMasker version open-4.0.7 on bionc and am running it now.  Repeat masker is set to use HMMER 3.1b2 (February 2015) to find the repeats.  I downloaded the most recent library of RepeatMaskerLib.embl, which is apparentaly "Dfam_Consensus RELEASE 20170127;   RepBase RELEASE 20170127".  But I'm actually not using this library and instead am using the default library that comes with Repeatmasker (Dfam 2.0).  The commandline was this:

```
./RepeatMasker -qq -species macaque -gff /mnt/expressions/ben_evans/MacaM/MacaM_mt_y.fa
```

The `-qq` flag tells it to do a quick search, `-gff` tells it to output a gff annotation file, and `species macaque` is a recognized species in the taxonomy and this will limit searches to repeats shared with humans and other OWMs.

Repeatmasker is on bionc in this directory:
```
/mnt/expressions/ben_evans/bin/RepeatMasker
```

gff file is on bionc here:
```
/mnt/expressions/ben_evans/MacaM/MacaM_mt_y.fa.out.gff
```


