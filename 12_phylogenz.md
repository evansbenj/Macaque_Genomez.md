# Making phylogenz

I'm using iqtree again to make concatenated treez.  I'm gonna do this for mtDNA and yDNA using genotypes called based on depth.  Then I'll do this for aDNA and xDNA using genotypes with degenerate bases (or for xDNA, with degen bases for females only, maybe).  I generated the mtDNA and yDNA depth files using a script that is in the mtDNA page.

```
/work/ben/2017_SEAsian_macaques/bin/iqtree-1.5.0a-Linux/bin/iqtree -s ../SEAsian_macaques_bam/females_and_males/FandM_chrM_BSQR_jointgeno_allsites_filtered.vcf.gz_bydepth.nxs -m TEST -nt 1 -pre ../SEAsian_macaques_bam/females_and_males/FandM_chrM_BSQR_jointgeno_allsites_filtered.vcf.gz_bydepth.nxs_
```

based on `.iqtree` output from iqtree, HKY+I+G4 is the favored model.  Now do ultrafast bootstraps:

```
/work/ben/2017_SEAsian_macaques/bin/iqtree-1.5.0a-Linux/bin/iqtree -s ../SEAsian_macaques_bam/females_and_males/FandM_chrM_BSQR_jointgeno_allsites_filtered.vcf.gz_bydepth.nxs -m HKY+I+G4 -bb 1000
```
