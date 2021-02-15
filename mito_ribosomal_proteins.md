# Mitoribosomal proteins

I identified 71 annotated MRP genes:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep '=MRP' | egrep 'transcript_01' | wc -l
```
Note that we are not including URGCP-MRPS24 because this is a readthrough transcript that does not encode functional MRPS24 (which is a separate gene that is in our list.

I got all non-MRP genes by adding a `-v` invert flag to egrep:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep -v '=MRP' > temp.txt
```
