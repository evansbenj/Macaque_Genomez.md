# Sanity check and coverage

I'd like to make a PCA of all the data and also check the coverage of each sample. I plan to do a PCA using SNPRelate, maybe on individual chromosomes?

# Depth
```
samtools depth  *bamfile*  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
```
