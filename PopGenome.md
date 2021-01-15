# Testing for adaptive introgression

An interesting thing to do with the introgression tracks is to ask whether they have any evidence of adaptive evolution. One way to do this is to compare dN/dS in the introgressed tracks to the rest of the genome.  Also similar with polymorphism stats.  An R package called PopGenome may allow me to do this: https://cran.r-project.org/web/packages/PopGenome/vignettes/Whole_genome_analyses_using_VCF_files.pdf

# extract specific chr from ref file
```
awk '/chr19/,/chrM/' MacaM_mt_female.fa > MacaM_mt_female_chr19only.fa
```
* still need to remove last line though...


# Work with filtered SNPs only vcf files:

directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data
```
example file (SNPs only):
```
FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz
```

# R code (in progress)

```R
## Working directory
# https://cran.r-project.org/web/packages/PopGenome/vignettes/Whole_genome_analyses_using_VCF_files.pdf
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/PopGenome")
library(ggplot2)
library(PopGenome)

#bru_PF707
#hecki_PF505
#hecki_PF643
#hecki_PF644
#hecki_PF647
#hecki_PF648
#maura_PF615
#maura_PF713
#maura_PM613
#maura_PM614
#maura_PM616
#nem_GumGum_female
#nem_Ngsang_sumatra_female
#nem_PM1206
#nem_PM664
#nem_PM665
#nem_Sukai_male
#nigra_PF1001
#nigra_PF660
#nigra_PM1003
#nigrescens_PM1011
#nigrescens_PM654
#tog_PF549
#tonk_PF511
#tonk_PF559
#tonk_PF563
#tonk_PF597
#tonk_PF626
#tonk_PM592

# read the vcf file
GENOME.class <- readVCF("FandM_chr19_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz",
                        numcols=10000, 
                        tid="chr19",
                        from=1, 
                        to= 53113586, # for chr19
                        approx=FALSE, 
                        out="", 
                        parallel=FALSE, 
                        gffpath="MacaM_Rhesus_Genome_Annotation_v7.6.8.gff")

# define the populations
GENOME.class <- set.populations(GENOME.class,
                                list(
                                c("bru_PF707"),
                                c("hecki_PF505","hecki_PF643","hecki_PF644","hecki_PF647","hecki_PF648"),
                                c("maura_PF615","maura_PF713","maura_PM613","maura_PM614","maura_PM616"),
                                c("nem_Ngsang_sumatra_female"),
                                c("nem_GumGum_female","nem_PM1206","nem_PM664","nem_PM665","nem_Sukai_male"),
                                c("nigra_PF1001","nigra_PF660","nigra_PM1003"),
                                c("nigrescens_PM1011","nigrescens_PM654"),
                                c("tog_PF549"),
                                c("tonk_PF511","tonk_PF559","tonk_PF563","tonk_PF597","tonk_PF626","tonk_PM592")), 
                                diploid=TRUE)

# verify the SNPs (there is a typo on pg 8 of the manual about this command)
GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="MacaM_mt_female_chr19only.fa",save.codons=FALSE)

GENOME.class@region.data@synonymous
GENOME.class@region.data@CodingSNPS

# This is complex object, with several slots
get.sum.data(GENOME.class)
show.slots(GENOME.class)

# this works but takes forever
get_gff_info(object=GENOME.class,"MacaM_Rhesus_Genome_Annotation_v7.6.8.gff","chr19",
             position=1,feature=FALSE,extract.gene.names=FALSE)


split_data_into_GFF_attributes(GENOME.class,
                               "MacaM_Rhesus_Genome_Annotation_v7.6.8.gff", 
                               "chr19", 
                               "GID")
GENOME.class.split@region.names
GENOME.class.split@feature.names



# https://wurmlab.github.io/genomicscourse/2016-SIB/practicals/population_genetics/popgen
# Diversities and FST (by scaffold)
GENOME.class <- F_ST.stats(GENOME.class) # this does the calculations and 
                                          # adds the results to the appropriate slots

# Print FST
get.F_ST(GENOME.class) # each line is a scaffold
GENOME.class@nucleotide.F_ST

# Transform object into object divided by sliding window
win_snp <- sliding.window.transform(GENOME.class, 
                                    width=10000, jump=10000, 
                                    type=2,
                                    whole.data=FALSE)

# Measurements per window
win_snp <- F_ST.stats(win_snp)

win_snp@nucleotide.F_ST
win_snp@nuc.diversity.within

# A simple plot
win_fst <- win_snp@nucleotide.F_ST[,1]
bb_div  <- win_snp@nuc.diversity.within[,1] # diversity among B (bb = "big B")
lb_div  <- win_snp@nuc.diversity.within[,2] # diversity among B (lb = "little b")


plot(1:length(win_fst), win_fst)

par(mar=c(1,1,1,1))
win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(bb_div), bb_div)

win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(lb_div), lb_div)

GENOME.class@Coding.region

genes <- splitting.data(GENOME.class, subsites="coding")
is(genes)

```
