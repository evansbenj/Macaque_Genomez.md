# Admixfrog

Vcf files are here on graham:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas
```
These files should be used for ABBA_BABA, Fst, pi, etc for the genome.

For admixfrog, I need to filter these files so they can be processed in a manageable way.

On graham I installed admixfrog here `/home/ben/.local/bin/admixfrog` as follows:
```
module load scipy-stack/2019b
pip install cython scipy --upgrade --user
pip install git+https://github.com/benjaminpeter/admixfrog --user
```

Now I think it works based on this:
```
/home/ben/.local/bin/admixfrog --help
```


# Filtering

I'm going to agressively filter the data before analysis using admix frog as follows:

on graham, load vcftools:
```
module load nixpkgs/16.09
module load intel/2018.3
module load vcftools/0.1.14
```
Keep variants that have been successfully genotyped in 50% of individuals and a minimum quality score of 30
The --recode flag tells the program to write a new vcf file with the filters, 
The --recode-INFO-all keeps all the INFO flags from the old vcf file in the new one. 
```
vcftools --gzvcf ../FandM_chr01_BSQR_jointgeno_allsites_filtered_SNPsonly.vcf.gz --max-missing 0.5 --minQ 30 --recode --recode-INFO-all --out ./FandM_chr01_mm_0.5_minQ_30
```

The next step is to assess individual levels of missing data.
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30_minDP_3.vcf.gz --missing-indv
```

Examine the output file

```
cat out.imiss
```

This shows that, for chr01, papio, nigrescens_PM654, hecki_PF644, hecki_PF647, tonk_PF597 all have a proportion of missing sites respectively of 20.8%, 37.5%, 40.6%, 45.4%, 56.6%. All others have 10.6% (tonk_PF559) or <6.1%.

Now restrict the data to variants called in a high percentage of non-missing genotypes across individuals from every species.
Make a tab delimited file calleed 'species.txt' with the species assignment of each sample.
```
bru_PF707	brunnescens
hecki_PF505	hecki
hecki_PF643	hecki
hecki_PF644	hecki
hecki_PF647	hecki
hecki_PF648	hecki
maura_PF615	maura
maura_PF713	maura
maura_PM613	maura
maura_PM614	maura
maura_PM616	maura
nem_GumGum_female	nemestrina
nem_Ngsang_sumatra_female	nemestrina
nem_PM1206	nemestrina
nem_PM664	nemestrina
nem_PM665	nemestrina
nem_Sukai_male	nemestrina
nigra_PF1001	nigra
nigra_PF660	nigra
nigra_PM1003	nigra
nigrescens_PM1011	nigrescens
nigrescens_PM654	nigrescens
tog_PF549	togeanus
tonk_PF511	tonkeana
tonk_PF559	tonkeana
tonk_PF563	tonkeana
tonk_PF597	tonkeana
tonk_PF626	tonkeana
tonk_PM592	tonkeana
download  papio
```
```
mawk '$2 == "brunnescens"' species.txt > 1.keep && mawk '$2 == "hecki"' species.txt > 2.keep && mawk '$2 == "maura"' species.txt > 3.keep && mawk '$2 == "nemestrina"' species.txt > 4.keep && mawk '$2 == "nigra"' species.txt > 5.keep && mawk '$2 == "nigrescens"' species.txt > 6.keep && mawk '$2 == "togeanus"' species.txt > 7.keep && mawk '$2 == "tonkeana"' species.txt > 8.keep && mawk '$2 == "papio"' species.txt > 9.keep 
```
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 1.keep --missing-site --out 1_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 2.keep --missing-site --out 2_chr01 
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 3.keep --missing-site --out 3_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 4.keep --missing-site --out 4_chr01 
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 5.keep --missing-site --out 5_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 6.keep --missing-site --out 6_chr01 
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 7.keep --missing-site --out 7_chr01
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --keep 8.keep --missing-site --out 8_chr01 
```
Now make a list of the bad loci we want to filter:
```
cat 1_chr01.lmiss 2_chr01.lmiss 3_chr01.lmiss 4_chr01.lmiss 5_chr01.lmiss 6_chr01.lmiss 7_chr01.lmiss 8_chr01.lmiss | mawk '!/CHR/' | mawk '$6 > 0.2' | cut -f1,2 >> bad_chr01_loci
```
Now filter these loci:
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30.recode.vcf --exclude-positions bad_chr01_loci --recode --recode-INFO-all --out FandM_chr01_mm_0.5_minQ_30_exclude_missingness
```

After filtering, this reeduced SNPs by about half for chr01:
```
After filtering, kept 3178562 out of a possible 7229006 Sites
```
Now filter these loci (this level reduces the number of SNPs to 10% of the original):
```
vcftools --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness.recode.vcf --thin 500 --recode --recode-INFO-all --out FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned
```
now, for chr01: 
```
After filtering, kept 322613 out of a possible 3178562 Sites
```
I now have these in sbatch scrips: `vcftools_1.sh` and `vcftools_2.sh`

# Admixfrog

OK first load baftools and index the filtered file
```
module load bcftools/1.9
bgzip -c FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.recode.vcf > FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
bcftools index FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
tabix -p vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz
```


# Admixfrog seems to be working
I am in this directory on graham:
`/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/filtered_for_admixfrog`

Now create a reference file for each individual.  Before running admixfrog, do this:
```
module load scipy-stack/2019b
```
then for each sample do this:
```
admixfrog-ref [-h] --outfile OUTFILE [--states [STATES [STATES ...]]]
                     [--state-file STATE_FILE] [--cont-id CONT_ID]
                     [--ancestral ANCESTRAL]
                     [--random-read-samples [RANDOM_READ_SAMPLES [RANDOM_READ_SAMPLES ...]]]
                     [--vcf-ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz

```

This works
make the ref file
```
/home/ben/.local/bin/admixfrog-ref --vcf FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.vcf.gz --out FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --states TON HEC NEM --pop-file pops.yaml 
```
or just do this with `sbatch admixfrog_make_refs.sh chrX` (and modify the chr)

make the target (input) file
```
admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF511sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF511.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF559sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF559.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF563sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF563.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF597sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF597.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/tonk_PF626sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF626.in.xz
/home/ben/.local/bin/admixfrog-bam --bam /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/tonk_PM592sorted_ddedup_rg_realigned.bamBSQR.bam --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PM592.in.xz
```
or just do this with `sbatch admixfrog_make_target.sh chrX` (and modify the chr)

run the analysis
```
admixfrog --infile tonk_PF511.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF511_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF559.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF559_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF563.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF563_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF597.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF597_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PF626.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PF626_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination
admixfrog --infile tonk_PM592.in.xz --ref FandM_chr01_mm_0.5_minQ_30_exclude_missingness_thinned.ref.xz --out tonk_PM592_ref_TONK_HEC_NEM -b 10000 --states TON HEC NEM --c0 0 --dont-est-contamination

```

# Strategy for admixfrog
** within Sulawesi:
I think the best strategy is to focus on adjacent species when possible.  So for nigrescens - nigra and hecki; for hecki - tonk and nigresc; for tonk - hecki, maura; for maura - tonk and tog.

** betweeh Sulawesi and Borneo
I think the best strategy is to compare Borneo to Sumatra plus either hec, tonk, or mau.

# Problems overcome
Something was wrong with nem_PM1206 chr17 but I fixed it by changing the - 10000 flag to -b 20000 as follows:
```
admixfrog --infile nem_PM1206.chr17.in.xz --ref FandM_chr17_mm_0.5_minQ_30_exclude_missingness_thinned.recode.xz --out nem_PM1206_chr17_NEM_SUM_HEC.out -b 20000 --states NEM SUM HEC --c0 0 --dont-est-contamination
```


# Plotting circular plots

```R
## Working directory
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/admixfrog/TON_MAU_TOG")
library(tidyverse)
library(ggplot2)
rm(list=ls()) # removes all variables
sample <-"maura_PM616"
analysis <-"_TON_MAU_TOG"
chrs <- factor(c("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09",
         "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"))

# loop through chrs
for(i in levels(chrs)){
  print(eval(i))
  a <- read_csv(paste(eval(sample), "_",i,eval(analysis),".out.bin.xz", sep=""))
  assign(i,a)
}  

#for(i in levels(chrs)){
#  print(names(eval(as.name(i)))[1])
#  print(i[[names(eval(as.name(i)))[1]]])
 # assign(paste(i,"$",eval(names(eval(as.name(i)))[1],sep="")),444)
#}  

# rename chromosome column
chr01$chrom <- "chr01"
chr02a$chrom <- "chr02a"
chr02b$chrom <- "chr02b"
chr03$chrom <- "chr03"
chr04$chrom <- "chr04"
chr05$chrom <- "chr05"
chr06$chrom <- "chr06"
chr07$chrom <- "chr07"
chr08$chrom <- "chr08"
chr09$chrom <- "chr09"
chr10$chrom <- "chr10"
chr11$chrom <- "chr11"
chr12$chrom <- "chr12"
chr13$chrom <- "chr13"
chr14$chrom <- "chr14"
chr15$chrom <- "chr15"
chr16$chrom <- "chr16"
chr17$chrom <- "chr17"
chr18$chrom <- "chr18"
chr19$chrom <- "chr19"
chrX$chrom <- "chrX"

Allchrs <- rbind(chr01,chr02a,chr02b,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10,
                 chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX)

# make the chr a factor so we can use it for faceting
Allchrs$chrom <- as.factor(Allchrs$chrom)

#View(Allchrs)

Allchrs %>% gather(k, v, -chrom:-n_snps) %>% 
 # filter(v>.1) %>%
  ggplot(aes(x=map, y=v, fill=k)) + geom_col() + 
  facet_wrap(~chrom, ncol=1, strip='l')

ggsave(paste(eval(sample), "_",i,eval(analysis),".png", sep=""), width = 8, height = 8, dpi = 200)






##################
##################
# circular plotting
##################
##################

# make a dataframe for circular plotting
Allchr_circular <- as.data.frame(Allchrs[,c(1,3,8,9,10,11,12,13)])
# create end coordinates based on next start site and have NA for last start site
Allchr_circular$end <- c(Allchr_circular$pos[-1]-1, NA)
#View(Allchr_circular)
# check for changes in chr at the last window
temp <- ifelse(Allchr_circular$end != Allchr_circular$pos + 999999,
                              Allchr_circular$pos+999999,
                              Allchr_circular$end)
# add this to the dataframe
Allchr_circular$end <- temp
# make an entry for last end positioin 
Allchr_circular$end[nrow(Allchr_circular)]<-Allchr_circular$pos[nrow(Allchr_circular)]+999999
#View(Allchr_circular)
# Now fix the first entry of each chr
Allchr_circular$end <- ifelse(Allchr_circular$pos < 999999,
                              999999,
                              Allchr_circular$end)

# reorder the columns
Allchr_circular <- Allchr_circular[, c(1,2,9,3,4,5,6,7,8)]
names(Allchr_circular)[names(Allchr_circular) == "pos"] <- "start"
# ok looks good

# this is a list of chr lengths
chr_length_list = list("1" = 225002135,
                "2a" = 108967917,
                "2b" = 131519175,
                "3" = 198060209,
                "4" = 190369981,
                "5" = 179725205,
                "6" = 171866349,
                "7" = 185267708,
                "8" = 144034664,
                "9" = 111027318,
                "10" = 129655328,
                "11" = 127102482,
                "12" = 133317794,
                "13" = 95368959,
                "14" = 169736342,
                "15" = 92674614,
                "16" = 74750809,
                "17" = 76917410,
                "18" = 70128972,
                "19" = 53113586,
                "X" = 148935249
)

# this is a vector of chr lenths
begins <- rep(1,21)
ends = c(225002135,
                108967917,
                131519175,
                198060209,
                190369981,
                179725205,
                171866349,
                185267708,
                144034664,
                111027318,
                129655328,
                127102482,
                133317794,
                95368959,
                169736342,
                92674614,
                74750809,
                76917410,
                70128972,
                53113586,
                148935249)

# this is needed to initalize the graph below
chr_lengths <- cbind(begins,ends)

# chr names
chr_names <- c("chr1","chr2a","chr2b","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
              "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
               "chrX")
chr_names_ordered <- factor(chr_names, ordered = TRUE, 
                                levels = c("chr1","chr2a","chr2b","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                           "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                           "chrX"))

chr_names_ordered_simple <- factor(chr_names, ordered = TRUE, 
                            levels = c("1","2a","2b","3","4","5","6","7","8","9","10",
                                       "11","12","13","14","15","16","17","18","19",
                                       "X"))

# https://jokergoo.github.io/circlize_book/book/introduction.html
# Initialize library
library(circlize)
library(dplyr)
library(stringr)
# https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
# https://cran.r-project.org/web/packages/circlize/circlize.pdf
# initilize
circos.clear()

# set track size
circos.par("track.height" = 0.1)
# initializes the plot by setting the chr as a factor and 
# the position on each chr is defined as the x coord for each chr
# this is for plotting a histogram
circos.initialize(factors = chr_names_ordered, xlim = chr_lengths)

# This generates colored sectors 
circos.track(ylim = c(0, 1), bg.col = c(rep(c("light grey", "white"), 10),"grey40"))

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = c(rep(c("light grey", "white"), 10),"grey40"))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)


# this makes a nice scale
brk <- c(0,5,10,15,20,25)*10^7
col_text <- "grey40"
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",
              major.at=brk,
              labels=round(brk/10^7,1),
              sector.index = get.cell.meta.data("sector.index"),
              labels.cex=0.4,
              col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
},bg.border=F)







# trying ideogram (this works nicely!)
# https://mran.microsoft.com/snapshot/2014-12-11/web/packages/circlize/vignettes/genomic_plot.pdf
circos.clear()
# start at the top
# insert a space after chrX

circos.par("gap.degree" = c(rep(3, 20),10), 
            "cell.padding" = c(0, 0, 0, 0), 
            "start.degree" = 90)


circos.initializeWithIdeogram(Allchr_circular,sort.chr = FALSE,
                              plotType = c("axis", "labels"),
                              tickLabelsStartFromZero = FALSE,
                              major.by = 50000000,
                              axis.labels.cex = 0.3)


circos.genomicTrackPlotRegion(data=Allchr_circular,panel.fun=function(region,value,...) {
  circos.genomicLines(region,value,type="l",
                      col=c("gray","blue","light blue","red","yellow","green"),
                      area = TRUE,
                      border = c("gray","blue",NA,"red","yellow","green"),
                      lwd=2)
                      }, track.height = 0.1)

# circos.info(plot = TRUE) # this prints sector and track names
set.current.cell(sector.index = "chr01", track.index = 1)
# label the track
circos.text(1, -1.5,
            'PF660', 
            facing = "clockwise",
            cex = 0.7,
            pos = 3,
            offset = 4.5)

# label the species in the center
text(0, 0, "M. nigra", cex = 2.0)


ggsave(paste("nigra.png", sep=""), width = 8, height = 8, dpi = 200)



                                
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top")
})

circos.genomicAxis( )















highlight.sector(sector.index = c("chrX", "chr01"), track.index = 1, 
                 text = "fooo", col = NA)








circos.text(sector.index="chr1",track.index = 2,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "Sample_1",facing = "clockwise", 
            niceFacing = TRUE, adj = c(0,0),cex = 0.5)




circos.genomicTrackPlotRegion(Allchr_circular, numeric.column = 5,
                              panel.fun = function(region, value, ...) {
                                circos.genomicPoints(region, value, ...)
                                circos.genomicPoints(region, value)
                                # here `numeric.column` is measured in `value`
                                circos.genomicPoints(region, value, numeric.column = 1)
                              })

  circos.track(Allchr_circular,
               ylim = c(0, 1), 
               panel.fun = function(x, y) {
    value = matrix(Allchr_circular)
    circos.barplot(value, 0:1, col = 1:10)
  })

  


circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",
              major.at=brk,labels=round(brk/10^7,1),labels.cex=0.4,
              col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
},bg.border=F)



# this labels chrs
col_text <- "grey40"
circos.track(track.index = get.current.track.index(),
             panel.fun=function(x,y) {
               chr=CELL_META$sector.index
               xlim=CELL_META$xlim
               ylim=CELL_META$ylim
               circos.text(mean(xlim),mean(ylim),chr)
             })





circos.track(
             track.index = get.current.track.index(), 
             panel.fun = function(x, y) {
  circos.axis(h="top",get.current.sector.index())
},bg.border=F)

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",chr_names_ordered)
},bg.border=F)

circos.axis(h, sector.index, track.index)



col_text <- "grey40"


circos.track(track.index = get.current.track.index(),panel.fun=function(x,y) {
  chr=chr_names_ordered
  xlim=chr_lengths
  ylim=CELL_META$ylim
  circos.text(mean(chr_lengths),mean(ylim),chr,cex=0.5,col=col_text,
              facing="bending.inside",niceFacing=TRUE)
},bg.col="grey90",bg.border=F,track.height=0.06)





circos.track(factors = Allchr_circular$chrom, x = Allchr_circular$start, 
             y = Allchr_circular$NGA,
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = "·")
             })


circos.track(ylim = (0,1.01), panel.fun = function(x, y) {
  value = as.matrix(Allchr_circular[4:8], ncol = 5)
  circos.barplot(value, Allchr_circular$start, 
                 col = 1:5,
                 lwd = 0.1)
})



circos.initialize(fa = letters[1:4], xlim = c(0, 10))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = runif(10)
  circos.barplot(value, 1:10 - 0.5, col = 1:10)
})
circos.clear()

factors = c("a", "a", "a", "b", "b")
x = 1:5
y = 5:1
circos.track(factors = factors, x = x, y = y,
             panel.fun = function(x, y) {
               circos.points(x, y)
             })



circos.track(y = c(0, 1), panel.fun = function(x, y) {
  #value = Allchr_circular[4:8], ncol = 5)
  circos.barplot(as.matrix(Allchr_circular[4:8]), 0:1, col = 1:5)
})



# first set the color of the chrs to alternate (we need 21 colors for the 20 autosomes and the X chr)
bgcol = c(rep(c("#FF0000", "#00FF00"), 10),"#FF0000")
circos.barplot(Allchr_circular$NGA, Allchr_circular$start, bar_width = 0.6,
               col = NA, border = "black", lwd = par("lwd"), lty = par("lty"))



circos.initialize(fa = letters[1:4], xlim = c(0, 10))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = runif(10)
  circos.barplot(value, 1:10 - 0.5, col = 1:10)
})
circos.track(ylim = c(-1, 1), panel.fun = function(x, y) {
  value = runif(10, min = -1, max = 1)
  circos.barplot(value, 1:10 - 0.5, col = ifelse(value > 0, 2, 3))
})
circos.clear()

circos.initialize(factors = Allchr_circular$chrom, x = Allchr_circular$start)
circos.track(ylim = c(0, 4), panel.fun = function(x, y) {
  value = matrix(as.matrix(Allchr_circular[4:8]), ncol = 5)
  circos.barplot(value, Allchr_circular$start, col = 1:5)
})
circos.clear()
# }



circos.trackHist(Allchr_circular$chrom, Allchr_circular$NGA, 
                 bin.size = 0.02, bg.col = bgcol, col = NA)



circos.trackHist(Allchr_circular$chrom, Allchr_circular$NGA, bg.col = "blue", col = "#69b3a2")

circos.trackHist(Allchr_circular$chrom, Allchr_circular$pos, 
                 bin.size = 0.2,
                 bg.col = "white", 
                 col = "#69b3a2")

#set.seed(999)
#n = 1000
#df = data.frame(factors = sample(letters[1:8], n, replace = TRUE),
#                x = rnorm(n), y = runif(n)) # this is some random data




circos.trackHist(Allchr_circular$chrom, Allchr_circular$NGA, bg.col = "white", col = "#69b3a2")

circos.trackHist(factors = Allchr_circular$chrom, x = Allchr_circular$NGA, 
                 bin.size = 0.2,
                 col = "#FF0000",
                 border = "#00FF00")

# this is for plotting dots

circos.track(factors = Allchr_circular$chrom, y = Allchr_circular$NGA,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$factors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
# add histogram
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(df$factors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)

# add heatmap
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 0.1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = rand_color(n_breaks), border = NA)
})

circos.barplot(Allchr_circular, pos, bar_width = 0.6,
               col = NA, border = "black", lwd = par("lwd"), lty = par("lty"))

circos.track(ylim = c(0, 4), panel.fun = function(x, y) {
  circos.barplot(value = Allchr_circular, 1:10 - 0.5, col = 2:6)
})


library(BioCircos)
myGenome = list("1" = 225002135,
                "2a" = 108967917,
                "2b" = 131519175,
                "3" = 198060209,
                "4" = 190369981,
                "5" = 179725205,
                "6" = 171866349,
                "7" = 185267708,
                "8" = 144034664,
                "9" = 111027318,
                "10" = 129655328,
                "11" = 127102482,
                "12" = 133317794,
                "13" = 95368959,
                "14" = 169736342,
                "15" = 92674614,
                "16" = 74750809,
                "17" = 76917410,
                "18" = 70128972,
                "19" = 53113586,
                "X" = 148935249
                )

m <- myGenome[c("1","2a","2b","3","4","5","6","7","8","9","10","11","13","14","15","16","17","18","19","X")]

BioCircos(genome = m, yChr = FALSE, genomeFillColor = "Reds", chrPad = 0, 
          displayGenomeBorder = FALSE, genomeTicksDisplay = FALSE, genomeLabelDy = 0)
```
    

randomstuff
```
        --states AFR BNEM=nem_GumGum_female,nem_Ngsang_sumatra_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male DEN=Denisova \
        --pop-file data.yaml \
        (no need for this until we have separate files by chr):
        --rec-file rec.{CHROM}
        
--states MA=MA1+MA2+MA3+MA4+MA5 MB=MB1+MB2
        

