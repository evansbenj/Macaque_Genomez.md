# Testing for LD between nonsynonymous SNPs in the mtDNA and variation in the nuclear genome

Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data
```

First I added a Papio anubis mtDNA genome (GenBank accessionm NC 020006) to the alignment and separated the mtDNA genomes by the 13 protein coding genes.  I exported the protein coding sequences, opened them each in paup and excluded gapped, constant and uninf (parsimony non-informative) positions.  The last category are autapomorphic sites.

Then I further filtered these positions as follows in order to remove variation that is solely attributable to geographic structure:
removed autapomorphic variation
removed variation fixed within species
removed gapped and constant positions
removed variation that was only present within only one species
removed variation that was present only in papio and one macaque individual, population, or species
removed variation that was present in only tonk + tog, even if not fixed in tonk (IBD)
kept variation that was present in some heck plus tonk 511
removed variation where papio was divergent and remaining polymorphism was in only one macaque population or species
removed variation fixed on Sulawesi, fixed in borneo, or fixed in nem, including if this variation was shared with papio.

In other words, only variation that was present in more than one macaque species was considered


Plink
Associations between SNPs and a phenotype (such as phenotypic sex) are easily calculated using plink.
First make plink files out of the vcf files:
```
module load nixpkgs/16.09
module load plink/1.9b_5.2-x86_64
plink --vcf FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz --recode --const-fid 0 --set-missing-var-ids @:# --allow-extra-chr --out chr19_plink
```
where this flag `--set-missing-var-ids @:# ` tells plink to set the SNP ids as "chr:bp".  This was important because the permutation test reports significant SNPs based on SNPID, which was blank without this flag.  The '--allow-extra-chr' flag was needed for chr02a and chr02b.

or, for chrX:
```
plink --vcf all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode.vcf.gz --recode --const-fid 0 --set-missing-var-ids @:# --allow-extra-chr --out chrX_plink
```

Now test for associations for each SNP for each chromosome
```
plink --file chr01_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr01_ND1_
plink --file chr02a_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr02a_ND1_
plink --file chr02b_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr02b_ND1_
plink --file chr03_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr03_ND1_
plink --file chr04_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr04_ND1_
plink --file chr05_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr05_ND1_
plink --file chr06_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr06_ND1_
plink --file chr07_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr07_ND1_
plink --file chr08_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr08_ND1_
plink --file chr09_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr09_ND1_
plink --file chr10_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr10_ND1_
plink --file chr11_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr11_ND1_
plink --file chr12_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr12_ND1_
plink --file chr13_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr13_ND1_
plink --file chr14_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr14_ND1_
plink --file chr15_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr15_ND1_
plink --file chr16_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr16_ND1_
plink --file chr17_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr17_ND1_
plink --file chr18_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr18_ND1_
plink --file chr19_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr19_ND1_
plink --file chrX_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chrX_ND1_

```
where the "ND1_base.txt" file is a tab-delimited file that looks like this:
```
0	bru_PF707	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	download	2	1	2	1	1	0	0	1	1	1	1	2	2	2	1	1	1	2
0	hecki_PF505	1	1	2	1	1	1	1	1	1	1	1	1	0	1	0	1	2	1
0	hecki_PF643	1	1	1	1	1	1	2	2	1	1	2	1	2	1	2	1	1	1
0	hecki_PF644	1	1	1	1	1	1	2	2	1	1	2	1	2	1	2	1	1	1
0	hecki_PF647	1	1	2	1	1	1	1	1	1	1	1	1	0	1	0	1	2	1
0	hecki_PF648	1	1	2	1	1	1	1	1	1	1	1	1	0	1	0	1	2	1
0	maura_PF615	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	maura_PF713	2	1	1	1	1	1	1	1	1	2	1	2	1	2	1	1	1	1
0	maura_PM613	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	maura_PM614	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	maura_PM616	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	nem_GumGum_female	1	2	1	2	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nem_Ngsang_sumatra_female	1	1	1	0	2	1	1	1	1	2	1	2	2	2	1	2	2	1
0	nem_PM1206	1	1	2	1	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nem_PM664	1	1	1	1	1	2	1	1	2	1	1	1	1	1	1	1	1	2
0	nem_PM665	1	1	2	1	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nem_Sukai_male	1	2	1	2	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nigra_PF1001	1	1	1	1	1	1	1	1	1	1	1	1	2	0	1	1	1	1
0	nigra_PF660	1	1	1	1	1	1	1	1	1	1	1	1	2	0	1	1	2	1
0	nigra_PM1003	1	1	1	1	1	1	1	1	1	1	1	1	2	0	1	1	1	1
0	nigrescens_PM1011	1	1	1	1	1	1	1	1	1	1	1	1	1	0	1	1	1	2
0	nigrescens_PM654	1	1	1	1	1	1	1	1	1	1	1	1	1	0	1	1	1	2
0	tog_PF549	1	1	1	2	2	2	1	1	1	1	1	1	1	1	1	2	1	1
0	tonk_PF511	1	1	1	1	1	1	2	2	1	1	2	1	2	1	2	1	1	1
0	tonk_PF559	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	tonk_PF563	1	1	1	1	2	1	1	1	1	1	1	1	1	1	2	1	1	1
0	tonk_PF597	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	tonk_PF626	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	tonk_PM592	2	2	1	1	1	1	1	1	2	1	1	1	1	1	1	1	1	1
```
The first column is the family ID (just zeros here). The second column is the sample name - this is the same as in the vcf file. The third column is the phenotype - should use 1 and 2, NOT 1 and 0, because 0 might be interpreted as a missing phenotye.  In positions where I have three amino acids, I encoded the 3rd one with a zero if it was population specific.  In some cases I had to split the variation up when three amino acids were present but none seemed to be based on population structure (e.g. ATP8 amino acid 9).

Also this flag "--const-fid 0" sets the family id to zero and tells plink to use the vcf sample name as the sample ID irrespective of whether there is an underscore in the name.

This flag "--chr-set 36" allows extra chrs. They will be numbers in the order they are encountered in the vcf file (I think - this will need to be confirmed...)

This --all-pheno tells plink to iterate over all phenotypes

This flag --pfilter 1e-4 tells plink to only output SNPs with Pvals <= 0.0001 (otherwise the files are huge!)

Results can be presented for each mtDNA polymorphism by gene like this (for ND1):
```
# for computecandada, load these modules
# module load StdEnv/2020  gcc/9.3.0 r-bundle-bioconductor/3.12
# module load nixpkgs/16.09  gcc/7.3.0 r/4.0.2

# I have the input files for ND1 here on graham: /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/
with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/ND1
# but I couldn't get the karyoploteR to load on graham so I moved them to cedar scratch instead

# for computecandada, set this dir on cedar
# setwd('/home/ben/scratch/ATP8_aa_phenotypes')
setwd('/home/ben/scratch/ND1_aa_phenotypes/ND1')
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("karyoploteR")




library(regioneR)
library(GenomicRanges)
library(karyoploteR)
library(plyr)
# loop through all of the SNPs and make a graph for each one
for (i in 1:18) {  # for ND1 there are 18 interesting SNPs

# for computecandada, set this dir on cedar
#mydir = "/home/ben/scratch/ATP8_aa_phenotypes"
mydir = "/home/ben/scratch/ND1_aa_phenotypes/ND1"
# for computecandada, change the file name
myfiles1 = list.files(path=mydir, pattern=paste("chr.+_ND1_.P",i,".assoc$",sep=""), full.names=TRUE)
#myfiles2 = list.files(path=mydir, pattern="popgenWindows_chr.+_nem_nig_windowstats", full.names=TRUE)
#myfiles3 = list.files(path=mydir, pattern="popgenWindows_chr.+_nem_ton_windowstats", full.names=TRUE)
#myfiles4 = list.files(path=mydir, pattern="popgenWindows_chr.+_nem_mau_windowstats", full.names=TRUE)


# have to change order so that chr02a and chr02b are read properly
myfiles1 <- myfiles1[c(2,3,1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
# now it is ok
dat1 = ldply(myfiles1, read.table, sep = "", fill=TRUE, header = TRUE, row.names = NULL)
#dat2 = ldply(myfiles2, read.table, sep = ",", fill=TRUE, header = TRUE)
#dat3 = ldply(myfiles3, read.table, sep = ",", fill=TRUE, header = TRUE)
#dat4 = ldply(myfiles4, read.table, sep = ",", fill=TRUE, header = TRUE)

dat1$index <- paste(dat1$CHR,"_",dat1$BP, sep="")
#dat2$index <- paste(dat2$scaffold,"_",dat2$start, sep="")
#dat3$index <- paste(dat3$scaffold,"_",dat3$start, sep="")
#dat4$index <- paste(dat4$scaffold,"_",dat4$start, sep="")


#temp_dat1 <- merge(dat1, dat2, by = c("index","index"))
#temp_dat2 <- merge(dat3, dat4, by = c("index","index"))
#all_dat <- merge(temp_dat1, temp_dat2, by = c("index","index"))


# get rid of rows with NAs
#new_dat <- all_dat[complete.cases(all_dat[ ,]),]
new_dat <- dat1[complete.cases(dat1[ ,]),]
# now calculate the root mean square for each row
# important to use normalized values, otherwise the RMS may be dominated
# by one comparison

#######
#new_dat$RMS <- sqrt((0.25)*((new_dat$Fst_nem_hec**2)+(new_dat$Fst_nem_nig**2)+
#                              (new_dat$Fst_nem_ton**2)+(new_dat$Fst_nem_mau**2)))
#new_dat$mean_xpclr_norm <- (0.25)*(new_dat$Fst_nem_hec+new_dat$Fst_nem_nig+
#                                     new_dat$Fst_nem_ton+new_dat$Fst_nem_mau)

# rename some columns
#######
colnames(new_dat)[1] <- "chrom"
colnames(new_dat)[9] <- "pval"
#colnames(new_dat)[3] <- "start"
#colnames(new_dat)[4] <- "stop"
# make the start begin in the middle of the window
#new_dat$start <- new_dat$start + 50000

# I need to rename chr in a very specific way
new_dat$chrom <- gsub('\\b23\\b', 'chrX', new_dat$chrom)
new_dat$chrom <- gsub('\\b19\\b', 'chr19', new_dat$chrom)
new_dat$chrom <- gsub('\\b18\\b', 'chr18', new_dat$chrom)
new_dat$chrom <- gsub('\\b17\\b', 'chr17', new_dat$chrom)
new_dat$chrom <- gsub('\\b16\\b', 'chr16', new_dat$chrom)
new_dat$chrom <- gsub('\\b15\\b', 'chr15', new_dat$chrom)
new_dat$chrom <- gsub('\\b14\\b', 'chr14', new_dat$chrom)
new_dat$chrom <- gsub('\\b13\\b', 'chr13', new_dat$chrom)
new_dat$chrom <- gsub('\\b12\\b', 'chr12', new_dat$chrom)
new_dat$chrom <- gsub('\\b11\\b', 'chr11', new_dat$chrom)
new_dat$chrom <- gsub('\\b10\\b', 'chr10', new_dat$chrom)
new_dat$chrom <- gsub('\\b9\\b', 'chr09', new_dat$chrom)
new_dat$chrom <- gsub('\\b8\\b', 'chr08', new_dat$chrom)
new_dat$chrom <- gsub('\\b7\\b', 'chr07', new_dat$chrom)
new_dat$chrom <- gsub('\\b6\\b', 'chr06', new_dat$chrom)
new_dat$chrom <- gsub('\\b5\\b', 'chr05', new_dat$chrom)
new_dat$chrom <- gsub('\\b4\\b', 'chr04', new_dat$chrom)
new_dat$chrom <- gsub('\\b3\\b', 'chr03', new_dat$chrom)
new_dat$chrom <- gsub('\\b1\\b', 'chr01', new_dat$chrom)

custom.genome <- toGRanges(data.frame(chrom=c("chr01","chr02a","chr02b",
                                            "chr03","chr04","chr05",
                                            "chr06","chr07","chr08",
                                            "chr09","chr10","chr11",
                                            "chr12","chr13","chr14",
                                            "chr15","chr16","chr17",
                                            "chr18","chr19","chrX"), 
                                      start=rep(1,21), 
                                      end=c(225002135,108967917,131519175,
                                            198060209,190369981,179725205,
                                            171866349,185267708,144034664,
                                            111027318,129655328,127102482,
                                            133317794,95368959,169736342,
                                            92674614,74750809,76917410,
                                            70128972,53113586,148935249)))


GR_data <- makeGRangesFromDataFrame(new_dat,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field="chrom",
                                    start.field="BP",
                                    end.field="BP",
                                    #strand.field=NULL,
                                    starts.in.df.are.0based=FALSE)



plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$dataideogramheight <- 0
plot.params$data1height <- 6
plot.params$data2height <- 0
plot.params$data1inmargin <- 4
plot.params$data2inmargin <- 0
plot.params$data1outmargin <- 80
plot.params$data2outmargin <- 20
plot.params$bottommargin <- 10
plot.params$leftmargin <- 0.12

# where is the lowest P value?
#new_dat[(-log10(new_dat$P) == max(-log10(new_dat$P))),]

png(filename = paste("ND1_plink_LD_",i,".png",sep=""),w=800, h=800,
    units = "px", bg="white")
kp <- plotKaryotype(genome = custom.genome, plot.type=1, plot.params = plot.params,
                    ideogram.plotter = NULL, labels.plotter = NULL)
kpAddChromosomeNames(kp, xoffset=70, yoffset= 0.01)
kpDataBackground(kp, data.panel = 1, r0 = 0 , r1 = 20)
kpAxis(kp, r0 = 0 , r1 = 20, ymin=0, ymax = 20, numticks = 3, cex=0.5, data.panel=1)  
kpPoints(kp, chr = new_dat$chrom, 
         x=new_dat$BP, y=-log10(new_dat$pval),
         col=colByChr(GR_data))

# COMPLEX 1 (final)
# 3 genes interact directly with ND1:
kpPlotMarkers(kp, chr="chr03", x=104846075, y=c(20), labels="", label.dist = 0.01, line.color="red", lwd=2) # NDUFAF3
kpPlotMarkers(kp, chr="chr06", x=95211059, y=c(20), labels="", label.dist = 0.01, line.color="red", lwd=2) # NDUFAF4
kpPlotMarkers(kp, chr="chr03", x=157411293, y=c(20), labels="", label.dist = 0.01, line.color="red", lwd=2) # TIMMDC1

# many more may not
kpPlotMarkers(kp, chr="chr05", x=138178553, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr08", x=93359100, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr02b", x=93979633, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr11", x=6690185, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr14", x=92789615, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr19", x=9634493, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr09", x=94025131, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr19", x=19317399, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr15", x=46796983, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr02a", x=37732212, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr01", x=135484416, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr11", x=18260513, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr19", x=1153942, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr11", x=6532924, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr19", x=36878373, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr09", x=59806087, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr11", x=118312703, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr16", x=1930096, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chrX", x=47194137, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr03", x=156309713, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr03", x=84594254, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr08", x=72064967, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr03", x=148084830, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr07", x=64839692, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr19", x=11316389, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr14", x=17717208, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr04", x=139226852, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr11", x=69474342, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr11", x=77465139, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr16", x=8788484, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr09", x=34632092, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr16", x=22236971, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr07", x=166649885, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr02b", x=88692897, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr19", x=14340049, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr10", x=95956681, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr08", x=123481402, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)

kpPlotMarkers(kp, chr="chr15", x=84086730, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr19", x=8302280, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr12", x=94150466, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr05", x=54356765, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr05", x=1601204, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
kpPlotMarkers(kp, chr="chr07", x=3770437, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)

# COMPLEX 5 (final)
# 3 genes interact directly with mt-ATP6 and mt-ATP8:
#kpPlotMarkers(kp, chr="chr01", x=26938483, y=c(20), labels="", label.dist = 0.01, line.color="red", lwd=2)  # ATPIF1 (IF1)
#kpPlotMarkers(kp, chr="chr10", x=99032999, y=c(20), labels="", label.dist = 0.01, line.color="red", lwd=2) # USMG5 (ATP5MD
, DAPIT) # this is really chr 10
#kpPlotMarkers(kp, chr="chr14", x=166432920, y=c(20), labels="", label.dist = 0.01, line.color="red", lwd=2) # C14orf2 (ATP
5MPL, MP68, PLPM, MLQ,6,8PL) # this is really chr 15

# many more may not
#kpPlotMarkers(kp, chr="chr01", x=45953218, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr08", x=72064967, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr17", x=18045894, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr07", x=12722004, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr10", x=7723482, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr19", x=1002804, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr15", x=5385192, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr12", x=54910245, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr18", x=35035298, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr17", x=44211923, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr12", x=51828101, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr02b", x=62111164, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr01", x=112030533, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr17", x=68443001, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr07", x=21199179, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr04", x=634919, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr11", x=110459798, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)
#kpPlotMarkers(kp, chr="chr07", x=124582469, y=c(20), labels="", label.dist = 0.01, line.color="blue", lwd=2)

dev.off()

} # end of loop
```

# Association with multiple phenotypes (mtDNA SNPs)

Make bed/bim/fam files from map/bed files

```
plink --file chr19_plink --make-bed --out chr19_plink
```
