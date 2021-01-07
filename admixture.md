Great directions are here: https://speciationgenomics.github.io/ADMIXTURE/


working directory: 
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data
```

```
mkdir ADMIXTURE
cd ADMIXTURE
```
use bcftools to combine phased SNPs into one file to feed into plink
```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.10.2
bcftools concat -o autosomes.vcf ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz 
```
compress and index
```
bgzip -c autosomes.vcf > autosomes.vcf.gz
tabix -p vcf autosomes.vcf.gz
```
convert to geno format using plink
first make a bed file and remove any SNP with no data for autosomes:

```
module load nixpkgs/16.09  intel/2016.4 plink/1.9b_5.2-x86_64

plink --vcf ./autosomes.vcf.gz --make-bed --geno 0.999 --out ./autosomes --allow-extra-chr --const-fid
```
and then separately for chrX
```
plink --vcf ../all_diploid_haploid_chrX_phased.vcf.gz.vcf.gz --make-bed --geno 0.999 --out ./chrX --allow-extra-chr --const-fid
```
for autosomes only, we need to change the chr names in the .bim file because these cause problems for admixture:
```
awk -v OFS='\t' '{$1=0;print $0}' autosomes.bim > autosomes.bim.tmp
mv autosomes.bim.tmp autosomes.bim
```


now run admixture for k=2
```
module load StdEnv/2020 nixpkgs/16.09 admixture/1.3.0
admixture --cv chrX.bed 2 > chrXlog2.out
```
or try K=2 to 9
```
for i in {2..9}
do
 admixture --cv autsomes.bed $i > autosomeslog${i}.out
done
```
For plotting, first make a file with individual and species names
```
awk '{split($2,name,"_"); print name[2],name[1]}' chr01.nosex > chr01.list
```
edit to fix papio



download R script to plot admixture 
```
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r
chmod +x plotADMIXTURE.r
```
plot it
```
module load StdEnv/2020 r/4.0.2
Rscript plotADMIXTURE.r -p chr01 -i chr01.list -k 2 -l papio,nem,nigra,nigrescens,hecki,tonk,tog,maura,bru 
```
I made a custom plot that I like better:
```
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/Admixture/chrX')
library (ggplot2)
library(reshape2) # used to make tall data for plotting stacked barplot
library(egg) # used for making grid plots

species <- factor(species, levels= c("nemestrina","nigra",
                                         "nigrescens","hecki","tonkeana",
                                         "togeanus","maura","brunnescens"), ordered = T)
color <- factor(color, levels= c("blue","purple",
                                     "pink","maroon","black",
                                     "gray","orange","lime"), ordered = T)


dat<-read.table("./chrX.2.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
  )
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
dat$Species <- c(
  "Papio",
  "nemstrina",
  "nemstrina","nemstrina","nemstrina","nemstrina","nemstrina",
  "nigra","nigra","nigra",
  "nigrescens","nigrescens",
  "hecki","hecki","hecki","hecki","hecki",
  "tonkeana","tonkeana","tonkeana","tonkeana","tonkeana","tonkeana",
  "togeanus",
  "maura","maura","maura","maura","maura",
  "brunnescens")
tall_dat <- melt(dat, id.vars=c("Sample","Species"))

k2<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "Sample") +
  labs(y="2") +
  scale_x_discrete(position = "top") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust=-2, hjust=0))+ 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  theme(axis.title.x = element_text(size=18)) +
  # custom colors
  scale_fill_manual(values=c("blue","purple")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none") 


dat<-read.table("./chrX.3.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k3<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="3") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("red","purple","blue")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")


dat<-read.table("./chrX.4.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k4<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="4") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("purple","blue","red","orange")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")


dat<-read.table("./chrX.5.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4","anc5")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k5<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="5") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("blue","orange","red","purple","yellow")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")


dat<-read.table("./chrX.6.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4","anc5",
                   "anc6")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k6<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="6") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("red","steel blue","yellow","orange","purple","blue")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")


dat<-read.table("./chrX.7.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4","anc5",
                   "anc6","anc7")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k7<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="7") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("yellow","lightgoldenrod4","orange","purple","blue","steel blue","red")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")


dat<-read.table("./chrX.8.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4","anc5",
                   "anc6","anc7","anc8")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k8<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="8") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("yellow","lightgoldenrod4","purple","red","blue","plum2","steel blue","orange")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")

dat<-read.table("./chrX.9.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4","anc5",
                   "anc6","anc7","anc8","anc9")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k9<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="9") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("yellow","red","lightgoldenrod4","blue","steel blue","orange",
                             "sienna","plum2","purple")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")


dat<-read.table("./chrX.10.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4","anc5",
                   "anc6","anc7","anc8","anc9","anc10")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)
tall_dat <- melt(dat, id.vars="Sample")

k10<- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(y="10") +
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # custom colors
  scale_fill_manual(values=c("plum2","yellow","blue","red","purple","hotpink",
                             "steel blue","orange","sienna","lightgoldenrod4")) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) +
  # get rid of the legend
  theme(legend.position = "none")


dat<-read.table("./chrX.11.Q",header=FALSE)
colnames(dat) <- c("anc1","anc2","anc3","anc4",
                   "anc5","anc6","anc7","anc8",
                   "anc9","anc10","anc11")
dat$Sample <- c(
  "PF707",
  "Papio",
  "PF505","PF643","PF644","PF647","PF648",
  "PF615","PF713","PM613","PM614","PM616",
  "GumGum","Ngsang","PM1206","PM664","PM665","Sukai",
  "PF1001","PF660","PM1003",
  "PM1011","PM654",
  "PF549",
  "PF511","PF559","PF563","PF597","PF626","PM592"
)
dat$Sample <- factor(dat$Sample, levels=c(
  "Papio",
  "Ngsang",
  "Sukai","GumGum","PM1206","PM664","PM665",
  "PF660","PF1001","PM1003",
  "PM1011","PM654",
  "PF647","PF648","PF505","PF643","PF644",
  "PF511","PF563","PF559","PM592","PF597","PF626",
  "PF549",
  "PF615","PM616","PM614","PF713","PM613",
  "PF707"), ordered = T)

tall_dat <- melt(dat, id.vars="Sample")


k11 <- ggplot(tall_dat, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  labs(x = "Species") +
  labs(y="11") +
  theme(axis.text.x = element_text(size=12)) +  
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
  theme(axis.title.x = element_text(size=18)) +
  #theme(axis.text.y=element_text(size=12)) +
  scale_fill_manual(values=c("red","lightgoldenrod4","yellow","steel blue","plum2","sienna",
                             "purple","lightseagreen","orange","blue","hotpink")) +
  # get rid of the legend
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("Papio" = expression(alpha),
                            "GumGum" = expression(gamma),"Ngsang" = expression(beta),"PM1206" = expression(gamma),"PM664" = expression(gamma),
                            "PM665" = expression(gamma), "Sukai" = expression(gamma),
                            "PF1001" = expression(epsilon),"PF660" = expression(epsilon),"PM1003" = expression(epsilon),
                            "PM1011" = expression(zeta),"PM654" = expression(zeta),
                            "PF505" = expression(eta), "PF643" = expression(eta),"PF644" = expression(eta),"PF647" = expression(eta),"PF648" = expression(eta),
                            "PF511" = expression(theta),"PF559" = expression(theta),"PF563" = expression(theta),"PF597" = expression(theta), "PF626" = expression(theta), "PM592" = expression(theta),
                            "PF549" = expression(iota),
                            "PF615" = expression(kappa),"PF713" = expression(kappa),"PM613" = expression(kappa),"PM614" = expression(kappa),"PM616" = expression(kappa),
                            "PF707" = expression(lambda))) +
  # custom y axis tick labels
  scale_y_continuous(breaks=seq(0,1,1)) #+ 
 # theme(plot.margin = unit(c(0.05,0.05,2,0.25), "cm"))


pdf("./2020_chrX_admixture.pdf",w=8, h=6, version="1.4", bg="transparent")
ggarrange(k2, k3, k4, k5, k6,
          k7, k8, k9, k10, k11, 
          nrow = 10,
          ncol = 1)
dev.off()

```

(From the admixture manual), we can also check the cross validation standard error like this:
```
grep -h CV chrXlog*.out

CV error (K=2): 0.68009
CV error (K=3): 0.56954
CV error (K=4): 0.55752
CV error (K=5): 0.50140
CV error (K=6): 0.48437
CV error (K=7): 0.55287
CV error (K=8): 0.35196
CV error (K=9): 0.69166
CV error (K=10): 0.57981
CV error (K=11): 0.73068
```
This suggests the best one is k=8 for chrX; this is substantially better than the others.

for autosomes it is k=9:
```
grep -h CV autosomeslog*.out
CV error (K=2): 0.63166
CV error (K=3): 0.63444
CV error (K=4): 0.52661
CV error (K=5): 0.49589
CV error (K=6): 0.49198
CV error (K=7): 0.59790
CV error (K=8): 0.58931
CV error (K=9): 0.43763
CV error (K=10): 0.68583
CV error (K=11): 0.63831
```
