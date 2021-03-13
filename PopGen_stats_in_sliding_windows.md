# PopGenome

This R script generates stat files for each population with at least 3 genome (bor, mau, ton, hec, nig):
```R
# This script will calculate lots of neutrality stats in sliding windows
# and print them to file for each of 5 populations
# bor, mau, ton, hec, and nig

## Working directory
# https://cran.r-project.org/web/packages/PopGenome/vignettes/Whole_genome_analyses_using_VCF_files.pdf
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/PopGenome")
library(ggplot2)
library(PopGenome)
library(stringr)
options(scipen = 999)

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
# chr lengths
# chr01	225002135
# chr02a	108967917
# chr02b	131519175
# chr03	198060209
# chr04	190369981
# chr05	179725205
# chr06	171866349
# chr07	185267708
# chr08	144034664
# chr09	111027318
# chr10	129655328
# chr11	127102482
# chr12	133317794
# chr13	95368959
# chr14	169736342
# chr15	92674614
# chr16	74750809
# chr17	76917410
# chr18	70128972
# chr19	53113586
# chrX	148935249

chr <- c("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09",
                 "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX")
# this is a vector of chr lenths
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
chr_lengths <- cbind(chr,ends)

# initialize some variables
bor_final = data.frame()
mau_final = data.frame()
ton_final = data.frame()
hec_final = data.frame()
nig_final = data.frame()

# read the vcf file (download from graham: 
# /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males
# also need tabix tbi file from same directory)
for(i in 1:20){ # cycle through each autosome, don't go to 21, because this is the X and it has a different name
  # and may need different settings for popgenome
  #open the vcf file for this autosome
  GENOME.class <- readVCF(paste("FandM_",chr_lengths[i,1],"_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz", sep=""),
                          numcols=10000, 
                          tid=chr_lengths[i,1], #the chr name
                          from=1, 
                          to= chr_lengths[i,2], # the chr length
                          approx=FALSE, 
                          out="", 
                          parallel=FALSE, 
                          )

  # define the populations
  bru <- c("bru_PF707")
  pap <- c("download")
  hec <- c("hecki_PF505","hecki_PF643","hecki_PF644","hecki_PF647","hecki_PF648")
  mau <- c("maura_PF615","maura_PF713","maura_PM613","maura_PM614","maura_PM616")
  sum <- c("nem_Ngsang_sumatra_female")
  bor <- c("nem_GumGum_female","nem_PM1206","nem_PM664","nem_PM665","nem_Sukai_male")
  nig <- c("nigra_PF1001","nigra_PF660","nigra_PM1003")
  nge <- c("nigrescens_PM1011","nigrescens_PM654")
  tog <- c("tog_PF549")
  ton <- c("tonk_PF511","tonk_PF559","tonk_PF563","tonk_PF597","tonk_PF626","tonk_PM592")
  
  GENOME.class <- set.populations(GENOME.class,
                                  list(bru,pap,hec,mau,sum,bor,nig,nge,tog,ton),
                                  diploid=TRUE)
  
  # define outgroup
  GENOME.class <- set.outgroup(GENOME.class,pap)
  
  # Transform object into object divided by sliding window
  win_snp <- sliding.window.transform(GENOME.class, 
                                      width=100000, jump=100000, 
                                      type=2, # if type=1, only biallelic positions, type=2 all positions
                                      whole.data=FALSE)
  
  # calculate the neutrality stats in the windows
  win_snp <- neutrality.stats(win_snp,
                   subsites=FALSE,detail=FALSE, FAST=FALSE, do.R2=FALSE)
  
  ## S4 method for signature 'GENOME'
  #get.neutrality(win_snp,theta=TRUE,stats=TRUE)
  
  
  # hec is pop3, mau is pop4,sum is pop5,bor is pop6,nig is pop7, ton is pop10
  
  ###########
  # this is the bit for bor
  bor <- get.neutrality(win_snp)[[6]]
  bor_df <- cbind(rownames(bor), data.frame(bor, row.names=NULL))
  colnames(bor_df)[1] <- "rowname"
  # split first column
  position<- str_split(gsub("-:", "", bor_df$rowname), " ", simplify=TRUE)
  bor_df <- cbind(position[,4], bor_df)
  bor_df$chr <-chr_lengths[i,1] # this is the chr we are on
  colnames(bor_df)[1] <- "position"
  new_bor_df<-bor_df[,c(12,1,3,4,6,7,9,10)]
  new_bor_df$pop <-"bor"
  bor_final <-rbind(bor_final,new_bor_df)
  # this is the end of the bit for ton
  ###########
  
  ###########
  # this is the bit for mau
  mau <- get.neutrality(win_snp)[[4]]
  mau_df <- cbind(rownames(mau), data.frame(mau, row.names=NULL))
  colnames(mau_df)[1] <- "rowname"
  # split first column
  position<- str_split(gsub("-:", "", mau_df$rowname), " ", simplify=TRUE)
  mau_df <- cbind(position[,4], mau_df)
  mau_df$chr <-chr_lengths[i,1] # this is the chr we are on
  colnames(mau_df)[1] <- "position"
  new_mau_df<-mau_df[,c(12,1,3,4,6,7,9,10)]
  new_mau_df$pop <-"mau"
  mau_final <-rbind(mau_final,new_mau_df)
  # this is the end of the bit for ton
  ###########
  
  ###########
  # this is the bit for ton
  tonk <- get.neutrality(win_snp)[[10]]
  tonk_df <- cbind(rownames(tonk), data.frame(tonk, row.names=NULL))
  colnames(tonk_df)[1] <- "rowname"
  # split first column
  position<- str_split(gsub("-:", "", tonk_df$rowname), " ", simplify=TRUE)
  tonk_df <- cbind(position[,4], tonk_df)
  tonk_df$chr <-chr_lengths[i,1] # this is the chr we are on
  colnames(tonk_df)[1] <- "position"
  new_tonk_df<-tonk_df[,c(12,1,3,4,6,7,9,10)]
  new_tonk_df$pop <-"ton"
  ton_final <-rbind(ton_final,new_tonk_df)
  # this is the end of the bit for ton
  ###########
  
  ###########
  # this is the bit for hec
  hec <- get.neutrality(win_snp)[[3]]
  hec_df <- cbind(rownames(hec), data.frame(hec, row.names=NULL))
  colnames(hec_df)[1] <- "rowname"
  # split first column
  position<- str_split(gsub("-:", "", hec_df$rowname), " ", simplify=TRUE)
  hec_df <- cbind(position[,4], hec_df)
  hec_df$chr <-chr_lengths[i,1] # this is the chr we are on
  colnames(hec_df)[1] <- "position"
  new_hec_df<-hec_df[,c(12,1,3,4,6,7,9,10)]
  new_hec_df$pop <-"hec"
  hec_final <-rbind(hec_final,new_hec_df)
  # this is the end of the bit for ton
  ###########
  
  ###########
  # this is the bit for nigra
  nig <- get.neutrality(win_snp)[[7]]
  nig_df <- cbind(rownames(nig), data.frame(nig, row.names=NULL))
  colnames(nig_df)[1] <- "rowname"
  # split first column
  position<- str_split(gsub("-:", "", nig_df$rowname), " ", simplify=TRUE)
  nig_df <- cbind(position[,4], nig_df)
  nig_df$chr <-chr_lengths[i,1] # this is the chr we are on
  colnames(nig_df)[1] <- "position"
  new_nig_df<-nig_df[,c(12,1,3,4,6,7,9,10)]
  new_nig_df$pop <-"nig"
  nig_final <-rbind(nig_final,new_nig_df)
  # this is the end of the bit for ton
  ###########

}  # end of chr loop

#View(bor_final)
#View(mau_final)
#View(ton_final)
#View(hec_final)
#View(nig_final)


write.table(bor_final, file = "stats_in_windows_bor.csv", sep = "\t",
            row.names = F, col.names = T, quote=F)
write.table(mau_final, file = "stats_in_windows_mau.csv", sep = "\t",
            row.names = F, col.names = T, quote=F)
write.table(ton_final, file = "stats_in_windows_ton.csv", sep = "\t",
            row.names = F, col.names = T, quote=F)
write.table(hec_final, file = "stats_in_windows_hec.csv", sep = "\t",
            row.names = F, col.names = T, quote=F)
write.table(nig_final, file = "stats_in_windows_nig.csv", sep = "\t",
            row.names = F, col.names = T, quote=F)
```

And this script identifies outliers:
```perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program reads in the output of the the general genomics popgenWindows.py script
# and find the top 99.9% (or whatever percentile) of the stat windows.  It then reads in a simplified gff file 
# that has only the mRNA and tells what genes are in these windows.

# must have a file called "all_mRNA.gff" in the same directory

# first use an R script to generate an root mean square file

# then run like this:
# Finds_genes_near_stat_peaks.pl stats_in_windows_nig.csv 7
# where 7 is the column of the stats
# TajD is 3
# seg sites is 4
# Fu.Li.F is 5
# Fu.Li.D is 6
# Fay.Wu.H is 7
# Zeng.E is 8

# to get only the OXPHOS genes add this pipe:
# ./Finds_genes_near_stat_window_peaks.pl stats_in_windows_nig.csv 7 | egrep 'transcript_01' | egrep 'NDUFA2|NDUFAF6|NDUFS1|NDUFV1|NDUFV2|NUBPL|NDUFA3|NDUFA8|NDUFA13|NDUFAF5|NDUFAF7|NDUFS2|NDUFS3|NDUFS7|NDUFS8|NDUFAF3|NDUFAF4|TIMMDC1|ATP5SL|NDUFB6|FOXRED1|NDUFB10|NDUFB11|NDUFB4|NDUFB5|TMEM70|ACAD9|=COA1|ECSIT|NDUFAF1|NDUFC1|NDUFC2|TMEM126B|TMEM186|C9orf123|NDUFAB1|NDUFB2|NDUFB3|NDUFB7|NDUFB8|NDUFB9|NDUFA6|NDUFA7|NDUFA12|NDUFS4|NDUFS6|NDUFV3|SDHA|SDHAF1|SDHAF2|SDHB|SDHC|SDHD|ACN9|C6orf57|UQCRB|UQCRQ|UQCRC1|UQCRC2|CYC1|UQCRH_|UQCR10|UQCR11|BCS1L|TTC19|UQCC|MNF1|C11orf83|TACO1|COX14|=COA3|CMC1|SURF1|C4orf52|COX17|COX11|COX19|COX10|COX15|COX4I1|COX5A|HIGD1A|COX17|=COX16|=COA6|=SCO2|=SCO1|COX20|COX18|TMEM177|COX5B|COX6C|COX7B|COX7C|COX8A|=MR1_|PET100|PET117|COX6A1|COX6B1|COX7A1|NDUFA4|ATPA|TMEM70|ATPIF1|USMG5|C14orf2_|ATP5O|ATP5F1|ATP5I|ATP5L|ATP5J2|ATP5C1|ATP5D|ATP5E|ATP5B|ATP5A1|ATP5G|ATP5J|ATP5H' 

# to get only the N_interact genes add this pipe:
# ./Finds_genes_near_stat_window_peaks.pl stats_in_windows_nig.csv 7 | egrep 'transcript_01'| egrep 'NDUFA2|NDUFAF6|NDUFS1|NDUFV1|NDUFV2|NUBPL|NDUFA3|NDUFA8|NDUFA13|NDUFAF5|NDUFAF7|NDUFS2|NDUFS3|NDUFS7|NDUFS8|NDUFAF3|NDUFAF4|TIMMDC1|ATP5SL|NDUFB6|FOXRED1|NDUFB10|NDUFB11|NDUFB4|NDUFB5|TMEM70|ACAD9|=COA1|ECSIT|NDUFAF1|NDUFC1|NDUFC2|TMEM126B|TMEM186|C9orf123|NDUFAB1|NDUFB2|NDUFB3|NDUFB7|NDUFB8|NDUFB9|NDUFA6|NDUFA7|NDUFA12|NDUFS4|NDUFS6|NDUFV3|SDHA|SDHAF1|SDHAF2|SDHB|SDHC|SDHD|ACN9|C6orf57|UQCRB|UQCRQ|UQCRC1|UQCRC2|CYC1|UQCRH_|UQCR10|UQCR11|BCS1L|TTC19|UQCC|MNF1|C11orf83|TACO1|COX14|=COA3|CMC1|SURF1|C4orf52|COX17|COX11|COX19|COX10|COX15|COX4I1|COX5A|HIGD1A|COX17|=COX16|=COA6|=SCO2|=SCO1|COX20|COX18|TMEM177|COX5B|COX6C|COX7B|COX7C|COX8A|=MR1_|PET100|PET117|COX6A1|COX6B1|COX7A1|NDUFA4|ATPA|TMEM70|ATPIF1|USMG5|C14orf2_|ATP5O|ATP5F1|ATP5I|ATP5L|ATP5J2|ATP5C1|ATP5D|ATP5E|ATP5B|ATP5A1|ATP5G|ATP5J|ATP5H|=MRP|C6orf203|GADD45GIP1|ARS2|POLRMT|TFAM|TFB1M|TFB2M' 

# to get only the REP genes:
# ./Finds_genes_near_stat_window_peaks.pl stats_in_windows_nig.csv 7 | egrep 'POLRMT|TFAM|TFB1M|TFB2M'

my $inputfile = $ARGV[0];
my $stat_column = $ARGV[1];
my $cutoffpercentile = 0.999;

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my @temp;
my @temp1;
my %stat;
my @stat;

while ( my $line = <DATAINPUT>) {
	@temp=split('\t',$line);
	chomp(@temp);
	if(($temp[0] ne 'chr')&&($temp[0] ne 'chrX')){ # deliberately exclude chrX
		if($temp[$stat_column-1] ne 'NA'){
			$stat{$temp[0]."_".$temp[1]} = $temp[$stat_column-1];
			push(@stat,$temp[$stat_column-1]);
		}	
	}	
}		

my @stat_sorted;
#@stat_sorted = sort { $a <=> $b } @stat;
@stat_sorted = sort { $b <=> $a } @stat; # this is reverse sorted because we care about negative ones


#print "@stat_sorted";
my $stat_cutoff_value = $stat_sorted[int($cutoffpercentile*$#stat_sorted+1)-1];
print "\n";
print "length ",$#stat_sorted+1,"\n";
print "cutoff ",$stat_cutoff_value,"\n";

# now find windows
my @high_stat_windows;
#my %high_stat_windows;
foreach my $key (keys %stat){
	if($stat{$key} < $stat_cutoff_value){ # looking for values below the cutoff
		@temp1=split('_',$key);
		push(@high_stat_windows,$key)
		#print $key,"\n";
		#$high_stat_windows{$key} = $stat{$key};
	}
}

#print "high_stat_windows @high_stat_windows \n\n";

# OK now find genes in these windows

$inputfile = "all_mRNA.gff";

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my $start;
my $value;
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	# round down the start of the mRNA to nearest 100,000, which is where the stat window begins
	if (length($temp[3]) >=5 ) {
		$temp[3] =~ s/\d{5}$//; 
		$start = $temp[3]."00001";
	}	
    else{
     	$start = 1;
    }
    # print this line if it is in the high stat windown
	#if ( $temp[0]."_".$start ~~ @high_stat_windows ) {
	$value=	$temp[0]."_".$start;
	for (@high_stat_windows) {
    	if ($_ eq $value) {
       		print $line." stat ".$stat{$value}."\n";
       		last;
    	}
	}
}

```
and this script does permutations:
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# all N-mt interact genes, their acronyms, and whether or not (1 or 0) they interact
# directly with mt genes.

# The other file is a file with stats in windows with coordinates.  First
# the mean stat of N-mt interacting (1) and non-interacting (0) genes will be calculated
# then permutations will be performed where the difference between these categories is 
# recalculated after the interaction is permuted n times. This will allow a p value of the
# stat value to be estimated.

# this program deliberately ignores all genes in chrX.  A separate script will be generated
# that is only for chrX

# execute like this:
# ./All_N_mt_allinteract_stat_column_permutation.pl FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez.txt stats_in_windows_nig.csv 7

# where the last argument is the 1-based column of the statistic being tested.
# TajD is 3
# seg sites is 4
# Fu.Li.F is 5
# Fu.Li.D is 6
# Fay.Wu.H is 7
# Zeng.E is 8

my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $stat_column = $ARGV[2];

my $window_size=100000;
my @windowsites;
my @stat_values;
my $sumsites=0;
my $counter=0;
my @temp;
my $y;
my $x;
my %OXPHOS;

# first open up the OXPHOS gene info (# OXPHOS_ARP2_MRP_andallother_genez.txt)
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(($temp[0] ne 'gene')&&($temp[2] ne 'chrX')){ # deliberately blocks with no stat value
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0];
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"complex"} = $temp[1];
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5];
	}	
}		
close DATAINPUT;


# now open up the stat data
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

my @temp1;
my $N_interact_window=0;

while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split('\t',$line);
		$N_interact_window=0;
		# first check if the window has a stat, if not ignore that window
		#if($temp[$stat_column-1] ne 'NA'){
			# cycle through each gene
			foreach my $key (keys %OXPHOS){
				@temp1=split('_',$key);
				# check if this window contains one or more N_mt genes
				if(($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= ($temp[1]+$window_size))){
						$OXPHOS{$key}{"start_stat"} = $temp[$stat_column-1];
						if($OXPHOS{$key}{"mt_interact"} == 1){
							$N_interact_window=1;
						}	
				} # start is in block
				if(($temp1[0] eq $temp[0])&&($temp1[2] >= $temp[1])&&($temp1[2] <= ($temp[1]+$window_size))) {  # end is in block
						$OXPHOS{$key}{"end_stat"} = $temp[$stat_column-1];
						if($OXPHOS{$key}{"mt_interact"} == 1){
							$N_interact_window=1;
						}	
				}
			}
		#}	
		#if($temp[0] ne 'chr'){
		#	print OUTFILE $temp[$stat_column-1],"\t",$N_interact_window,"\n";
		#}	
}

# Now the OXPHOS hash has stat and coordinates of all blocks that have genes

my @weighted_ave_stat_for_perms; # this has only the mtinteractors and all genes
my @weighted_ave_stat_for_perms_all_N_mt; # this has only the mtinteractors and the non-interactors that are still OXPHOS, MRP, or ARP2
my @weighted_ave_stat_for_perms_complex1; # this has only the mtinteracters for complex1 and the non-interactors for all complexes
my @weighted_ave_stat_for_perms_complex3; # this has only the mtinteracters for complex3 and the non-interactors for all complexes
my @weighted_ave_stat_for_perms_complex4; # this has only the mtinteracters for complex4 and the non-interactors for all complexes
my @weighted_ave_stat_for_perms_complex5; # this has only the mtinteracters for complex5 and the non-interactors for all complexes

# now calculate the weighted averages for each gene
# this is necessary because the block that has the start site may be different
# from the one that has the end site, and these may have different numbers of 
# sites.  A weighted average will upweight the one with more sites
foreach my $key (keys %OXPHOS){
	if((exists($OXPHOS{$key}{"end_stat"}))&&(exists($OXPHOS{$key}{"start_stat"}))
		&&($OXPHOS{$key}{"start_stat"} ne 'NA')
		&&($OXPHOS{$key}{"end_stat"} ne 'NA')
		){
		# build @weighted_ave_stat_for_perms that includes all genes
		$OXPHOS{$key}{"weighted_ave_stat"}=($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2;
		push(@weighted_ave_stat_for_perms,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		# now build @weighted_ave_stat_for_perms_all_N_mt 
		# this includes only N_mt genes (including direct and non-direct interactions: OXPHOS, ARP2, MRP)
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} ne 0)){ # this is a non-interacting N_mt gene
			push(@weighted_ave_stat_for_perms_all_N_mt,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} ne 0)){ 
			push(@weighted_ave_stat_for_perms_all_N_mt,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		# build @weighted_ave_stat_for_perms_complex1
		# this includes only proteins in complex 1
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '1')){ 
			push(@weighted_ave_stat_for_perms_complex1,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '1')){ 
			push(@weighted_ave_stat_for_perms_complex1,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		# build @weighted_ave_stat_for_perms_complex3
		# this includes only proteins in complex 3
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '3')){ 
			push(@weighted_ave_stat_for_perms_complex3,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '3')){ 
			push(@weighted_ave_stat_for_perms_complex3,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		# build @weighted_ave_stat_for_perms_complex4
		# this includes only proteins in complex 4
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '4')){ 
			push(@weighted_ave_stat_for_perms_complex4,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '4')){ 
			push(@weighted_ave_stat_for_perms_complex4,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		# build @weighted_ave_stat_for_perms_complex5
		# this includes only proteins in complex 5
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '5')){ 
			push(@weighted_ave_stat_for_perms_complex5,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '5')){ 
			push(@weighted_ave_stat_for_perms_complex5,($OXPHOS{$key}{"start_stat"}+$OXPHOS{$key}{"end_stat"})/2);
		}
	}
	#else{
	#	print "Problem with stat value for ",$key,"\n";
	#	print "The problem is that there is no stat value for the gene\n";
	#}							
}	

close DATAINPUT2;

my $stat_associated=0; # all OXPHOS,MRP, ARP2 genes
my $stat_associated_complex1=0;
my $stat_associated_complex3=0;
my $stat_associated_complex4=0;
my $stat_associated_complex5=0;

my $stat_non_associated=0;
my $stat_non_associated_only_N_mt=0;
my $stat_non_associated_complex1=0;
my $stat_non_associated_complex3=0;
my $stat_non_associated_complex4=0;
my $stat_non_associated_complex5=0;


my $n_stat_associated=0; # this also works for only_N_mt
my $n_stat_associated_complex1=0;
my $n_stat_associated_complex3=0;
my $n_stat_associated_complex4=0;
my $n_stat_associated_complex5=0;

my $n_stat_non_associated=0; # this includes all non associated genes, including those anywhere 
my $n_stat_non_associated_only_N_mt=0; # this includes only non associated N-mt genes 
									  # (i.e., OXPHOS, MRP, ARP2 that don't directly interact with mt genes)
my $n_stat_non_associated_complex1=0;
my $n_stat_non_associated_complex3=0;
my $n_stat_non_associated_complex4=0;
my $n_stat_non_associated_complex5=0;

my $N_interact_counter=0;

# now calculate the average stat for associated and non-associated OXPHOS genes
foreach my $key (keys %OXPHOS){
	if(exists($OXPHOS{$key}{"weighted_ave_stat"})){
		if($OXPHOS{$key}{"mt_interact"} == 1){
			$N_interact_counter+=1;
			$stat_associated += $OXPHOS{$key}{"weighted_ave_stat"};
			$n_stat_associated += 1;	
			if($OXPHOS{$key}{"complex"} eq 1){
				$stat_associated_complex1+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_associated_complex1 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 3){
				$stat_associated_complex3+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_associated_complex3 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 4){
				$stat_associated_complex4+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_associated_complex4 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 5){
				$stat_associated_complex5+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_associated_complex5 +=1;
			}	
		}
		elsif($OXPHOS{$key}{"mt_interact"} == 0){
			$stat_non_associated += $OXPHOS{$key}{"weighted_ave_stat"};
			$n_stat_non_associated += 1;
			if($OXPHOS{$key}{"complex"} ne 0){
				$stat_non_associated_only_N_mt += $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_non_associated_only_N_mt += 1;
			}
			if($OXPHOS{$key}{"complex"} eq 1){
				$stat_non_associated_complex1+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_non_associated_complex1 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 3){
				$stat_non_associated_complex3+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_non_associated_complex3 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 4){
				$stat_non_associated_complex4+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_non_associated_complex4 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 5){
				$stat_non_associated_complex5+= $OXPHOS{$key}{"weighted_ave_stat"};
				$n_stat_non_associated_complex5 +=1;
			}	
		}
	}
}


# now report values
print "Number of associated blocks ",$n_stat_associated,"\n";
print "Number of N_interact genes ",$N_interact_counter,"\n";
print "Mean stat of all associated N_mt genes ",$stat_associated/$n_stat_associated,"\n";
print "Mean stat complex 1 associated ",$stat_associated_complex1/$n_stat_associated_complex1,"\n";
print "Mean stat complex 3 associated ",$stat_associated_complex3/$n_stat_associated_complex3,"\n";
print "Mean stat complex 4 associated ",$stat_associated_complex4/$n_stat_associated_complex4,"\n";
print "Mean stat complex 5 associated ",$stat_associated_complex5/$n_stat_associated_complex5,"\n";
print "Mean stat non-associated all genes ",$stat_non_associated/$n_stat_non_associated,"\n";
#print "Mean stat non-associated only N_mt ",$stat_non_associated_only_N_mt/$n_stat_non_associated_only_N_mt,"\n";


#################
# ALL COMPLEXES perms including all genes
#################

# calculate test statistic for all genes
my $test_stat = ($stat_associated/$n_stat_associated) - ($stat_non_associated/$n_stat_non_associated);

# do permutations for all_complex_mt_interact vs all_other_genes
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
my @associated_or_not_array = (('1') x $n_stat_associated, ('0') x $n_stat_non_associated);

my $perms=1000;
my $stat_associated_perm=0;
my $stat_not_associated_perm=0;
my $n_stat_associated_perm=0;
my $n_stat_not_associated_perm=0;

# first check if the length of the permutation array is the same as the stat array
# for this analysis
if($#associated_or_not_array ne $#weighted_ave_stat_for_perms){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length stat array ",$#weighted_ave_stat_for_perms,"\n";
}

my @perm_diffs;
@perm_diffs={};
for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$stat_associated_perm=0;
	$stat_not_associated_perm=0;
	$n_stat_associated_perm=0;
	$n_stat_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$stat_associated_perm+= $weighted_ave_stat_for_perms[$x];
			$n_stat_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$stat_not_associated_perm+= $weighted_ave_stat_for_perms[$x];
			$n_stat_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($stat_associated_perm/$n_stat_associated_perm) - ($stat_not_associated_perm/$n_stat_not_associated_perm));	
}

my @perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
my $switch=0;
my $pval=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for test including all genes (if negative, then not significant):",$test_stat,"\n";
print "P = ",1-($pval/$perms),"\n";




# fisher_yates_shuffle( \@array ) : 
    # generate a random permutation of @array in place
    sub fisher_yates_shuffle {
        my $array = shift;
        my $i;
        for ($i = @$array; --$i; ) {
            my $j = int rand ($i+1);
            next if $i == $j;
            @$array[$i,$j] = @$array[$j,$i];
        }
    }

```

