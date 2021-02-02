Sliding Windows Analyses

Simon Martin has nice software that calculates D and fdm stats in sliding windows.  I also wrote a program to do this but it makes more sense to use his since it is probably better vetted than mine (!).

I am working in this directory on graham:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing
```
Phase vcf files like this:
```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=def-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dept
h_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/beagle.18May20.d20.jar gt=${1} out=${1}_phased.vcf.gz impu
te=true
```

Convert my filtered vcf files to geno format like this:
```
python parseVCF.py -i ../../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz -o chr1.geno.gz
```

Then it is necessary to swap any astrisks with Ns:
```
gunzip chr18.geno.gz
sed -i 's/\*/N/g' chr18.geno 
gzip -c chr18.geno > chr18.geno.gz
```

for autosomes:
```
#!/bin/sh
#SBATCH --job-name=abba
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=8gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

# sbatch ABBABABA.sh chr H1 H2 H3 O
# sbatch ABBABABA.sh chr01 nig nge hec papio

# populations
# bru papio hec mau nem sum nig nge tog ton


module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2

echo python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5
}.csv -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFai
ledWindows --windType coordinate

python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5}.csv
 -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFailedWi
ndows --windType coordinate
```

for chrX:
```
#!/bin/sh
#SBATCH --job-name=abba
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

# sbatch ABBABABA.sh chr H1 H2 H3 O
# sbatch ABBABABA_chrX.sh chrX nga nge hec papio

# populations
# bru papio hec mau nem sum nig nge tog ton


module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2

echo python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5
}.csv -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFai
ledWindows --windType coordinate --haploid maura_PM613,maura_PM614,maura_PM616,nem_PM1206,nem_PM664,nem_PM665,nem_Suka
i_male,nigra_PM1003,nigrescens_PM1011,nigrescens_PM654,tonk_PM592

python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5}.csv
 -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFailedWi
ndows --windType coordinate --haploid maura_PM613,maura_PM614,maura_PM616,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_mal
e,nigra_PM1003,nigrescens_PM1011,nigrescens_PM654,tonk_PM592
```

This is working well now. I can concatenate autosomes with the script below called "Makes_inputfile_for_jackknife.pl" and do the block jackknife with this file using the script below "Does_block_jackknife.pl".

```perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program globs a bunch of chromosome output files from the script called Performs_ABBA_BABA_on_populations.pl
# and concatenates them for analysis by "Does_block_jackknife.pl"

# Run like this "Makes_inputfile_for_jackknife.pl inputprefix"

# The output prefix will be: inputprefix.concat

my $inputfile = $ARGV[0];


my @files = glob("'*${inputfile}*'");

print "hi @files\n";

# open an outputfile
unless (open(OUTFILE, ">$inputfile.concat"))  {
	print "I can\'t write to $inputfile.concat\n";
	exit;
}
print "Creating output file: $inputfile.concat\n";


unless (open DATAINPUT, $files[0]) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	print OUTFILE $line;
}		
close DATAINPUT;

my $counter=0;

foreach my $infile (1..$#files){
	unless (open DATAINPUT, $files[$infile]) {
		print "Can not find the input file.\n";
		exit;
	}
	$counter=0;
	while ( my $line = <DATAINPUT>) {
		if($counter>0){
			print OUTFILE $line;
		}
		$counter+=1;
	}		
	close DATAINPUT;
}	

close OUTFILE
```

```perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program reads in the output of the script called Performs_ABBA_BABA_on_populations.pl
# and calculates the standard error of the weighted mean value of fDM with weightings based 
# on the sum of the number of abba and baba sites in each window.

my $inputfile = $ARGV[0];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my @sumsites;
my $fDM_numsites;
my @fDM_numsites;
my @D_values;
my @weighted_D_values;
my @fdm_values;
my @weighted_fdm_values;
my $counter=0;
my $y;
my $x;
my @temp; 
my $weighted_fdm_values;
my $ABBA=0;
my $BABA=0;
my $BBAA=0;

while ( my $line = <DATAINPUT>) {
	@temp=split(',',$line);
	if($temp[0] ne 'scaffold'){
		# if fDM is not numeric, it shouldn't contribute to the weighted average
		if(($temp[6] !~ /nan/)&&($temp[7] !~ /nan/)&&($temp[8] !~ /nan/)&&($temp[10] !~ /nan/)){
			# load the number of ABBA and BABA sites for each window
			$sumsites[$counter]=$temp[6]+$temp[7];
			# load the number of fDM sites for each window
			$fDM_numsites[$counter]=$temp[6]+$temp[7];
			# add them up to eventually get the average
			$fDM_numsites+=$temp[6]+$temp[7];
			# also load the fDM stat for each window
			$D_values[$counter]=$temp[8];
			# keep track of how many windows are loaded
			$fdm_values[$counter]=$temp[10];
			# keep track of how many windows are loaded
			$counter+=1;
			$ABBA+=$temp[6];
			$BABA+=$temp[7];
		}
	}	
}		

# Calculate real stat
# make $fDM_numsites the average per site
$fDM_numsites=$fDM_numsites/($counter);
my $weighted_D_average=0;
my $weighted_fdm_average=0;
# calculate the weighted average of fDM
for ($y = 0 ; $y <= $#fdm_values; $y++ ) {
	$weighted_D_average+=($fDM_numsites[$y]*$D_values[$y]/$fDM_numsites);
	$weighted_fdm_average+=($fDM_numsites[$y]*$fdm_values[$y]/$fDM_numsites);
}	
print "The number of ABBA sites is ",$ABBA,"\n";
print "The number of BABA sites is ",$BABA,"\n";
print "The average number of ABBABABA sites per window is ",$fDM_numsites,"\n";
#print "The non-weighted average is ",sprintf("%.6f",$non_weighted_fdm_value/($counter)),"\n";
$weighted_D_average=$weighted_D_average/($counter);
$weighted_fdm_average=$weighted_fdm_average/($counter);
print "The weighted D average is ",sprintf("%.6f",$weighted_D_average),"\n";
print "The weighted fDM average is ",sprintf("%.6f",$weighted_fdm_average),"\n";

# now calculate the standard error.
my @jack_sumsites;
my $jack_averagesites;
my @jack_fdm_values;
my $jack_weighted_average=0;
my @jack_weighted_fdm_values;
my @jackarray;
my $counter2=0;

for ($y = 0 ; $y < $counter; $y++ ) {
	for ($x = 0 ; $x < $counter; $x++ ) {
		# leave out one row for each jackknfe replicates
		if($y != $x){
			# load the number of fDM sites for each window
			$jack_sumsites[$counter2]=$fDM_numsites[$x];
			# add them up to eventually get the average
			$jack_averagesites+=$fDM_numsites[$x];
			# also load the fDM stat for each window
			$jack_fdm_values[$counter2]=$fdm_values[$x];
			# keep track of the number of windows
			$counter2+=1;			
		}
	}
	# make the $jack_averagesites an average
	$jack_averagesites=$jack_averagesites/($counter2);
	# calculate the weighted average
	for ($x = 0 ; $x < $counter2; $x++ ) {
		$jack_weighted_average+=($jack_sumsites[$x]*$jack_fdm_values[$x]/$jack_averagesites);
	}
	push(@jackarray,($jack_weighted_average/$counter2));
	# reset variables
	$jack_weighted_average=0;
	$jack_averagesites=0;
	@jack_fdm_values=();
	$counter2=0;
}


# now calculate the standard error.
my @jack_Dsumsites;
my $jack_Daveragesites;
my @jack_D_values;
my $jack_Dweighted_average=0;
my @jack_Dweighted_fdm_values;
my @jackDarray;
$counter2=0;

for ($y = 0 ; $y < $counter; $y++ ) {
	for ($x = 0 ; $x < $counter; $x++ ) {
		# leave out one row for each jackknfe replicates
		if($y != $x){
			# load the number of fDM sites for each window
			$jack_Dsumsites[$counter2]=$fDM_numsites[$x];
			# add them up to eventually get the average
			$jack_Daveragesites+=$fDM_numsites[$x];
			# also load the fDM stat for each window
			$jack_D_values[$counter2]=$D_values[$x];
			# keep track of the number of windows
			$counter2+=1;			
		}
	}
	# make the $jack_averagesites an average
	$jack_Daveragesites=$jack_Daveragesites/($counter2);
	# calculate the weighted average
	for ($x = 0 ; $x < $counter2; $x++ ) {
		$jack_Dweighted_average+=($jack_Dsumsites[$x]*$jack_D_values[$x]/$jack_Daveragesites);
	}
	push(@jackDarray,($jack_Dweighted_average/$counter2));
	# reset variables
	$jack_Dweighted_average=0;
	$jack_Daveragesites=0;
	@jack_D_values=();
	$counter2=0;
}

# now calculate the variance of the jackknife replicates

# first we need the mean
my $jack_Dmean=0;
for ($x = 0 ; $x < $counter; $x++ ) {
	$jack_Dmean+=$jackDarray[$x];
}
$jack_Dmean=$jack_Dmean/($counter);

my $jack_Dvar=0;
for ($x = 0 ; $x < $counter; $x++ ) {
	$jack_Dvar+=($jack_Dmean-$jackDarray[$x])**2;
}

# for the sample variance, divide by (n-1)
print "jackDvar ",sprintf("%.9f",$jack_Dvar/($counter-1)),"\n";
my $sterr = sqrt($counter*($jack_Dvar/($counter-1)));
print "The standard error of the weighted D is ",sprintf("%.5f",$sterr),"\n";
print "The 95\%CI of the weighted D is ",
sprintf("%.6f",($weighted_D_average-1.96*$sterr))," - ",sprintf("%.6f",($weighted_D_average+1.96*$sterr)),"\n";

close DATAINPUT;

```

# 2020 plotting D and fDM across the Makassar Strait:
```
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/abbababa_stats/fDM_means_and_CIs')
library (ggplot2)
library(gridExtra)

nemsuladat<-read.table("./2020_Makassar_D_fDM.txt",header=TRUE)
# make a color column
nemsuladat$H3 <- factor(nemsuladat$H3, 
                 levels = c("all","nga","nge","hec","ton","tog","mau","bru"),
                 ordered = is.ordered(nemsuladat$H3))
nemsuladat$color <- "black"
nemsuladat$color[which(nemsuladat$XorA == 'X')]  <- 'red'

# make the grouping variable a factor so that the dodge works
nemsuladat$XorA <- factor(nemsuladat$XorA, 
                                 levels = c("Autosomes","X"))


# modify the facet labels
nemsuladat$H1 <- factor(nemsuladat$H1, levels = c("Sum"), 
                  labels = c("H1: Sumatra; H2: Borneo"))

pd <- position_dodge(0.5) # move them .05 to the left and right
nemsulaDstat_plot<- ggplot(nemsuladat, aes(x=H3, y=D, group=H3, colour=XorA)) +
  facet_grid(~ H1, scales = "free_x", space = "free_x") + 
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  geom_point(position = pd, size=4, alpha = 0.8, aes(group = XorA)) +
  geom_errorbar(aes(ymin=Dlower, ymax=Dupper, group = XorA), alpha = 0.8, width=.1, position=pd) +
  labs(x = expression(paste(italic(H3))), y=expression(paste(italic(D)))) +
  theme(axis.text=element_text(size=5), axis.title=element_text(size=18)) +
  geom_hline(yintercept=0) +
  scale_color_manual(breaks = c("Autosomes", "X"),values=c("black", "red")) +
  # get rid of legend title
  theme(legend.title = element_blank()) +
  theme(legend.position="none")+
  # put labels on outside
  theme(strip.placement = "outside")+
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(aspect.ratio=8)
nemsulaDstat_plot


nemsulafDMstat_plot<- ggplot(nemsuladat, aes(x=H3, y=fDM, group=H3, colour=XorA)) +
  facet_grid(~ H1, scales = "free_x", space = "free_x") + 
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  geom_point(position = pd, size=4, alpha = 0.8, aes(group = XorA)) +
  geom_errorbar(aes(ymin=fDMlower, ymax=fDMupper, group = XorA), alpha = 0.8, width=.1, position=pd) +
  labs(x = expression(paste(italic(H3))), y=expression(paste(italic(f[DM])))) +
  theme(axis.text=element_text(size=5), axis.title=element_text(size=18)) +
  geom_hline(yintercept=0) +
  scale_color_manual(breaks = c("Autosomes", "X"),values=c("black", "red")) +
  # get rid of legend title
  theme(legend.title = element_blank()) +
  theme(legend.position="none")+
  # put labels on outside
  theme(strip.placement = "outside") +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(aspect.ratio=8)


sulanemdat<-read.table("./2020_MakassarSula_D_fDM.txt",header=TRUE)
# make a color column
sulanemdat$H3 <- factor(sulanemdat$H3, 
                        levels = c("all","nga","nge","hec","ton","tog","mau","bru"),
                        ordered = is.ordered(sulanemdat$H3))
sulanemdat$color <- "black"
sulanemdat$color[which(sulanemdat$XorA == 'X')]  <- 'red'

# modify the order of the H1 taxa
sulanemdat$H1 <- factor(sulanemdat$H1, 
                        levels = c("nga","nge","tog","bru"),
                        ordered = is.ordered(nemsuladat$H1))


# modify the facet labels
sulanemdat$H2 <- factor(sulanemdat$H2, levels = c("hec","ton","mau"), 
                        ordered = is.ordered(sulanemdat$H2),
                        labels = c("H2: hec",
                                   "H2: ton",
                                   "H2: mau"))


pd <- position_dodge(0.5) # move them .05 to the left and right
sulanemDstat_plot<- ggplot(sulanemdat, aes(x=H1, y=D, group=H2, colour=XorA)) +
  facet_grid(~ H2, scales = "free_x", space = "free_x") + 
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  geom_point(position = pd, size=4, alpha = 0.8, aes(group = XorA)) +
  geom_errorbar(aes(ymin=Dlower, ymax=Dupper, group = XorA), 
                alpha = 0.8, 
                width=.1, 
                position=pd) +
  labs(x = expression(paste(italic(H1))), y=element_blank()) +
  theme(axis.text=element_text(size=5), axis.title=element_text(size=18)) +
  geom_hline(yintercept=0) +
  scale_color_manual(breaks = c("Autosomes", "X"),values=c("black", "red")) +
  # get rid of legend title
  theme(legend.title = element_blank()) +
  theme(legend.position="none")+
  # put labels on outside
  theme(strip.placement = "outside")+
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18)) +
  theme(axis.text.y=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

sulanemfDMstat_plot<- ggplot(sulanemdat, aes(x=H1, y=fDM, group=H3, colour=XorA)) +
  facet_grid(~ H2, scales = "free_x", space = "free_x") + 
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  geom_point(position=pd, size=4, alpha = 0.8, aes(group = XorA)) +
  geom_errorbar(aes(ymin=fDMlower, ymax=fDMupper, group = XorA), alpha = 0.8, width=.1, position=pd) +
  labs(x = expression(paste(italic(H1))), y=element_blank()) +
  theme(axis.text=element_text(size=5), axis.title=element_text(size=18)) +
  geom_hline(yintercept=0) +
  scale_color_manual(breaks = c("Autosomes", "X"),values=c("black", "red")) +
  # get rid of legend title
  theme(legend.title = element_blank()) +
  theme(legend.position="none")+
  # put labels on outside
  theme(strip.placement = "outside")+
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18)) +
  theme(axis.text.y=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) 

# make the y axes the same
nemsulaDstat_plot <- nemsulaDstat_plot + scale_y_continuous(limits=c(-0.06,0.16))
sulanemDstat_plot <- sulanemDstat_plot + scale_y_continuous(limits=c(-0.06,0.16))
nemsulafDMstat_plot <- nemsulafDMstat_plot + scale_y_continuous(limits=c(-0.03,0.03))
sulanemfDMstat_plot <- sulanemfDMstat_plot + scale_y_continuous(limits=c(-0.03,0.03))

pdf("./D_means_and_CIs_Makassar.pdf",w=10, h=4, version="1.4", bg="transparent")
grid.arrange(nemsulaDstat_plot, sulanemDstat_plot, ncol=2)
dev.off()

pdf("./fDM_means_and_CIs_Makassar.pdf",w=10, h=4, version="1.4", bg="transparent")
grid.arrange(nemsulafDMstat_plot, sulanemfDMstat_plot, ncol=2)
dev.off()

```

