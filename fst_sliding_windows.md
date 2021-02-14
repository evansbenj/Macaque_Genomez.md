# Fst sliding windows

I used Martin Simon's general genomics scripts to generate Fst values in 100,000 bp nonoverlapping sliding windows. I then used a per script below to read in a concatenated file of these values for all chromosomes (generated with Makes_inputfile_for_jackknife.pl), and also a gff file with coordinates for all mRNA (all_mRNA.gff) which was made from the original gff file using grep for 'mRNA'.  This file then prints out genes in the top 99.9th percentile for Fst.

As with Twist, first I used Beagle to phase the data:
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

Then I made geno files:
```
python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i XXX.phased.vcf.gz.vcf.gz | gzip > phased_genos/XXX.geno.gz
```
then I used the popwindows script:
```
#!/bin/sh
#SBATCH --job-name=popgenWindows
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=2gb
#SBATCH --output=popgenWindows.%J.out
#SBATCH --error=popgenWindows.%J.err
#SBATCH --account=def-ben

# sbatch 2020_popgenWindows.sh  pathandname_of_geno_file name_of_pop_or_sample1 name_of_pop_or_sample2
# sbatch 2020_popgenWindows.sh ../raw_data/XTgenomez_Chr7.vcf.gz_20mil.geno.gz XT10_WZ XT11_WW

# XT10_WZ	XT11_WW	XT7_WY

#module load StdEnv/2020
#module load scipy-stack/2020b
#module load python/3.8.2
module --force purge
module load StdEnv/2020 scipy-stack/2020b python/3.8.2

echo python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3si
gmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/popgenWindows.py -g ${1} -o ${1}_${2}_${3}.csv -m 1 -
p ${2} -p ${3} -f phased -T 10 --popsFile pops.txt --writeFailedWindows -w 100000 -s 100000 -m 10 --windType coordinate

python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/
final_data_including_sites_with_lots_of_missing_data/genomics_general/popgenWindows.py -g ${1} -o ${1}_${2}_${3}.csv -m 1 -p ${2
} -p ${3} -f diplo -T 10 --popsFile pops.txt --writeFailedWindows -w 100000 -s 100000 -m 10 --windType coordinate
```
this was used to parse the results
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program reads in the output of the the general genomics popgenWindows.py script
# and find the top 99.9% of the Fst windows.  It then reads in a simplified gff file 
# that has only the mRNA and tells what genes are in these windows.

# first concatenate data from all chrs like this:
# ../Makes_inputfile_for_jackknife.pl ton_mau_windowstats

# then run like this:
# Finds_genes_near_Fst_peaks.pl ./fst_ton/ton_mau_windowstats.concat

# to get only the NON_INTERACTING OXPHOS genes add this pipe:
# Finds_genes_near_Fst_peaks.pl ./fst_ton/ton_mau_windowstats.concat | egrep 'NDUFA2|NDUFAF6|NDUFS1|NDUFV1|NUBPL|NDUFA3|NDUFA8|NDUFA13|NDUFAF5|NDUFAF7|NDUFS2|NDUFS3|NDUFS7|NDUFS8|NDUFA6|NDUFA7|NDUFA12|NDUFS4|NDUFS6|NDUFV3|SDHA|SDHAF1|SDHAF2|ACN9|C6orf57|SDHB|SDHC|SDHD|UQCRC1|UQCRC2|CYC1|UQCRH_|UQCR10|UQCR11|BCS1L|TTC19|COX4I1|COX5A|HIGD1A|=MR1_|PET100|PET117|ATPAF1|TMEM70|ATPAF2|ATP5O|ATP5C1|ATP5D|ATP5E|ATP5B|ATP5A1|ATP5G1|ATP5G2|ATP5G3|ATP5F1|ATP5H|ATP5J|ATP5I|ATP5L|ATP5J2'

# to get only the INTERACTING OXPHOS genes add this pipe:
# Finds_genes_near_Fst_peaks.pl ./fst_ton/ton_mau_windowstats.concat | egrep 'NDUFAF3|NDUFAF4|TIMMDC1|ATP5SL|NDUFB6|FOXRED1|NDUFB10|NDUFB11|NDUFB4|NDUFB5|TMEM70|ACAD9|"=COA1"|ECSIT|NDUFAF1|NDUFC1|NDUFC2|TMEM126B|TMEM186|C9orf123|NDUFAB1|NDUFB2|NDUFB3|NDUFB7|NDUFB8|NDUFB9|UQCC|MNF1|C11orf83|UQCRB|UQCRQ|TACO1|COX14|COA3|CMC1|"=COA1"|SURF1|C4orf52|COX17|COX11|COX19|COX10|COX15|COX16|COA6|SCO2|SCO1|COX20|COX18|TMEM177|COX5B|COX6C|COX7B|COX7C|COX8A|COX6A1|COX6B1|COX7A1|NDUFA4|ATPIF1|USMG5|C14orf2_'


my $inputfile = $ARGV[0];
my $cutoffpercentile = 0.999;

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my @temp;
my @temp1;
my %fst;
my @fst;

while ( my $line = <DATAINPUT>) {
	@temp=split(',',$line);
	chomp(@temp);
	if($temp[0] ne 'scaffold'){
		$fst{$temp[0]."_".$temp[1]} = $temp[8];
		if($temp[8] ne 'nan'){
			push(@fst,$temp[8]);
		}	
	}	
}		

my @fst_sorted;
@fst_sorted = sort { $a <=> $b } @fst;


print "@fst_sorted";
my $fst_cutoff_value = $fst_sorted[int($cutoffpercentile*$#fst_sorted+1)-1];
print "\n";
print "length ",$#fst_sorted+1,"\n";
print "cutoff ",$fst_cutoff_value,"\n";

# now find windows
my @high_fst_windows;
#my %high_fst_windows;
foreach my $key (keys %fst){
	if($fst{$key} > $fst_cutoff_value){
		@temp1=split('_',$key);
		push(@high_fst_windows,$key)
		#print $key,"\n";
		#$high_fst_windows{$key} = $fst{$key};
	}
}

print "high_fst_windows @high_fst_windows \n\n";

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
	# round down the start of the mRNA to nearest 10,000, which is where the fst window begins
	if (length($temp[3]) >=5 ) {
		$temp[3] =~ s/\d{5}$//; 
		$start = $temp[3]."00001";
	}	
    else{
     	$start = 1;
    }
    # print this line if it is in the high fst windown
	#if ( $temp[0]."_".$start ~~ @high_fst_windows ) {
	$value=	$temp[0]."_".$start;
	for (@high_fst_windows) {
    	if ($_ eq $value) {
       		print $line." Fst ".$fst{$value}."\n";
       		last;
    	}
	}
}


```
