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

permutations:
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# all N-mt interact genes, their acronyms, and whether or not (1 or 0) they interact
# directly with mt genes.

# The other file is a file with Fst (or pi) in windows with coordinates.  First
# the mean Fst of N-mt interacting (1) and non-interacting (0) genes will be calculated
# then permutations will be performed where the difference between these categories is 
# recalculated after the interaction is permuted n times. This will allow a p value of the
# Fst value to be estimated.

# this program deliverately ignores all genes in chrX.  A separate script will be generated
# that is only for chrX

# execute like this:
# ./All_N_mt_allinteract_fst_permutation.pl FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez_orientation.txt fst_bor/nem_mau_windowstats.concat

my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];

my @windowsites;
my @Fst_values;
my $sumsites=0;
my $counter=0;
my @temp;
my $y;
my $x;
my %OXPHOS;

# first open up the OXPHOS gene info (# FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez_orientation.txt)
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(($temp[0] ne 'gene')&&($temp[2] ne 'chrX')){ # deliberately ignores chrX
		if($temp[6] eq '+'){ # the gene is in the forward orientation
			$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0];
			$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"complex"} = $temp[1];
			$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5];
		}
		elsif($temp[6] eq '-'){ # the gene is in the reverse orientation
			$OXPHOS{$temp[2]."_".$temp[4]."_".$temp[3]}{"gene"} = $temp[0];
			$OXPHOS{$temp[2]."_".$temp[4]."_".$temp[3]}{"complex"} = $temp[1];
			$OXPHOS{$temp[2]."_".$temp[4]."_".$temp[3]}{"mt_interact"} = $temp[5];
		}
		else{
			print "something wrong with gene orientation $line\n";
		}
	}	
}		
close DATAINPUT;

# this will print out a file that I may use for plotting later

my $outputfile = $inputfile2."_FST__density.txt"; # the name of the output file is from the commandline
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}

# now open up the Fst data
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

my @temp1;
my $N_interact_window=0;
my $gene_containing_window=0;
my $number_of_genes_in_this_window=0;
my $number_of_Ninteract_genes_in_this_window=0;
my $Ninteract_acronym="";
my $number_of_Ninteract_genes_spanning_a_window=0;

print OUTFILE "chr\tpos\tFst\tcontainsgenes\tcontainsNinteractgenez\tnumber_of_genes\tnumber_of_Ninteractgenez\tNinteract_acronym\n";

while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split(',',$line);
		$N_interact_window=0;
		$gene_containing_window=0;
		$number_of_genes_in_this_window=0;
		$number_of_Ninteract_genes_in_this_window=0;
		$Ninteract_acronym="-";
		# cycle through each gene
		foreach my $key (keys %OXPHOS){
			@temp1=split('_',$key);
			# check if this window contains one or more N_mt genes
			if(($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= $temp[2])){
					$gene_containing_window=1; # only count the start of genes in windows
					$number_of_genes_in_this_window+=1; # only count the start of genes in windows
					$OXPHOS{$key}{"start_fst"} = $temp[8];
					#$OXPHOS{$key}{"start_fst_sites"} = $temp[4];
					if($OXPHOS{$key}{"mt_interact"} == 1){
						$N_interact_window=1;
						$number_of_Ninteract_genes_in_this_window+=1; # only count the start of genes in windows
						if($number_of_Ninteract_genes_in_this_window == 1){
							$Ninteract_acronym=$OXPHOS{$key}{"gene"};
						}	
						else{
							$Ninteract_acronym=$Ninteract_acronym.",".$OXPHOS{$key}{"gene"};
						}	
					}	
			} # start is in block
		}
		if($temp[0] ne 'scaffold'){
			print OUTFILE $temp[0],"\t",$temp[1],"\t",$temp[8],"\t",$gene_containing_window,"\t",$N_interact_window,
			"\t",$number_of_genes_in_this_window,"\t",$number_of_Ninteract_genes_in_this_window,"\t",$Ninteract_acronym,"\n";
		}	
}

close OUTFILE;
close DATAINPUT2;

# Now the OXPHOS hash has Fst and coordinates of all blocks that have genes
my @fst_for_perms; # this has only the mtinteractors and all genes
my $Fst_associated=0; # all OXPHOS,MRP, ARP2 genes
my $Fst_non_associated=0;
my $n_Fst_associated=0; # this also works for only_N_mt
my $n_Fst_non_associated=0; # this includes all non associated genes, including those anywhere 
my $N_interact_counter=0;

# now calculate the average fst for associated and non-associated OXPHOS genes
foreach my $key (keys %OXPHOS){
	if((exists($OXPHOS{$key}{"start_fst"})) &&
		($OXPHOS{$key}{"start_fst"} ne 'nan')){

		if($OXPHOS{$key}{"mt_interact"} == 1){
			$N_interact_counter+=1;
			$Fst_associated += $OXPHOS{$key}{"start_fst"};
			$n_Fst_associated += 1;	
		}
		elsif($OXPHOS{$key}{"mt_interact"} == 0){
			$Fst_non_associated += $OXPHOS{$key}{"start_fst"};
			$n_Fst_non_associated += 1;
		}
		push(@fst_for_perms,$OXPHOS{$key}{"start_fst"});
	}
}


# now report values
print "Number of associated blocks ",$n_Fst_associated,"\n";
print "Number of N_interact genes ",$N_interact_counter,"\n";
print "Mean Fst all associated N_mt genes ",$Fst_associated/$n_Fst_associated,"\n";
print "Mean Fst non-associated all genes ",$Fst_non_associated/$n_Fst_non_associated,"\n";
#print "Mean Fst non-associated only N_mt ",$Fst_non_associated_only_N_mt/$n_Fst_non_associated_only_N_mt,"\n";


#################
# ALL COMPLEXES perms including all genes
#################

# calculate test statistic for all genes
my $test_stat = ($Fst_associated/$n_Fst_associated) - ($Fst_non_associated/$n_Fst_non_associated);

# do permutations for all_complex_mt_interact vs all_other_genes
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
my @associated_or_not_array = (('1') x $n_Fst_associated, ('0') x $n_Fst_non_associated);

my $perms=1000;
my $Fst_associated_perm=0;
my $Fst_not_associated_perm=0;
my $n_Fst_associated_perm=0;
my $n_Fst_not_associated_perm=0;

# first check if the length of the permutation array is the same as the fst array
# for this analysis
if($#associated_or_not_array ne $#fst_for_perms){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#fst_for_perms,"\n";
}

my @perm_diffs;
@perm_diffs={};
for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$Fst_associated_perm=0;
	$Fst_not_associated_perm=0;
	$n_Fst_associated_perm=0;
	$n_Fst_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$Fst_associated_perm+= $fst_for_perms[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $fst_for_perms[$x];
			$n_Fst_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($Fst_associated_perm/$n_Fst_associated_perm) - ($Fst_not_associated_perm/$n_Fst_not_associated_perm));	
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
