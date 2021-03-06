# Abbababa tests  

I have some files on goblin here:
```
/work/ben/2017_SEAsian_macaques/
```

* Genotype chrX by depth

from here `/work/ben/2017_SEAsian_macaques/ben_scripts`

```
sqsub -r 2d --mpp 16G -o chrX_geno_dept.log bash -c "perl 10_Genotypes_only_male_chrX_based_on_allelic_depth.pl /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered.vcf.gz 00000000000000000000000 /work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM_chrX_BSQR_jointgeno_allsites_filtered_females_and_males_haploid.vcf.gz.tab"
```

* Run abbababa on autosomes

```
./Wrapper_for_Performs_ABBA_BABA_populations_H3_hecki_H1_maura_H2_tonk.pl
```
script "Performs_ABBA_BABA_on_populations_withH3pi.pl":

```
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and performs the population ABBA-BABA test using one outgroup
#  sequence and one or more individual sequences from three other
#  taxa (H3, H2, and H1). 

#  This analysis will include all positions that have data from at least 
#  one individual from species H3, H2, and H1, including those that that 
#  missing data in some individuals.

# to execute type Performs_ABBA_BABA_on_populations.pl inputfile.tab 1111100110000111100011100110010100000000 
# 3_6_14-18-19-20_2-3-4-5-6-7_32-33-34-35-36-37-38-39-40 output1.txt output2.txt
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_6 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case for rhesus) and (ii) the column number of the first individual in the ingroup 
# (6 in this case), and 14-18-19-20, 2-3-4-5-6-7, and 32-33-34-35-36-37-38-39-40 refer to H3, H1, and H2 samples as 
# itemized below (here, they are Borneo nemestrina, hecki, and tonkeana, respectively). 

# output1.txt is the one that is fed into the jacknife.R script
# output2.txt has more information for each window, including BBAA sites, fdom and dxy

# IMPORTANT: (i) and (ii) are columns beginning with 1 but (iii) is based on the individual samples such
# as enumerated below

# Notes: 
# ##### nigra_PM1000 is actually nigrescens_PM1000 #####
# ##### these samples have very low coverage (<10x): ####
# ##### nigrescens_PF654_sorted (7.33X) ####
# ##### maura_PM613_sorted (8.65X) ####
# ##### ochreata_PM596_sorted (9.02X)####
# ##### nigra_660_sorted (9.61X) ####
# ##### togeanus_PF549 (9.63X) ####

# Here is the order of the samples for the SulaRad project:

#	1	brunescens_PF707_stampy_sorted
#	2	hecki_PF643_stampy_sorted
#	3	hecki_PF644_stampy_sorted
#	4	hecki_PF648_stampy_sorted
#	5	hecki_PF651_stampy_sorted
#	6	hecki_PM639_stampy_sorted
#	7	hecki_PM645_stampy_sorted
#	8	maura_PF615_stampy_sorted
#	9	maura_PF713_stampy_sorted
#	10	maura_PM613_stampy_sorted
#	11	maura_PM614_stampy_sorted
#	12	maura_PM616_stampy_sorted
#	13	maura_PM618_stampy_sorted
#	14	nem_Gumgum_stampy_sorted
#	15	nem_Kedurang_stampy_sorted
#	16	nem_Malay_stampy_sorted
#	17	nem_Ngasang_stampy_sorted
#	18	nem_PM664_stampy_sorted
#	19	nem_PM665_stampy_sorted
#	20	nem_Sukai_male_stampy_sorted
#	21	nem_pagensis_stampy_sorted
#	22	nigra_PF1001_stampy_sorted
#	23	nigra_PF660_stampy_sorted
#	24	nigrescens_PM1000_stampy_sorted
#	25	nigra_PM1003_stampy_sorted
#	26	nigrescens_PF654_stampy_sorted
#	27	ochreata_PF625_stampy_sorted
#	28	ochreata_PM571_stampy_sorted
#	29	ochreata_PM596_stampy_sorted
#	30	togeanus_PF549_stampy_sorted
#	31	togeanus_PM545_stampy_sorted
#	32	tonk_PF515_stampy_sorted
#	33	tonk_PM561_stampy_sorted
#	34	tonk_PM565_stampy_sorted
#	35	tonk_PM566_stampy_sorted
#	36	tonk_PM567_stampy_sorted
#	37	tonk_PM582_stampy_sorted
#	38	tonk_PM584_stampy_sorted
#	39	tonk_PM592_stampy_sorted
#	40	tonk_PM602_stampy_sorted

# born
# 14-18-19-20
# sum
# 15-16-17-21
# tonk
# 32-33-34-35-36-37-38-39-40
# heck
# 2-3-4-5-6-7_
# maura
# 8-9-10-11-12-13
# nigrescens
# 24-26
# nigra
# 22-23-25
# togeanus
# 30-31
# och_brun
# 27-28-29-1


# for example, with a tab file with the rhesus, baboon and human sequence , here is the input command using the rhesus as the outgroup:

# tonk_heck
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_2-3-4-5-6-7_32-33-34-35-36-37-38-39-40 heck_tonk_born_pop.abbababa
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_32-33-34-35-36-37-38-39-40_2-3-4-5-6-7 tonk_heck_born_pop.abbababa

# nigra_och/brun
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_22-23-25_27-28-29-1 nigra_ochbrun_borneo_pop.abbababa
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_15-16-17-21_22-23-25_27-28-29-1 nigra_ochbrun_sum_pop.abbababa

# heck_maura
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_2-3-4-5-6-7_8-9-10-11-12-13 hecki_maura_borneo_pop.abbababa
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_15-16-17-21_2-3-4-5-6-7_8-9-10-11-12-13 hecki_maura_sum_pop.abbababa



my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile = $ARGV[3];
my $outputfile2 = $ARGV[4];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}


my @temp;
my $range;
my $line_number=0;
my $counter=0;



unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2\n";
	exit;
}
print "Creating output file: $outputfile2\n";


my @sexes = split("",$ARGV[1]);
my @whotoinclude = split("_",$ARGV[2]);
my @H3=split("-",$whotoinclude[2]);
my @H1=split("-",$whotoinclude[3]);
my @H2=split("-",$whotoinclude[4]);

if($#whotoinclude != 4){
	print "Problem with number of taxa in commandline\n";
}



my @temp1;
my $previous= 0;
my $string;
my $m;
my $n;
my $w;
my $y;
my $x;
my @unique;
my $x_uniq;

my $number_of_H3_genotyped=($#H3 + 1);
my $number_of_H2_genotyped=($#H2 + 1);
my $number_of_H1_genotyped=($#H1 + 1);

my $number_of_individuals_genotyped=$#H3+$#H2+$#H1+3;

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped=0;
for ($y = 0 ; $y <= $#H1 ; $y++ ) {
	if($sexes[$H1[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
for ($y = 0 ; $y <= $#H2 ; $y++ ) {
	if($sexes[$H2[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	
for ($y = 0 ; $y <= $#H3 ; $y++ ) {
	if($sexes[$H3[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";

my $sliding_window=5000000;
my $current_window=0;
my $window_counter=0;
my $current_chromosome="blah";


my %ABBA_hash;
my $ABBA_peak_hash;
my %BABA_hash;
my $BABA_peak_hash;
my %BBAA_hash;
my $A;
my $B;
my @allelez;
my $derived;
my $ancestral;
my $H3_derived_freq;
my $H1_derived_freq;
my $H2_derived_freq;
my $H3_ancestral_freq;
my $H1_ancestral_freq;
my $H2_ancestral_freq;
my $peak_H2_H3_derived_freq;
my $peak_H1_H3_derived_freq;
my @H3allelez;
my @H2allelez;
my @H1allelez;
my $diff;
my $num_comparisons;
my %H2_H3_pairwise_divergence_per_window;
my $diffH1H3;
my $diffH1H2;
my $num_comparisonsH1H3;
my $num_comparisonsH1H2;
my %H1_H3_pairwise_divergence_per_window;
my %H1_H2_pairwise_divergence_per_window;
my %number_of_sites_per_window;
my %number_of_sites_per_windowH1H3;
my %number_of_sites_per_windowH1H2;
my $ABBA_peak_hashH1H3;
my $BABA_peak_hashH1H3;
my %H1_pairwise_nucleotide_diversity_per_window;
my %H2_pairwise_nucleotide_diversity_per_window;
my %H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window;
my %H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window;
my %H3_pairwise_nucleotide_diversity_per_window;
my %H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window;


my $diffH2;
my $num_comparisonsH2;
my @H2_first_allelez;
my @H2_second_allelez;
my $aa;
my $bb;

my $diffH1;
my $num_comparisonsH1;
my @H1_first_allelez;
my @H1_second_allelez;

my $diffH3;
my $num_comparisonsH3;
my @H3_first_allelez;
my @H3_second_allelez;


my %f_H2_H3;
my %f_H1_H3;
my %f_dm;
my %f_dm_counter;

my %f_H2_H3_counter;
my %f_H1_H3_counter;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		print "The outgroup sequence is ",$temp[$whotoinclude[0]-1],"\n";
		print "H3 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H3; $y++ ) {
				print $temp[$H3[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
		print "H1 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H1; $y++ ) {
				print $temp[$H1[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
		print "H2 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H2; $y++ ) {
				print $temp[$H2[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
	}
	else{	
		if($temp[0] ne $current_chromosome){
			$current_chromosome = $temp[0];
			$current_window = 0;
			$window_counter+=1;
			$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$number_of_sites_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$number_of_sites_per_windowH1H3{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$number_of_sites_per_windowH1H2{$window_counter."_".$current_chromosome."_".$current_window}=0;
		}
		until($temp[1] < ($current_window+$sliding_window)){
			$current_window = $current_window+$sliding_window;
			$window_counter+=1;
			$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$number_of_sites_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$number_of_sites_per_windowH1H3{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$number_of_sites_per_windowH1H2{$window_counter."_".$current_chromosome."_".$current_window}=0;
		}
		if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")&&($temp[0] ne "chrM")){
			$string=();
			if((uc $temp[$whotoinclude[0]-1] eq "A")||
			   (uc $temp[$whotoinclude[0]-1] eq "C")||
			   (uc $temp[$whotoinclude[0]-1] eq "T")||
			   (uc $temp[$whotoinclude[0]-1] eq "G")){
				# the outgroup is defined
				$string=();
				$A = uc $temp[$whotoinclude[0]-1];
				$string=$string.$A;
				# now calculate the frequency of the derived allele in the H3 data
				$H3_derived_freq=0;
				$H3_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H3; $y++ ) {
					if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)){
						@allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
					}	
				}
				if(($derived+$ancestral)>0){
					$H3_derived_freq=$derived/($derived+$ancestral);
					$H3_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H1 data
				$H1_derived_freq=0;
				$H1_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H1; $y++ ) {
					if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3)){
						@allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
					}
				}
				if(($derived+$ancestral)>0){
					$H1_derived_freq=$derived/($derived+$ancestral);
					$H1_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H2 data
				$H2_derived_freq=0;
				$H2_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H2; $y++ ) {
					if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)){
						@allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
					}	
				}
				if(($derived+$ancestral)>0){
					$H2_derived_freq=$derived/($derived+$ancestral);
					$H2_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# only consider sites for which there are data for at least one genotype for H1, H2, and H3
				if((defined($string))&&(($H1_derived_freq>0)||($H1_ancestral_freq>0))&&(($H2_derived_freq>0)||($H2_ancestral_freq>0))&&(($H3_derived_freq>0)||($H3_ancestral_freq>0))){
					# Now calculate the average pairwise divergence between H2 and H3
					$diff=0;
					$num_comparisons=0;
					@H3allelez=();
					@H2allelez=();
					for ($y = 0 ; $y <= $#H3; $y++ ) {
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'T/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'G/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'C/*')
							&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)){
							for ($w = 0 ; $w <= $#H2; $w++ ) {
								if(($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'A/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'T/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'G/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'C/*')
									&&(length($temp[$H2[$w]+$whotoinclude[1]-2]) == 3)){
									# there are data for H2 and H3
									@H3allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
									@H2allelez=split('/',$temp[$H2[$w]+$whotoinclude[1]-2]);
									foreach(@H3allelez){
										if($_ ne $H2allelez[0]){
											$diff+=1; 
										}
										if($_ ne $H2allelez[1]){
											$diff+=1;
										}
										$num_comparisons+=2;	
									}
								}
							}
						}					
					}
					# Now calculate the average pairwise divergence between H1 and H3
					$diffH1H3=0;
					$num_comparisonsH1H3=0;
					@H3allelez=();
					@H1allelez=();
					for ($y = 0 ; $y <= $#H3; $y++ ) {
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'T/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'G/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'C/*')
							&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)){
							for ($w = 0 ; $w <= $#H1; $w++ ) {
								if(($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'A/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'T/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'G/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'C/*')
									&&(length($temp[$H1[$w]+$whotoinclude[1]-2]) == 3)){
									# there are data for H1 and H3
									@H3allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
									@H1allelez=split('/',$temp[$H1[$w]+$whotoinclude[1]-2]);
									foreach(@H3allelez){
										if($_ ne $H1allelez[0]){
											$diffH1H3+=1; 
										}
										if($_ ne $H1allelez[1]){
											$diffH1H3+=1;
										}
										$num_comparisonsH1H3+=2;	
									}
								}
							}
						}					
					}
					# Now calculate the average pairwise divergence between H1 and H2
					$diffH1H2=0;
					$num_comparisonsH1H2=0;
					@H2allelez=();
					@H1allelez=();
					for ($y = 0 ; $y <= $#H2; $y++ ) {
						if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'A/*')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'T/*')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'G/*')
							&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'C/*')
							&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)){
							for ($w = 0 ; $w <= $#H1; $w++ ) {
								if(($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'A/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'T/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'G/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'C/*')
									&&(length($temp[$H1[$w]+$whotoinclude[1]-2]) == 3)){
									# there are data for H1 and H2
									@H2allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
									@H1allelez=split('/',$temp[$H1[$w]+$whotoinclude[1]-2]);
									foreach(@H2allelez){
										if($_ ne $H1allelez[0]){
											$diffH1H2+=1; 
										}
										if($_ ne $H1allelez[1]){
											$diffH1H2+=1;
										}
										$num_comparisonsH1H2+=2;	
									}
								}
							}
						}					
					}
					# now calculate the average pairwise divergence for H2 and H3
					if($num_comparisons>0){
						# this does not depend on the position being polymorphic
							$number_of_sites_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
							$H2_H3_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diff/$num_comparisons);
							# later we will standardize this by the number of sites per window
					}
					# now calculate the average pairwise divergence for H1 and H3
					if($num_comparisonsH1H3>0){
						# this does not depend on the position being polymorphic
							$number_of_sites_per_windowH1H3{$window_counter."_".$current_chromosome."_".$current_window}+=1;
							$H1_H3_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH1H3/$num_comparisonsH1H3);
							# later we will standardize this by the number of sites per window
					}
					# now calculate the average pairwise divergence for H1 and H2
					if($num_comparisonsH1H2>0){
						# this does not depend on the position being polymorphic
							$number_of_sites_per_windowH1H2{$window_counter."_".$current_chromosome."_".$current_window}+=1;
							$H1_H2_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH1H2/$num_comparisonsH1H2);
							# later we will standardize this by the number of sites per window
					}
					# Now calculate the average pairwise nucleotide diversity within H2
					# this approach is justified here: http://binhe.org/2011/12/29/calculate-nucleotide-diversity-per-base-pair-summation-method/
					# but with a formula that is probably much quicker then the stuff below pi = (2j(n-j) / n(n-1) ), where j is the number of minor
					# alleles out of n alleles total.
					$diffH2=0;
					$num_comparisonsH2=0;
					@H2_first_allelez=();
					@H2_second_allelez=();
					for ($y = 0 ; $y < $#H2; $y++ ) {
						for ($w = ($y+1) ; $w <= $#H2; $w++ ) {
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'G/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'C/*')
								&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3) 
								&& ($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'C/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'G/*')
								&&(length($temp[$H2[$w]+$whotoinclude[1]-2]) == 3)){
								# there are data for both H2 alleles
								@H2_first_allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
								@H2_second_allelez=split('/',$temp[$H2[$w]+$whotoinclude[1]-2]);
								# combine the alleles into one array
								@H2_first_allelez = (@H2_first_allelez, @H2_second_allelez);
								# check that this array has 4 elements
								if($#H2_first_allelez != 3){
									print "Problem with number of alleles @H2_first_allelez @H2_second_allelez\n";
								}
								else{ # calculate the average pairwise diversity for this pair of genotypes
									for ($bb = 0 ; $bb < $#H2_first_allelez; $bb++) {
										for ($aa = ($bb+1); $aa <= $#H2_first_allelez; $aa++) {
											if($H2_first_allelez[$bb] ne $H2_first_allelez[$aa]){
												$diffH2+=1; 
											}
											$num_comparisonsH2+=1;
										}
									}		
								}
							}
						}					
					}
					# for some sites, there may be only one individual with a genotype
					if($num_comparisonsH2==0){
						for ($y = 0 ; $y <= $#H2; $y++ ) {
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'C/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'G/*')
								&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)){
								@H2_first_allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
								if($H2_first_allelez[0] ne $H2_first_allelez[1]){
									$diffH2+=1;
								}
								$num_comparisonsH2+=1;	
							}
						}	
					}	
					# now tabulate the average pairwise nucleotide diversity within H2
					if($num_comparisonsH2>0){
							# this does not depend on the position being polymorphic
								$H2_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH2/$num_comparisonsH2);
								$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
								# later we will standardize this by the number of sites per window
					}
					# Now calculate the average pairwise nucleotide diversity within H1
					$diffH1=0;
					$num_comparisonsH1=0;
					@H1_first_allelez=();
					@H1_second_allelez=();
					for ($y = 0 ; $y < $#H1; $y++ ) {
						for ($w = ($y+1) ; $w <= $#H1; $w++ ) {
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'G/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3) 
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'G/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'C/*')
								&&(length($temp[$H1[$w]+$whotoinclude[1]-2]) == 3)){
								# there are data for both H1 alleles
								@H1_first_allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
								@H1_second_allelez=split('/',$temp[$H1[$w]+$whotoinclude[1]-2]);
								# combine the alleles into one array
								@H1_first_allelez = (@H1_first_allelez, @H1_second_allelez);
								# check that the array has 4 elements
								if($#H1_first_allelez != 3){
									print "Problem with number of alleles @H1_first_allelez @H1_second_allelez\n";
								}
								else{
									for ($bb = 0 ; $bb < $#H1_first_allelez; $bb++) {
										for ($aa = ($bb+1); $aa <= $#H1_first_allelez; $aa++) {
											if($H1_first_allelez[$bb] ne $H1_first_allelez[$aa]){
												$diffH1+=1; 
											}
											$num_comparisonsH1+=1;
										}
									}		
								}
							}
						}					
					}
					# for some sites, there may be only one individual with a genotype
					if($num_comparisonsH1==0){
						for ($y = 0 ; $y <= $#H1; $y++ ) {
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'C/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'G/*')
								&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3)){
								@H1_first_allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
								if($H1_first_allelez[0] ne $H1_first_allelez[1]){
									$diffH1+=1;
								}
								$num_comparisonsH1+=1;	
							}
						}	
					}	
					# now tabulate the average pairwise nucleotide diversity within H1
					if($num_comparisonsH1>0){
							# this does not depend on the position being polymorphic
								$H1_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH1/$num_comparisonsH1);
								$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
								# later we will standardize this by the number of sites per window
					}
					# Now calculate the average pairwise nucleotide diversity within H3
					$diffH3=0;
					$num_comparisonsH3=0;
					@H3_first_allelez=();
					@H3_second_allelez=();
					for ($y = 0 ; $y < $#H3; $y++ ) {
						for ($w = ($y+1) ; $w <= $#H3; $w++ ) {
							if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'G/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3) 
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne 'G/*')
								&&($temp[$H3[$w]+$whotoinclude[1]-2] ne 'C/*')
								&&(length($temp[$H3[$w]+$whotoinclude[1]-2]) == 3)){
								# there are data for both H3 alleles
								@H3_first_allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
								@H3_second_allelez=split('/',$temp[$H3[$w]+$whotoinclude[1]-2]);
								# combine the alleles into one array
								@H3_first_allelez = (@H3_first_allelez, @H3_second_allelez);
								# check that the array has 4 elements
								if($#H3_first_allelez != 3){
									print "Problem with number of alleles @H3_first_allelez @H3_second_allelez\n";
								}
								else{
									for ($bb = 0 ; $bb < $#H3_first_allelez; $bb++) {
										for ($aa = ($bb+1); $aa <= $#H3_first_allelez; $aa++) {
											if($H3_first_allelez[$bb] ne $H3_first_allelez[$aa]){
												$diffH3+=1; 
											}
											$num_comparisonsH3+=1;
										}
									}		
								}
							}
						}					
					}
					# for some sites, there may be only one individual with a genotype
					if($num_comparisonsH3==0){
						for ($y = 0 ; $y <= $#H3; $y++ ) {
							if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'T/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'C/*')
								&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'G/*')
								&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)){
								@H3_first_allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
								if($H3_first_allelez[0] ne $H3_first_allelez[1]){
									$diffH3+=1;
								}
								$num_comparisonsH3+=1;	
							}
						}	
					}	
					# now tabulate the average pairwise nucleotide diversity within H3
					if($num_comparisonsH3>0){
							# this does not depend on the position being polymorphic
								$H3_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH3/$num_comparisonsH3);
								$H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
								# later we will standardize this by the number of sites per window
					}					
					# calculate the ABBA BABBA stats, plus more, but only for positions with only two variants
					@temp1=split('',$string);
					$x_uniq = uniq @temp1;
					if($x_uniq == 2){
						$BBAA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*$H2_derived_freq*(1-$H3_derived_freq));
						if((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)>0) || (($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)>0)){
							# this is a polymorphic position and it is an ABBA_BABA site   
							#if($current_window ==10000000 ){
							#	print $H1_derived_freq,"\t",$H2_derived_freq,"\t",$H3_derived_freq,"\n";
							#}

							$ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq);
							$BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq);
							
							# define the doner population for the comparison between H2 and H3
							if($H2_derived_freq > $H3_derived_freq){ # H3 is the doner
								$peak_H2_H3_derived_freq = $H2_derived_freq;
							}
							else{ # H3 is the doner
								$peak_H2_H3_derived_freq = $H3_derived_freq;	
							}
							# The difference between these two will be the doner difference, which is the denominator if the site is an ABBA site
							# i.e. if $ABBA_hash{} > $BABA_hash{} 
							$ABBA_peak_hash=((1-$H1_derived_freq)*$peak_H2_H3_derived_freq*$peak_H2_H3_derived_freq);
							$BABA_peak_hash=($H1_derived_freq*(1-$peak_H2_H3_derived_freq)*$peak_H2_H3_derived_freq);
							# here we are calculating stats for f for H2 and H3.
							# we need to do this for each site because some of them can be undefined and should
							# not be included in the average
							if(($ABBA_peak_hash - $BABA_peak_hash) != 0){
								$f_H2_H3{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ ($ABBA_peak_hash - $BABA_peak_hash);
								$f_H2_H3_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}	
							# here we are calculating stats for f and the assignment of H1 and H2 needs to be switched.
							if($H1_derived_freq > $H3_derived_freq){
								$peak_H1_H3_derived_freq = $H1_derived_freq;
							}
							else{
								$peak_H1_H3_derived_freq = $H3_derived_freq;	
							}
							$BABA_peak_hashH1H3=($peak_H1_H3_derived_freq*(1-$H2_derived_freq)*$peak_H1_H3_derived_freq);
							$ABBA_peak_hashH1H3=((1-$peak_H1_H3_derived_freq)*$H2_derived_freq*$peak_H1_H3_derived_freq);
							if(($BABA_peak_hashH1H3 - $ABBA_peak_hashH1H3) != 0){
								$f_H1_H3{$window_counter."_".$current_chromosome."_".$current_window} += ((($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq))-(((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)))/ ($BABA_peak_hashH1H3 - $ABBA_peak_hashH1H3);
								$f_H1_H3_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}
							if(($H2_derived_freq >= $H1_derived_freq)&&(($ABBA_peak_hash - $BABA_peak_hash)!=0)){
								$f_dm{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ ($ABBA_peak_hash - $BABA_peak_hash);
								# only count the sites that contribute to f_dm
								if(((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))!=0){
									$f_dm_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
								}															
							}
							elsif(($ABBA_peak_hashH1H3 - $BABA_peak_hashH1H3)!=0){
								$f_dm{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ -($ABBA_peak_hashH1H3 - $BABA_peak_hashH1H3);	
								# only count the sites that contribute to f_dm
								if ( ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ -($ABBA_peak_hashH1H3 - $BABA_peak_hashH1H3) != 0){
									$f_dm_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
								}
							}
	
						}
					}
				}	
			}
		}
	}
}


close DATAINPUT;

# now merge the hash keys
my @common_keys = ();

foreach (keys %ABBA_hash) {
	push(@common_keys, $_);
}

foreach (keys %BABA_hash) {
	push(@common_keys, $_) unless exists $ABBA_hash{$_};
}

my @common_keys2 = ();

foreach (keys %ABBA_hash) {
	push(@common_keys2, $_);
}

foreach (keys %BBAA_hash) {
	push(@common_keys2, $_) unless exists $ABBA_hash{$_};
}


my @out = keys %{{map {($_ => 1)} (@common_keys, @common_keys2)}};

@out = map  { $_->[0] }
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
                 @out;



foreach (@out) {
	@temp1=split('_',$_);
	print OUTFILE $temp1[1],"\t",$temp1[2]+1,"\t",$temp1[2]+$sliding_window,"\t";
	if(defined($ABBA_hash{$_})){
		print OUTFILE $ABBA_hash{$_},"\t";
	}
	else{
		print OUTFILE "0\t";
	}
	if(defined($BABA_hash{$_})){
		print OUTFILE $BABA_hash{$_},"\t0\t0\t0\t0\n";
	}
	else{
		print OUTFILE "0\t0\t0\t0\t0\n";
	}
}
if($#out == -1){
	print OUTFILE "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE;


print OUTFILE2 "chromosome\tbegin\tend\tABBA\tBABA\tBBAA\tD\tfd_H2H3\tfd_H1H3\tf_dm\tnumsites_f_dm\tdH2H3\tnum_sites_per_windowH2H3\tdH1H3\tnum_sites_per_windowH1H3\tdH1H2\tnum_sites_per_windowH1H2\tH3pi\tnumsitesH3pi\tH2pi\tnumsitesH2pi\tH1pi\tnumsitesH1pi\n";
foreach (@out) {
	@temp1=split('_',$_);
	print OUTFILE2 $temp1[1],"\t",$temp1[2]+1,"\t",$temp1[2]+$sliding_window,"\t";
	if(defined($ABBA_hash{$_})){
		print OUTFILE2 $ABBA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	if(defined($BABA_hash{$_})){
		print OUTFILE2 $BABA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	if(defined($BBAA_hash{$_})){
		print OUTFILE2 $BBAA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	#print D for this window
	if(($ABBA_hash{$_}+$BABA_hash{$_})>0){
		print OUTFILE2 ($ABBA_hash{$_}-$BABA_hash{$_}) / ($ABBA_hash{$_}+$BABA_hash{$_}),"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}

	#print fd H2H3 for this window
	if(defined($f_H2_H3_counter{$_})){
		print OUTFILE2 $f_H2_H3{$_}/$f_H2_H3_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}
	#print fd H1H3 for this window
	if(defined($f_H1_H3_counter{$_})){
		print OUTFILE2 $f_H1_H3{$_}/$f_H1_H3_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}	
	#print f_dm for this window
	if(defined($f_dm{$_})){
		print OUTFILE2 $f_dm{$_}/$f_dm_counter{$_},"\t",$f_dm_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}	
	#print average H2H3 pairwise divergence for this window and number of H2H3 sites
	if(defined($number_of_sites_per_window{$_})){
		if($number_of_sites_per_window{$_}>0){
			print OUTFILE2 ($H2_H3_pairwise_divergence_per_window{$_}/$number_of_sites_per_window{$_}),"\t",$number_of_sites_per_window{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$number_of_sites_per_window{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$number_of_sites_per_window{$_},"\t";
	}
	#print average H1H3 pairwise divergence for this window and number of H1H3 sites
	if(defined($number_of_sites_per_windowH1H3{$_})){
		if($number_of_sites_per_windowH1H3{$_} > 0){
			print OUTFILE2 ($H1_H3_pairwise_divergence_per_window{$_}/$number_of_sites_per_windowH1H3{$_}),"\t",$number_of_sites_per_windowH1H3{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$number_of_sites_per_windowH1H3{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$number_of_sites_per_windowH1H3{$_},"\t";
	}
	#print average H1H2 pairwise divergence for this window and number of H1H2 sites
	if(defined($number_of_sites_per_windowH1H2{$_})){
		if($number_of_sites_per_windowH1H2{$_} > 0){
			print OUTFILE2 ($H1_H2_pairwise_divergence_per_window{$_}/$number_of_sites_per_windowH1H2{$_}),"\t",$number_of_sites_per_windowH1H2{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$number_of_sites_per_windowH1H2{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$number_of_sites_per_windowH1H2{$_},"\t";
	}
	#print H3 pairwise nucleotide diversity per site for this window and number of H3 pairwise nucleotide diversity sites
	if(defined($H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_})){
		if($H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}>0){
			print OUTFILE2 ($H3_pairwise_nucleotide_diversity_per_window{$_}/$H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}),"\t",$H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$H3_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
	}
	#print H2 pairwise nucleotide diversity per site for this window and number of H2 pairwise nucleotide diversity sites
	if(defined($H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_})){
		if($H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}>0){
			print OUTFILE2 ($H2_pairwise_nucleotide_diversity_per_window{$_}/$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}),"\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
	}
	#print H1 pairwise nucleotide diversity per site for this window and number of H1 pairwise nucleotide diversity sites
	if(defined($H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_})){
		if($H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}>0){
			print OUTFILE2 ($H1_pairwise_nucleotide_diversity_per_window{$_}/$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}),"\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
		}
		else{
			print OUTFILE2 "NAN\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
	}
}
if($#out == -1){
	print OUTFILE2 "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE2;

```
