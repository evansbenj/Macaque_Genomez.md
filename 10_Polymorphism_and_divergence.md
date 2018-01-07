# Polymorphism and divergence

I have a new in progress script called "Boot_from_tab_diverge_poly_2018_allowmissingdata.pl" that calculates polymorphism and divergence statistics from a tab delimited file. 

For autosomes in this format: #CHROM	POS	REF	SAMN03083651	SAMN03264597	SAMN03264605	SAMN03264609	SAMN03264610	SAMN03264614	SAMN03264646	SAMN03264649	SAMN03264650	SAMN
03264651	SAMN03264652	SAMN03264653	SAMN03264658	SAMN03264659	SAMN03264660	SAMN03264661	SAMN03264662	SAMN03264666	SAMN03264667	SAMN03264673
	SAMN03264674	SAMN03264677	SAMN03264678	SAMN03264683	SAMN03264685	SAMN03264696	SAMN03264703	SAMN03264706	SAMN03264707	SAMN03264708	SAMN
03264709	SAMN03264710	SAMN03264712	SAMN03264713	SAMN03264715	SAMN03264718	SAMN03264719	SAMN03264726	SAMN03264733	SAMN03264734 use this:

```
Boot_from_tab_diverge_poly_2015.pl /work/ben/2017_rhesus_genomez/F_and_M/FandM_chr01_BSQR_jointgeno_allsites_filtered.vcf.gz.tab 1100000100011010101001111110110101001011 3_4_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39_40 rhesus_chr01_poly_and_diverge.txt
```
For all chrs this can be automated with this script: "30_runs_polydiverge.pl"
```
#!/usr/bin/env perl                                                                                                                                                           
use strict;
use warnings;

# this program will read in all tab delimited files in a folder and                                                                                                           
# calculate popgen stats plus boostraps for all species for each                                                                                          
# file                                                                                                                                                                        

my $status;
my @tabfiles = glob("/work/ben/2017_rhesus_genomez/F_and_M/FandM_chr*_BSQR_jointgeno_allsites_filtered.vcf.gz.tab");
my $bin_sex = "1100000100011010101001111110110101001011";
my $commandline;
my @what_what;
my @species = ("rhesus");
my @numberz = ("1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39_40");
my $y;

foreach(@tabfiles){
    @what_what = split(".vcf.gz.",$_);
    for ($y = 0 ; $y <= $#species ; $y++ ) {
        $commandline = "sqsub -r 2d --mpp 16G -o ".$what_what[0].".log bash -c \"perl Boot_from_tab_diverge_poly_2018_allowmissingdata.pl ".$_." ".$bin_sex." 3_4_".$numberz[$y]." ".$what_what[0]."_".$species[$y]."_boot.poly\"";
        print $commandline,"\n";
       $status = system($commandline);                                                                                                                                       
    }
}


```

Here is the preliminary script: "Boot_from_tab_diverge_poly_2018_allowmissingdata.pl":
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;


#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and bootstraps TajD and polymorphism statistics by autosomal DNA 
#  and by xDNA by resampling bases
#  with replacement. 

#	Because TajD does not require an outgroup, the data analyzed will be 
# 	different from other analyzes that require divergence data (such as 
#	pi/D or S/D and also analyses that require outgroup information (such
#	as the analysis of the derived AFS).

#	This analysis will allow for missing data but in doing so assumes that the missing data are randomly distributed.

# 	If this were not the case, for example if a diverged sample had lots of missing data, 
#	then the estimates would be biased (downward biased in that example)

# 	This program will be compatible with files with one or multiple outgroup columns

# 	It will also accomodate data from multiple species and calculate the stats only from selected columns

# to execute type Boot_from_tab_diverge_poly_2018_allowmissingdata.pl inputfile.tab 1111100110000111100011100110010100000000 
# 3_4_22_23_25_26 nigra_poly_and_diverge.txt  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_4_22_23_25_26 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case), (ii) the column number of the first individual in the ingroup 
# (4 in this case), and (iii) the sample number that contain the data from the individuals you want to 
# include (22, 23, 25, and 26 in this case), which are the four nigra samples itemized below.

# IMPORTANT: (i) and (ii) are columns beginning with 1 but (iii) is based on the individual samples such
# as enumerated below

# for example, with a tab file with only the baboon sequence in the 4th column, here is the input command:

# tonk
# Boot_from_tab_diverge_poly_2015.pl final_round2_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_32_33_34_35_36_37_38_39_40 tonk_poly_and_diverge.txt


my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile = $ARGV[3];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";


my @sexes = split("",$ARGV[1]);
my @whotoinclude = split("_",$ARGV[2]);


my @temp;
my @temp1;
my $previous= 0;
my @aDNA_segregating_sites=();
my @aDNA_singleton_sites=();
my @aDNA_sites=();
my $aDNA_total_number_of_sites=0; # this is the total number of sites with at least two alleles genotyped.
my @xDNA_segregating_sites=();
my @xDNA_singleton_sites=();
my @xDNA_sites=();
my $xDNA_total_number_of_sites=0; # this is the total number of sites with at least two alleles genotyped.
my @yDNA_segregating_sites=();
my @yDNA_singleton_sites=();
my @yDNA_sites=();
my $yDNA_total_number_of_sites=0; # this is the total number of sites with at least two alleles genotyped.

# variables for bootstrapping
my @d_aDNA=();
my @d_xDNA=();
my @d_yDNA=();
my $pi_aDNA=0;
my $pi_xDNA=0;
my $pi_yDNA=0;

my @pi_aDNA=();
my @pi_xDNA=();
my @pi_yDNA=();
my %S_aDNA=();
my %S_xDNA=();
my %S_yDNA=();

my $aDNA_divergence=0;
my $xDNA_divergence=0;
my $yDNA_divergence=0;
my $string;
my $m;
my $n;
my $number_of_bootstraps=0;
my $lower=int($number_of_bootstraps*0.025);
my $upper=int($number_of_bootstraps*0.975);
my $JC_divergence_aDNA;
my $JC_divergence_xDNA;
my $JC_divergence_yDNA;
my $RAD_tag_count_aDNA=0;
my $RAD_tag_count_xDNA=0;
my $RAD_tag_count_yDNA=0;
my $distance_between_RAD_tags=500;
my $previousone_aDNA=0-$distance_between_RAD_tags;
my $previousone_xDNA=0-$distance_between_RAD_tags;
my $previousone_yDNA=0-$distance_between_RAD_tags;

print "bootlower ", $lower," bootupper ",$upper,"\n";

my $w;
my $y;
my $x;
my @unique;

my $number_of_individuals_genotyped=($#whotoinclude - 1);

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped=0;
for ($y = 2 ; $y <= $#whotoinclude ; $y++ ) {
	if($sexes[$whotoinclude[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	

print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";


# initialize some arrays; don't bother with zero because this is the slot for only 1 allele and no diversity estiamtes
# are possible for one allele so it is ignored.

for ($y = 1 ; $y < ($number_of_individuals_genotyped*2) ; $y++ ) {
	$aDNA_segregating_sites[$y]=0; # This is the number of segregating sites for positions with genotype sample size $y+1
	$aDNA_singleton_sites[$y]=0; # This is the number of segregating sites that are singletons for positions with genotype sample size $y+1
	$aDNA_sites[$y]=0; # This is the number of sites for positions with genotype sample size $y+1
}	
for ($y = 1 ; $y < (2*$number_of_individuals_genotyped-($number_of_individuals_genotyped-$number_of_female_individuals_genotyped)) ; $y++ ) {
	$xDNA_segregating_sites[$y]=0; # This is the number of segregating sites for positions with genotype sample size $y+1
	$xDNA_singleton_sites[$y]=0; # This is the number of segregating sites that are singletons for positions with genotype sample size $y+1
	$xDNA_sites[$y]=0; # This is the number of sites for positions with genotype sample size $y+1
}
for ($y = 1 ; $y < ($number_of_individuals_genotyped-$number_of_female_individuals_genotyped) ; $y++ ) {
	$yDNA_segregating_sites[$y]=0; # This is the number of segregating sites for positions with genotype sample size $y+1
	$yDNA_singleton_sites[$y]=0; # This is the number of segregating sites that are singletons for positions with genotype sample size $y+1
	$yDNA_sites[$y]=0; # This is the number of sites for positions with genotype sample size $y+1
}

my $x;
my @a_W=();
# calculate the denominator for Watterson's estimator for all possible sample sizes
# the index is the sample size minus one (sample size of 2 has an index of 1)
for ($y = 1 ; $y < ($number_of_individuals_genotyped*2) ; $y++ ) {
	for ($x = 0; $x <= $y; $x++ ) {	
		$a_W[$y]+=(1/($x+1));
	}
}

my $asum=0;
my $xsum=0;
my $ysum=0;
my $asum_baboons=0;
my $xsum_baboons=0;
my $ysum_baboons=0;

my $pi_counter=0;
my $diff=0;
my @slice;
my $x_uniq=0;
my $y_uniq=0;


#print "Num automomal alleles: ",$n_aDNA,"\n";
#print "Num chrX alleles: ",$n_xDNA,"\n";
#print "Num chrY alleles: ",$n_yDNA,"\n";


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		# base the genomic location on the outgroup	
		if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")&&($temp[0] ne "chrM")){
				# load the autosomal data
				$string=();
				for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
					# load the first allele
					if((uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'A')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'T')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'C')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'G')){
						$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
						$string=$string.$w;
					}	
					# now load the second allele
					if((uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'A')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'T')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'C')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'G')){
						$w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
						$string=$string.$w;
					}	
				}
				if(defined($string)){
					if(length($string)>1){ # we need at least two alleles to estimate diversity
						@temp1=split('',$string);
						if(
							(uc $temp[$whotoinclude[0]-1] eq "A")||
							(uc $temp[$whotoinclude[0]-1] eq "C")||
							(uc $temp[$whotoinclude[0]-1] eq "T")||
							(uc $temp[$whotoinclude[0]-1] eq "G")
							){ # the outgroup must be defined
								$aDNA_sites[$#temp1]+=1; 		# this is the number of sites with this sample size that are analyzed
																# and the outgroup defined
								$aDNA_total_number_of_sites+=1; # this is the total number of sites with at least two alleles genotyped
																# and the outgroup defined
								$x_uniq = uniq @temp1; 			# this is the number of unique variants at this position
								if($x_uniq == 1){
									push(@pi_aDNA,0); 				# this will be used later for bootstrapping
									push(@{$S_aDNA{$#temp1}},0); 	# this will be used later for bootstrapping, stores number of 
																	# seg sites at each position for each sample size
									# invariant positions do not contribute to pi
								}
								elsif($x_uniq > 1){
									$aDNA_segregating_sites[$#temp1]+=1; # this is a segregating site; counted once if there are 2 or more SNPs
									push(@{$S_aDNA{$#temp1}},1); # this will be used later for bootstrapping, stores number of 
																 # seg sites at each position for each sample size
									# now calculate pi
									$diff=0; # this is for pairwise comparisons within a site
									for ($y = 0 ; $y < $#temp1 ; $y++ ) {
										for ($x = ($y+1) ; $x <= $#temp1 ; $x++ ) {
											if($temp1[$y] ne $temp1[$x]){
												$diff+=1;
											}
											$pi_counter+=1;
										}
									}
									$pi_aDNA+=$diff/$pi_counter; 	# this will be standardized later by dividing by the total number of sites
																	# $aDNA_total_number_of_sites
									push(@pi_aDNA,$diff/$pi_counter); # this will be used later for bootstrapping
									if($diff == $#temp1){ # there is a singleton SNP
										$aDNA_singleton_sites[$#temp1]+=1;
										# push(@{$aDNA_singleton_sites{$#temp1},1); #This could be used later for bootstrapping
									}
									#else{
									#	push(@{$aDNA_singleton_sites{$#temp1},0);#This could be used later for bootstrapping
									#}
									$diff=0;			
									$pi_counter=0;	
								}
								# now calculate divergence
								if(uc $temp[$whotoinclude[0]-1] ne uc $temp1[0]){ # divergence is calculated based on the first genotyped ingroup allele
									$aDNA_divergence+=1;
									push(@d_aDNA,1);
								}
								else{
									push(@d_aDNA,0);
								}

						}
					}
				}	
		}
		elsif($temp[0] eq "chrX"){
			$string=();
			# for chrX, load both female alleles but only one male allele
			# for each column, we need to check if the individual is a female
			# first the female
			for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
				# load both alleles if the individual is a female
				if($sexes[$whotoinclude[$y+2]-1] eq "1"){
					# load the first allele
					if((uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'A')||
					(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'T')||
					(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'G')||
					(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'C')){
						$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
						$string=$string.$w;
					}	
					# now load the second allele
					if((uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'A')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'T')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'G')||
						(uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] eq 'C')){
						$w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
						$string=$string.$w;
					}	
				}
				# load one allele if the individual is a male
				elsif($sexes[$whotoinclude[$y+2]-1] eq "0"){
					# load only the first allele
					if(($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'A')||
						($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'T')||
						($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'G')||
						($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'C')){
						$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
						$string=$string.$w;
					}	
				}
				else{
					print "Something is wrong with figuring out what sex each individual is X ",$sexes[$whotoinclude[$y+2]-1]," ",$whotoinclude[$y+2],"\n";
				}
			} # end of cycling through each individual for chrX	
			if(defined($string)){
				if(length($string)>1){ # we need at least two alleles to estimate diversity
					@temp1=split('',$string);			
					if(
					    (uc $temp[$whotoinclude[0]-1] eq "A")||
						(uc $temp[$whotoinclude[0]-1] eq "C")||
						(uc $temp[$whotoinclude[0]-1] eq "T")||
						(uc $temp[$whotoinclude[0]-1] eq "G")
						){ # the outgroup must be defined
							$xDNA_sites[$#temp1]+=1; 				# this is the number of sites with this sample size that are analyzed
																	# and the outgroup defined
							$xDNA_total_number_of_sites+=1; 		# this is the total number of sites with at least two alleles genotyped
																	# and the outgroup defined
							$x_uniq = uniq @temp1;					# this is the number of unique variants at this position	
							if($x_uniq == 1){
								push(@pi_xDNA,0);					# this will be used later for bootstrapping
								push(@{$S_xDNA{$#temp1}},0);		# this will be used later for bootstrapping, stores number of 
																	# seg sites at each position for each sample size
								# invariant positions do not contribute to pi									
							}
							elsif($x_uniq > 1){
								$xDNA_segregating_sites[$#temp1]+=1; 	# this is a segregating site; counted once if there are 2 or more SNPs
								push(@{$S_xDNA{$#temp1}},1); 			# this will be used later for bootstrapping, stores number of 
																		# seg sites at each position for each sample size

								# now calculate pi
								$diff=0;
								for ($y = 0 ; $y < $#temp1 ; $y++ ) {
									for ($x = ($y+1) ; $x <= $#temp1 ; $x++ ) {
										if($temp1[$y] ne $temp1[$x]){
											$diff+=1;
										}
										$pi_counter+=1;
									}
								}
								$pi_xDNA+=$diff/$pi_counter; 			# this will be standardized later by dividing by the total number of sites
																		# $xDNA_total_number_of_sites
								push(@pi_xDNA,$diff/$pi_counter);
								if($diff == $#temp1){
									$xDNA_singleton_sites[$#temp1]+=1;
									# push(@{$xDNA_singleton_sites{$#temp1}},1); #This could be used later for bootstrapping
								}
								#else{
								#	push(@{$xDNA_singleton_sites{$#temp1}},0);
								#}
								$diff=0;			
								$pi_counter=0;	
							}
							# now calculate divergence
							if(uc $temp[$whotoinclude[0]-1] ne uc $temp1[0]){
								$xDNA_divergence+=1;
								push(@d_xDNA,1);
							}
							else{
								push(@d_xDNA,0);
							}	
					}
				}
			}
		} # endifelse to check for aDNA, chrX	
		elsif($temp[0] eq "chrY"){
			$string=();
			# for chrX, load both female alleles but only one male allele
			# for each column, we need to check if the individual is a female
			# first the female
			for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
				# load one allele if the individual is a male
				if($sexes[$whotoinclude[$y+2]-1] eq "0"){
					# load only the first allele
					if(($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'A')||
						($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'T')||
						($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'G')||
						($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] eq 'C')){
						$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
						$string=$string.$w;
					}	
				}
				else{
					print "Something is wrong with figuring out what sex each individual is X ",$sexes[$whotoinclude[$y+2]-1]," ",$whotoinclude[$y+2],"\n";
				}
			} # end of cycling through each individual for chrY	
			if(defined($string)){
				if(length($string)>1){ # we need at least two alleles to estimate diversity
					@temp1=split('',$string);			
					if(
					    (uc $temp[$whotoinclude[0]-1] eq "A")||
						(uc $temp[$whotoinclude[0]-1] eq "C")||
						(uc $temp[$whotoinclude[0]-1] eq "T")||
						(uc $temp[$whotoinclude[0]-1] eq "G")
						){ # the outgroup must be defined
							$yDNA_sites[$#temp1]+=1; 				# this is the number of sites with this sample size that are analyzed
																	# and the outgroup defined
							$yDNA_total_number_of_sites+=1; 		# this is the total number of sites with at least two alleles genotyped
																	# and the outgroup defined
							$x_uniq = uniq @temp1;	# the number of unique variants at this position	
							if($x_uniq == 1){
								push(@pi_yDNA,0);
								push(@{$S_yDNA{$#temp1}},0);		# this will be used later for bootstrapping, stores number of 
																	# seg sites at each position for each sample size
								# invariant positions do not contribute to pi
							}
							elsif($x_uniq > 1){
								$yDNA_segregating_sites[$#temp1]+=1; # this is a segregating site; counted once if there are 2 or more SNPs
								push(@{$S_yDNA{$#temp1}},1); 		# this will be used later for bootstrapping, stores number of 
																	# seg sites at each position for each sample size
								# now calculate pi
								$diff=0;
								for ($y = 0 ; $y < $#temp1 ; $y++ ) {
									for ($x = ($y+1) ; $x <= $#temp1 ; $x++ ) {
										if($temp1[$y] ne $temp1[$x]){
											$diff+=1;
										}
										$pi_counter+=1;
									}
								}
								$pi_yDNA+=$diff/$pi_counter; 	# this will be standardized later by dividing by the total number of sites
																# $yDNA_total_number_of_sites
								push(@pi_yDNA,$diff/$pi_counter); # this will be used later for bootstrapping
								if($diff == $#temp1){
									$yDNA_singleton_sites[$#temp1]+=1;
									# push(@{$yDNA_singleton_sites{$#temp1}},1); #This could be used later for bootstrapping
								}
								#else{
								#	push(@{$yDNA_singleton_sites{$#temp1}},0);
								#}
								$diff=0;		
								$pi_counter=0;	
							}
							if(uc $temp[$whotoinclude[0]-1] ne uc $temp1[0]){
								$yDNA_divergence+=1;
								push(@d_yDNA,1);
							}
							else{
								push(@d_yDNA,0);
							}	
					}
				}
			}
		} # endifelse to check for aDNA, chrX, chrY
	}		
	elsif($temp[0] eq '#CHROM'){
		for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
			print "Individual ",$temp[$whotoinclude[$y+2]+$whotoinclude[1]-2]," is a ";
			if($sexes[$whotoinclude[$y+2]-1] == 1){
				print "female\n";
			} 
			elsif($sexes[$whotoinclude[$y+2]-1] == 0){
				print "male\n";
			} 
		}
	}	
} # end while

my $aDNA_weighted_thetaW=0;
my $xDNA_weighted_thetaW=0;
my $yDNA_weighted_thetaW=0;

if($#pi_aDNA>0){
	# Standardize the pi value
	$pi_aDNA=$pi_aDNA/$aDNA_total_number_of_sites;
	print OUTFILE "Segregating sites\n";
	print OUTFILE "#alleles\t#sites\t#segsites\tthetaW\n";
	for ($y = 1 ; $y < ($number_of_individuals_genotyped*2) ; $y++ ) {
		print OUTFILE $y+1,"\t",$aDNA_sites[$y],"\t",$aDNA_segregating_sites[$y],"\t",$aDNA_segregating_sites[$y]/$a_W[$y],"\n";
		$aDNA_weighted_thetaW+=$aDNA_sites[$y]*$aDNA_segregating_sites[$y]/$a_W[$y];
	}
	print OUTFILE "aDNA_weighted_thetaW\t",sprintf("%.5f",$aDNA_weighted_thetaW/$aDNA_total_number_of_sites),"\n";
	print OUTFILE "pi_aDNA\t",sprintf("%.5f",$pi_aDNA),"\n";
	print OUTFILE "Number of aDNA sites ",$aDNA_total_number_of_sites,"\n";
}
if($#pi_xDNA>0){
	# Standardize the pi value
	$pi_xDNA=$pi_xDNA/$xDNA_total_number_of_sites;
	print OUTFILE "Segregating sites\n";
	print OUTFILE "#alleles\t#sites\t#segsites\tthetaW\n";
	for ($y = 1 ; $y < (2*$number_of_individuals_genotyped-($number_of_individuals_genotyped-$number_of_female_individuals_genotyped)) ; $y++ ) {
		print OUTFILE $y+1,"\t",$xDNA_sites[$y],"\t",$xDNA_segregating_sites[$y],"\t",$xDNA_segregating_sites[$y]/$a_W[$y],"\n";
		$xDNA_weighted_thetaW+=$xDNA_sites[$y]*$xDNA_segregating_sites[$y]/$a_W[$y];
	}
	print OUTFILE "xDNA_weighted_thetaW\t",sprintf("%.5f",$xDNA_weighted_thetaW/$xDNA_total_number_of_sites),"\n";
	print OUTFILE "pi_xDNA\t",sprintf("%.5f",$pi_xDNA),"\n";
	print OUTFILE "Number of xDNA sites ",$xDNA_total_number_of_sites,"\n";
}
if($#pi_yDNA>0){
	# Standardize the pi value
	$pi_yDNA=$pi_yDNA/$yDNA_total_number_of_sites;
	print OUTFILE "Segregating sites\n";
	print OUTFILE "#alleles\t#sites\t#segsites\tthetaW\n";
	for ($y = 1 ; $y < ($number_of_individuals_genotyped-$number_of_female_individuals_genotyped) ; $y++ ) {
		print OUTFILE $y+1,"\t",$yDNA_sites[$y],"\t",$yDNA_segregating_sites[$y],"\t",$yDNA_segregating_sites[$y]/$a_W[$y],"\n";
		$yDNA_weighted_thetaW+=$yDNA_sites[$y]*$yDNA_segregating_sites[$y]/$a_W[$y];
	}
	print OUTFILE "yDNA_weighted_thetaW\t",sprintf("%.5f",$yDNA_weighted_thetaW/$yDNA_total_number_of_sites),"\n";
	print OUTFILE "pi_xDNA\t",sprintf("%.5f",$pi_yDNA),"\n";
	print OUTFILE "Number of yDNA sites ",$yDNA_total_number_of_sites,"\n";
}



=pod # beginning of a commented block
print $aDNA_sites," ",($#pi_aDNA+1)," ",($#S_aDNA+1),"\n";
print $xDNA_sites," ",($#pi_xDNA+1)," ",($#S_xDNA+1),"\n";
print $yDNA_sites," ",($#pi_yDNA+1)," ",($#S_yDNA+1),"\n";



my @aDNA_bootstrapped_indexes=();
my @xDNA_bootstrapped_indexes=();
# no bootstrap for the y - use coalescent simulations

my @S_aDNA_boot;
my @pi_aDNA_boot;
my @pi_aDNA_persite_boot;
my @d_aDNA_boot;

my @S_xDNA_boot;
my @pi_xDNA_boot;
my @pi_xDNA_persite_boot;
my @d_xDNA_boot;

my @aDNA_TajD_boot;
my @xDNA_TajD_boot;

my @JC_AP_divergence_aDNA_boot;
my @JC_AP_divergence_xDNA_boot;

my @thetapi_over_divergence_aDNA;
my @thetapi_over_divergence_xDNA;

my @aDNA_singleton_sites_boot;
my @xDNA_singleton_sites_boot;
my @aDNA_singleton_sites_indexes;
my @xDNA_singleton_sites_indexes;

for ($m = 0 ; $m < $number_of_bootstraps ; $m++ ) {
	# generate an array with bootstrapped indexes
	# first aDNA
	print "bootstrap ",$m,"\n";
	
	for ($n = 0 ; $n < $aDNA_sites ; $n++ ) {
		push(@aDNA_bootstrapped_indexes,int(rand($aDNA_sites)));
	} # end $n

	# same for segregating sites
	for ($n = 0 ; $n < $aDNA_segregating_sites ; $n++ ) {
		push(@aDNA_singleton_sites_indexes,int(rand($aDNA_segregating_sites)));
	} # end $n

	for ($n = 0 ; $n < $aDNA_sites ; $n++ ) {
		# calculate segregating sites and pi and d for bootstrap replicates
		$S_aDNA_boot[$m]+=$S_aDNA[$aDNA_bootstrapped_indexes[$n]];
		$pi_aDNA_boot[$m]+=$pi_aDNA[$aDNA_bootstrapped_indexes[$n]];
		$d_aDNA_boot[$m]+=$d_aDNA[$aDNA_bootstrapped_indexes[$n]];
	} # end $n
	
	push(@pi_aDNA_persite_boot,$pi_aDNA_boot[$m]/$aDNA_sites);


	# apply JC correction to divergence
	$JC_AP_divergence_aDNA_boot[$m]  = (-3/4)*log(1-(4/3)*($d_aDNA_boot[$m]/$aDNA_sites));
	# apply AP correction to divergence
	$JC_AP_divergence_aDNA_boot[$m]  = $JC_AP_divergence_aDNA_boot[$m] - ($pi_aDNA_boot[$m]/$aDNA_sites);

	# calculate pi/divergence for bootstrapped data
	push(@thetapi_over_divergence_aDNA,($pi_aDNA_boot[$m]/$aDNA_sites)/$JC_AP_divergence_aDNA_boot[$m]);

	# calculate TajD for bootstrapped data

	if((($e1_obs_aDNA*$S_aDNA_boot[$m])+($e2_obs_aDNA*$S_aDNA_boot[$m]*($S_aDNA_boot[$m]-1))) > 0){
		push(@aDNA_TajD_boot,($pi_aDNA_boot[$m] - ($S_aDNA_boot[$m]/$a1_obs_aDNA))/((($e1_obs_aDNA*$S_aDNA_boot[$m])+($e2_obs_aDNA*$S_aDNA_boot[$m]*($S_aDNA_boot[$m]-1)))**0.5));
 	}

	# singleton bootstrap
 	for ($n = 0 ; $n < $aDNA_segregating_sites ; $n++ ) {
		#calculate number of singletons for each replicate
		$aDNA_singleton_sites_boot[$m]+=$aDNA_singleton_sites[$aDNA_singleton_sites_indexes[$n]]
	} # end $n
	
	# now convert to a proportion
	$aDNA_singleton_sites_boot[$m]=($aDNA_singleton_sites_boot[$m]/$aDNA_segregating_sites);
  

	# now do xDNA
	for ($n = 0 ; $n < $xDNA_sites ; $n++ ) {
		push(@xDNA_bootstrapped_indexes,int(rand($xDNA_sites)));
	} # end $n

	for ($n = 0 ; $n < $xDNA_sites ; $n++ ) {
		# calculate segregating sites and pi for bootstrap replicates
		$S_xDNA_boot[$m]+=$S_xDNA[$xDNA_bootstrapped_indexes[$n]];
		$pi_xDNA_boot[$m]+=$pi_xDNA[$xDNA_bootstrapped_indexes[$n]];
		$d_xDNA_boot[$m]+=$d_xDNA[$xDNA_bootstrapped_indexes[$n]];
	} # end $n

	# same for segregating sites X
	for ($n = 0 ; $n < $xDNA_segregating_sites ; $n++ ) {
		push(@xDNA_singleton_sites_indexes,int(rand($xDNA_segregating_sites)));
	} # end $n

	if($xDNA_sites > 0){
		push(@pi_xDNA_persite_boot,$pi_xDNA_boot[$m]/$xDNA_sites);
	}
	else{
		push(@pi_xDNA_persite_boot,0);
	}
	# apply JC correction to divergence
	$JC_AP_divergence_xDNA_boot[$m]  = (-3/4)*log(1-(4/3)*($d_xDNA_boot[$m]/$xDNA_sites));
	# apply AP correction to divergence
	$JC_AP_divergence_xDNA_boot[$m]  = $JC_AP_divergence_xDNA_boot[$m] - ($pi_xDNA_boot[$m]/$xDNA_sites);

	# calculate pi/divergence for bootstrapped data
	push(@thetapi_over_divergence_xDNA,($pi_xDNA_boot[$m]/$xDNA_sites)/$JC_AP_divergence_xDNA_boot[$m]);
	# calculate TajD for bootstrapped data
	if((($e1_obs_xDNA*$S_xDNA_boot[$m])+($e2_obs_xDNA*$S_xDNA_boot[$m]*($S_xDNA_boot[$m]-1))) > 0){
		push(@xDNA_TajD_boot,($pi_xDNA_boot[$m] - ($S_xDNA_boot[$m]/$a1_obs_xDNA))/((($e1_obs_xDNA*$S_xDNA_boot[$m])+($e2_obs_xDNA*$S_xDNA_boot[$m]*($S_xDNA_boot[$m]-1)))**0.5));
	}

	# singleton bootstrap
 	for ($n = 0 ; $n < $xDNA_segregating_sites ; $n++ ) {
		#calculate number of singletons for each replicate
		$xDNA_singleton_sites_boot[$m]+=$xDNA_singleton_sites[$xDNA_singleton_sites_indexes[$n]]
	} # end $n
	
	# now convert to a proportion
	$xDNA_singleton_sites_boot[$m]=($xDNA_singleton_sites_boot[$m]/$xDNA_segregating_sites);


	# reset the indexes for next bootstrap
	@aDNA_bootstrapped_indexes=();
	@xDNA_bootstrapped_indexes=();
	@aDNA_singleton_sites_indexes=();
	@xDNA_singleton_sites_indexes=();



} # end $m bootstraps


#print "aDNA_boot @aDNA_TajD_boot \n";
#print "xDNA_boot @xDNA_TajD_boot \n";


# now get 95% CIs
@d_aDNA_boot = sort { $a <=> $b } @d_aDNA_boot;
@d_xDNA_boot = sort { $a <=> $b } @d_xDNA_boot;

@S_aDNA_boot = sort { $a <=> $b } @S_aDNA_boot;
@S_xDNA_boot = sort { $a <=> $b } @S_xDNA_boot;

@pi_aDNA_persite_boot = sort { $a <=> $b } @pi_aDNA_persite_boot;
@pi_xDNA_persite_boot = sort { $a <=> $b } @pi_xDNA_persite_boot;

@thetapi_over_divergence_aDNA = sort { $a <=> $b } @thetapi_over_divergence_aDNA;
@thetapi_over_divergence_xDNA = sort { $a <=> $b } @thetapi_over_divergence_xDNA;

@aDNA_TajD_boot = sort { $a <=> $b } @aDNA_TajD_boot;
@xDNA_TajD_boot = sort { $a <=> $b } @xDNA_TajD_boot;

@JC_AP_divergence_aDNA_boot  = sort { $a <=> $b } @JC_AP_divergence_aDNA_boot;
@JC_AP_divergence_xDNA_boot  = sort { $a <=> $b } @JC_AP_divergence_xDNA_boot;


@aDNA_singleton_sites_boot = sort { $a <=> $b } @aDNA_singleton_sites_boot;
@xDNA_singleton_sites_boot = sort { $a <=> $b } @xDNA_singleton_sites_boot;

my $pi_a=0;
my $pi_x=0;
my $pi_y=0;

if($aDNA_sites>0){
	print OUTFILE "aDNA\n";
	print OUTFILE "#_alleles\t",$n_aDNA,"\n";
	print OUTFILE "#_Sites\t",$aDNA_sites,"\n";
	print OUTFILE "RAD_tag_count\t",$RAD_tag_count_aDNA,"\n";
	print OUTFILE "S\t",$aDNA_segregating_sites," (",$S_aDNA_boot[$lower]," - ",$S_aDNA_boot[$upper],")\n";
	print OUTFILE "thetaW\t",sprintf("%.5f",$aDNA_segregating_sites/$a1_obs_aDNA/$aDNA_sites)," (",sprintf("%.5f",$S_aDNA_boot[$lower]/$a1_obs_aDNA/$aDNA_sites)," - ",sprintf("%.5f",$S_aDNA_boot[$upper]/$a1_obs_aDNA/$aDNA_sites),")\n";
	print OUTFILE "pi\t",sprintf("%.5f",$asum/($#pi_aDNA+1))," (",sprintf("%.5f",$pi_aDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_aDNA_persite_boot[$upper]),")\n";
	$pi_a = $asum/($#pi_aDNA+1);
	print OUTFILE "d\t",sprintf("%.5f",$aDNA_divergence/$aDNA_sites)," (",sprintf("%.5f",$d_aDNA_boot[$lower]/$aDNA_sites)," - ",sprintf("%.5f",$d_aDNA_boot[$upper]/$aDNA_sites),")\n";
	# apply JC correction
	$JC_divergence_aDNA = (-3/4)*log(1-(4/3)*($aDNA_divergence/$aDNA_sites));
	#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_aDNA),"\n";
	# apply correction for ancestral polymorphism
	$JC_divergence_aDNA = $JC_divergence_aDNA- $asum/($#pi_aDNA+1);
	print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_aDNA)," (",sprintf("%.5f",$JC_AP_divergence_aDNA_boot[$lower])," - ",sprintf("%.5f",$JC_AP_divergence_aDNA_boot[$upper]),")\n";
	print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($asum/($#pi_aDNA+1))/$JC_divergence_aDNA)," (",sprintf("%.5f",$thetapi_over_divergence_aDNA[$lower])," - ",sprintf("%.5f",$thetapi_over_divergence_aDNA[$upper]),")\n";
	# TajD aDNA
	if($n_aDNA > 2){
		$TajD_aDNA = ($asum - ($aDNA_segregating_sites/$a1_obs_aDNA))/((($e1_obs_aDNA*$aDNA_segregating_sites)+($e2_obs_aDNA*$aDNA_segregating_sites*($aDNA_segregating_sites-1)))**0.5);
		print OUTFILE "TajD\t",sprintf("%.3f",$TajD_aDNA)," (",sprintf("%.3f",$aDNA_TajD_boot[$lower])," - ",sprintf("%.3f",$aDNA_TajD_boot[$upper]),")\n";
	}
	else{
	 	$TajD_aDNA = "undefined";
	 	print OUTFILE "TajD\tundefined\n";
	}
	#singleton proportion
	print OUTFILE "Se/S\t",sprintf("%.3f",($aDNA_singleton_sites/$aDNA_segregating_sites))," (",sprintf("%.3f",($aDNA_singleton_sites_boot[$lower]))," - ",sprintf("%.3f",($aDNA_singleton_sites_boot[$upper])),")\n";
}

if($xDNA_sites>0){
	print OUTFILE "\n";
	print OUTFILE "xDNA\n";
	print OUTFILE "#_alleles\t",$n_xDNA,"\n";
	print OUTFILE "#_Sites\t",$xDNA_sites,"\n";
	print OUTFILE "RAD_tag_count\t",$RAD_tag_count_xDNA,"\n";
	print OUTFILE "S\t",$xDNA_segregating_sites," (",$S_xDNA_boot[$lower]," - ",$S_xDNA_boot[$upper],")\n";
	print OUTFILE "thetaW\t",sprintf("%.5f",$xDNA_segregating_sites/$a1_obs_xDNA/$xDNA_sites)," (",sprintf("%.5f",$S_xDNA_boot[$lower]/$a1_obs_xDNA/$xDNA_sites)," - ",sprintf("%.5f",$S_xDNA_boot[$upper]/$a1_obs_xDNA/$xDNA_sites),")\n";
	print OUTFILE "pi\t",sprintf("%.5f",$xsum/($#pi_xDNA+1))," (",sprintf("%.5f",$pi_xDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_xDNA_persite_boot[$upper]),")\n";
	my $pi_x = $xsum/($#pi_xDNA+1);
	print OUTFILE "d\t",sprintf("%.5f",$xDNA_divergence/$xDNA_sites)," (",sprintf("%.5f",$d_xDNA_boot[$lower]/$xDNA_sites)," - ",sprintf("%.5f",$d_xDNA_boot[$upper]/$xDNA_sites),")\n";
	# apply JC correction
	$JC_divergence_xDNA = (-3/4)*log(1-(4/3)*($xDNA_divergence/$xDNA_sites));
	#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_xDNA),"\n";
	# apply correction for ancestral polymorphism
	$JC_divergence_xDNA = $JC_divergence_xDNA- $xsum/($#pi_xDNA+1);
	print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_xDNA)," (",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$lower])," - ",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$upper]),")\n";
	print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)," (",sprintf("%.5f",$thetapi_over_divergence_xDNA[$lower])," - ",sprintf("%.5f",$thetapi_over_divergence_xDNA[$upper]),")\n";
	# TajD xDNA
	if(($n_xDNA > 2)&&(((($e1_obs_xDNA*$xDNA_segregating_sites)+($e2_obs_xDNA*$xDNA_segregating_sites*($xDNA_segregating_sites-1))))>0)){
		$TajD_xDNA = ($xsum - ($xDNA_segregating_sites/$a1_obs_xDNA))/((($e1_obs_xDNA*$xDNA_segregating_sites)+($e2_obs_xDNA*$xDNA_segregating_sites*($xDNA_segregating_sites-1)))**0.5);
		print OUTFILE "TajD\t",sprintf("%.3f",$TajD_xDNA)," (",sprintf("%.3f",$xDNA_TajD_boot[$lower])," - ",sprintf("%.3f",$xDNA_TajD_boot[$upper]),")\n";
	}
	else{
		$TajD_xDNA = "undefined";
	 	print OUTFILE "TajD\tundefined\n";
	}
	#singleton proportion
	print OUTFILE "Se/S\t",sprintf("%.3f",($xDNA_singleton_sites/$xDNA_segregating_sites))," (",sprintf("%.3f",($xDNA_singleton_sites_boot[$lower]))," - ",sprintf("%.3f",($xDNA_singleton_sites_boot[$upper])),")\n";
}

print OUTFILE "\n";
if(0){ # this comments out the yDNA stuff
if(($number_of_individuals_genotyped-$number_of_female_individuals_genotyped)>1){
	print OUTFILE "yDNA\n";
	print OUTFILE "#_alleles\t",$n_yDNA,"\n";
	print OUTFILE "#_Sites\t",$yDNA_sites,"\n";
	print OUTFILE "RAD_tag_count\t",$RAD_tag_count_yDNA,"\n";
	print OUTFILE "S\t",$yDNA_segregating_sites,"\n";
	print OUTFILE "thetaW\t",sprintf("%.5f",$yDNA_segregating_sites/$a1_obs_yDNA/$yDNA_sites),"\n";
	print OUTFILE "pi\t",sprintf("%.5f",$ysum/($#pi_yDNA+1)),"\n";
	my $pi_y = $ysum/($#pi_yDNA+1);
	print OUTFILE "d\t",sprintf("%.5f",$yDNA_divergence/$yDNA_sites),"\n";
	# apply JC correction
	$JC_divergence_yDNA = (-3/4)*log(1-(4/3)*($yDNA_divergence/$yDNA_sites));
	#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_yDNA),"\n";
	# apply correction for ancestral polymorphism
	$JC_divergence_yDNA = $JC_divergence_yDNA- $ysum/($#pi_yDNA+1);
	print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_yDNA),"\n";
	print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($ysum/($#pi_yDNA+1))/$JC_divergence_yDNA),"\n";
	# TajD yDNA
	if(($n_yDNA>2)&&($a1_obs_yDNA>0)&&(($e1_obs_yDNA*$yDNA_segregating_sites)+($e2_obs_yDNA*$yDNA_segregating_sites*($yDNA_segregating_sites-1))>0)){
		$TajD_yDNA = ($ysum - ($yDNA_segregating_sites/$a1_obs_yDNA))/((($e1_obs_yDNA*$yDNA_segregating_sites)+($e2_obs_yDNA*$yDNA_segregating_sites*($yDNA_segregating_sites-1)))**0.5);
		print OUTFILE "TajD\t",sprintf("%.3f",$TajD_yDNA),"\n";
	}
	else{
		$TajD_yDNA = "undefined";
		print OUTFILE "TajD\t",$TajD_yDNA,"\n\n";
	}
	#singleton proportion
	if($yDNA_segregating_sites>0){
		print OUTFILE "Se/S\t",sprintf("%.3f",($yDNA_singleton_sites/$yDNA_segregating_sites)),"\n\n";
	}
}

} # this is the end of commenting for the yDNA stuff

if(($aDNA_sites>0)&&($xDNA_sites>0)){

	print OUTFILE "ratio_of_X_to_A_pi_over_d\n";
	print OUTFILE "pi_X/d_jc_ad_X/pi_a/d_jc_ad_a\t",sprintf("%.3f",(($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA))," ";

	# now calculate and print the confidence intervals for this statistic using the delta method
	# as in Evans et al. 2014
	my $scaled_X_over_A = (($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA);
	my $var_dA_over_dX = (($JC_divergence_aDNA**2)/($JC_divergence_xDNA**4))*(sqrt($JC_divergence_xDNA)/sqrt($xDNA_sites))**2+(1/($JC_divergence_xDNA**2))*(sqrt($JC_divergence_aDNA)/sqrt ($aDNA_sites))**2;
	my $dA_over_dX = ($JC_divergence_aDNA/$JC_divergence_xDNA);
	my $var_piA = ((($n_aDNA+1)*$pi_a)/((3*$n_aDNA-1)*$aDNA_sites))+((2*($n_aDNA**2+$n_aDNA+3)*($pi_a**2))/(9*$n_aDNA*($n_aDNA-1)))/$RAD_tag_count_aDNA;
	my $var_piX = ((($n_xDNA+1)*$pi_x)/((3*$n_xDNA-1)*$xDNA_sites))+((2*($n_xDNA**2+$n_xDNA+3)*($pi_x**2))/(9*$n_xDNA*($n_xDNA-1)))/$RAD_tag_count_xDNA;
	my $var_piX_over_piA = ($pi_x**2/$pi_a**4)*$var_piA+(1/($pi_a**2))*$var_piX;
	my $stdev_scaled_X_over_A = sqrt ($scaled_X_over_A**2*$var_dA_over_dX+$dA_over_dX**2*$var_piX_over_piA);
	print OUTFILE " (",sprintf("%.3f",($scaled_X_over_A-1.96*$stdev_scaled_X_over_A))," - ",sprintf("%.3f",($scaled_X_over_A+1.96*$stdev_scaled_X_over_A)),")\n";

	if(($number_of_individuals_genotyped-$number_of_female_individuals_genotyped)>1){
		print OUTFILE "pi_Y/d_jc_ad_Y/pi_a/d_jc_ad_a is\t",(($ysum/($#pi_yDNA+1))/$JC_divergence_yDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA),"\n";
	}
}

=cut # comments out end part

close DATAINPUT;
close OUTFILE;
```
