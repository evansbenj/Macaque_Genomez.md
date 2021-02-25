I used this script to do the permutations for N_mt interact genes vs all genes, or vs only other non-direct interactors from OXPHOS, MRP, and ARG2:
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
# ./All_N_mt_interact_fst_permutation.pl OXPHOS_ARP2_MRP_andallother_genez.txt *concat

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

# first open up the OXPHOS gene info (# OXPHOS_ARP2_MRP_andallother_genez.txt)
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(($temp[0] ne 'gene')&&($temp[2] ne 'chrX')){ # deliberately blocks with no fst value
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0];
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"complex"} = $temp[1];
		$OXPHOS{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5];
	}	
}		
close DATAINPUT;

# now open up the Fst data
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

my @temp1;
while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split(',',$line);
		# cycle through each gene
		foreach my $key (keys %OXPHOS){
			@temp1=split('_',$key);
			# check if this window contains one or more N_mt genes
			if(($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= $temp[2])){
					$OXPHOS{$key}{"start_fst"} = $temp[8];
					$OXPHOS{$key}{"start_fst_sites"} = $temp[4];
			} # start is in block
			if(($temp1[0] eq $temp[0])&&($temp1[2] >= $temp[1])&&($temp1[2] <= $temp[2])) {  # end is in block
					$OXPHOS{$key}{"end_fst"} = $temp[8];
					$OXPHOS{$key}{"end_fst_sites"} = $temp[4];
			}
		}
}

# Now the OXPHOS hash has Fst and coordinates of all blocks that have genes

my @weighted_ave_fst_for_perms; # this has only the mtinteractors and all genes
my @weighted_ave_fst_for_perms_all_N_mt; # this has only the mtinteractors and the non-interactors that are still OXPHOS, MRP, or ARP2
my @weighted_ave_fst_for_perms_complex1; # this has only the mtinteracters for complex1 and the non-interactors for all complexes
my @weighted_ave_fst_for_perms_complex3; # this has only the mtinteracters for complex3 and the non-interactors for all complexes
my @weighted_ave_fst_for_perms_complex4; # this has only the mtinteracters for complex4 and the non-interactors for all complexes
my @weighted_ave_fst_for_perms_complex5; # this has only the mtinteracters for complex5 and the non-interactors for all complexes

# now calculate the weighted averages for each gene
# this is necessary because the block that has the start site may be different
# from the one that has the end site, and these may have different numbers of 
# sites.  A weighted average will upweight the one with more sites
foreach my $key (keys %OXPHOS){
	if((exists($OXPHOS{$key}{"end_fst"}))&&(exists($OXPHOS{$key}{"start_fst"}))
		&&($OXPHOS{$key}{"start_fst"} ne 'nan')
		&&($OXPHOS{$key}{"end_fst"} ne 'nan')
		){
		# build @weighted_ave_fst_for_perms that includes all genes
		$OXPHOS{$key}{"weighted_ave_fst"}=(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"});
		push(@weighted_ave_fst_for_perms,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		# now build @weighted_ave_fst_for_perms_all_N_mt 
		# this includes only N_mt genes (including direct and non-direct interactions: OXPHOS, ARP2, MRP)
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} ne 0)){ # this is a non-interacting N_mt gene
			push(@weighted_ave_fst_for_perms_all_N_mt,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} ne 0)){ 
			push(@weighted_ave_fst_for_perms_all_N_mt,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		# build @weighted_ave_fst_for_perms_complex1
		# this includes only proteins in complex 1
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '1')){ 
			push(@weighted_ave_fst_for_perms_complex1,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '1')){ 
			push(@weighted_ave_fst_for_perms_complex1,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		# build @weighted_ave_fst_for_perms_complex3
		# this includes only proteins in complex 3
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '3')){ 
			push(@weighted_ave_fst_for_perms_complex3,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '3')){ 
			push(@weighted_ave_fst_for_perms_complex3,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		# build @weighted_ave_fst_for_perms_complex4
		# this includes only proteins in complex 4
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '4')){ 
			push(@weighted_ave_fst_for_perms_complex4,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '4')){ 
			push(@weighted_ave_fst_for_perms_complex4,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		# build @weighted_ave_fst_for_perms_complex5
		# this includes only proteins in complex 5
		if(($OXPHOS{$key}{"mt_interact"} == 0)&&($OXPHOS{$key}{"complex"} eq '5')){ 
			push(@weighted_ave_fst_for_perms_complex5,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
		elsif(($OXPHOS{$key}{"mt_interact"} == 1)&&($OXPHOS{$key}{"complex"} eq '5')){ 
			push(@weighted_ave_fst_for_perms_complex5,(($OXPHOS{$key}{"start_fst"}*$OXPHOS{$key}{"start_fst_sites"})+
								($OXPHOS{$key}{"end_fst"}*$OXPHOS{$key}{"end_fst_sites"}))/
								($OXPHOS{$key}{"start_fst_sites"}+$OXPHOS{$key}{"end_fst_sites"}));
		}
	}
	else{
		print "Problem with fst value for ",$key,"\n";
		print "The problem is that there is no fst value for the gene\n";
	}							
}	

close DATAINPUT2;

my $Fst_associated=0; # all OXPHOS,MRP, ARP2 genes
my $Fst_associated_complex1=0;
my $Fst_associated_complex3=0;
my $Fst_associated_complex4=0;
my $Fst_associated_complex5=0;

my $Fst_non_associated=0;
my $Fst_non_associated_only_N_mt=0;
my $Fst_non_associated_complex1=0;
my $Fst_non_associated_complex3=0;
my $Fst_non_associated_complex4=0;
my $Fst_non_associated_complex5=0;


my $n_Fst_associated=0; # this also works for only_N_mt
my $n_Fst_associated_complex1=0;
my $n_Fst_associated_complex3=0;
my $n_Fst_associated_complex4=0;
my $n_Fst_associated_complex5=0;

my $n_Fst_non_associated=0; # this includes all non associated genes, including those anywhere 
my $n_Fst_non_associated_only_N_mt=0; # this includes only non associated N-mt genes 
									  # (i.e., OXPHOS, MRP, ARP2 that don't directly interact with mt genes)
my $n_Fst_non_associated_complex1=0;
my $n_Fst_non_associated_complex3=0;
my $n_Fst_non_associated_complex4=0;
my $n_Fst_non_associated_complex5=0;

# now calculate the average fst for associated and non-associated OXPHOS genes
foreach my $key (keys %OXPHOS){
	if(exists($OXPHOS{$key}{"weighted_ave_fst"})){
		if($OXPHOS{$key}{"mt_interact"} == 1){
			$Fst_associated += $OXPHOS{$key}{"weighted_ave_fst"};
			$n_Fst_associated += 1;	
			if($OXPHOS{$key}{"complex"} eq 1){
				$Fst_associated_complex1+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_associated_complex1 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 3){
				$Fst_associated_complex3+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_associated_complex3 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 4){
				$Fst_associated_complex4+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_associated_complex4 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 5){
				$Fst_associated_complex5+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_associated_complex5 +=1;
			}	
		}
		elsif($OXPHOS{$key}{"mt_interact"} == 0){
			$Fst_non_associated += $OXPHOS{$key}{"weighted_ave_fst"};
			$n_Fst_non_associated += 1;
			if($OXPHOS{$key}{"complex"} ne 0){
				$Fst_non_associated_only_N_mt += $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_non_associated_only_N_mt += 1;
			}
			if($OXPHOS{$key}{"complex"} eq 1){
				$Fst_non_associated_complex1+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_non_associated_complex1 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 3){
				$Fst_non_associated_complex3+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_non_associated_complex3 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 4){
				$Fst_non_associated_complex4+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_non_associated_complex4 +=1;
			}	
			if($OXPHOS{$key}{"complex"} eq 5){
				$Fst_non_associated_complex5+= $OXPHOS{$key}{"weighted_ave_fst"};
				$n_Fst_non_associated_complex5 +=1;
			}	
		}
	}
}


# now report values
print "Mean Fst all associated N_mt genes ",$Fst_associated/$n_Fst_associated,"\n";
print "Mean Fst complex 1 associated ",$Fst_associated_complex1/$n_Fst_associated_complex1,"\n";
print "Mean Fst complex 3 associated ",$Fst_associated_complex3/$n_Fst_associated_complex3,"\n";
print "Mean Fst complex 4 associated ",$Fst_associated_complex4/$n_Fst_associated_complex4,"\n";
print "Mean Fst complex 5 associated ",$Fst_associated_complex5/$n_Fst_associated_complex5,"\n";
print "Mean Fst non-associated all genes ",$Fst_non_associated/$n_Fst_non_associated,"\n";
print "Mean Fst non-associated only N_mt ",$Fst_non_associated_only_N_mt/$n_Fst_non_associated_only_N_mt,"\n";


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
if($#associated_or_not_array ne $#weighted_ave_fst_for_perms){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#weighted_ave_fst_for_perms,"\n";
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
			$Fst_associated_perm+= $weighted_ave_fst_for_perms[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $weighted_ave_fst_for_perms[$x];
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



#################
# ALL COMPLEXES perms including only N_mt genes
#################

# calculate test statistic for N_mt genes
# the first part ($Fst_associated/$n_Fst_associated) is the same as for all genes
$test_stat = ($Fst_associated/$n_Fst_associated) - ($Fst_non_associated_only_N_mt/$n_Fst_non_associated_only_N_mt);

# do permutations for all_complex_mt_interact vs all_complex_no_mt_interact
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
@associated_or_not_array = (('1') x $n_Fst_associated, ('0') x $n_Fst_non_associated_only_N_mt);

$perms=1000;
$Fst_associated_perm=0;
$Fst_not_associated_perm=0;
$n_Fst_associated_perm=0;
$n_Fst_not_associated_perm=0;


# first check if the length of the permutation array is the same as the fst array
# for this analysis
if($#associated_or_not_array ne $#weighted_ave_fst_for_perms_all_N_mt){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#weighted_ave_fst_for_perms_all_N_mt,"\n";
}

@perm_diffs={};
for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$Fst_associated_perm=0;
	$Fst_not_associated_perm=0;
	$n_Fst_associated_perm=0;
	$n_Fst_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$Fst_associated_perm+= $weighted_ave_fst_for_perms_all_N_mt[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $weighted_ave_fst_for_perms_all_N_mt[$x];
			$n_Fst_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($Fst_associated_perm/$n_Fst_associated_perm) - ($Fst_not_associated_perm/$n_Fst_not_associated_perm));	
}

@perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
$switch=0;
$pval=0;
$counter=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for test including only N_mt genes (if negative, then not significant):",$test_stat,"\n";
print "P = ",1-($pval/$perms),"\n";



#################
# COMPLEX 1 perms
#################


# calculate test statistic
$test_stat = ($Fst_associated_complex1/$n_Fst_associated_complex1) - ($Fst_non_associated_complex1/$n_Fst_non_associated_complex1);

# do permutations for complex1_mt_interact vs all_complexes_no_mt_interact
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
@associated_or_not_array={};
@associated_or_not_array = (('1') x $n_Fst_associated_complex1, ('0') x $n_Fst_non_associated_complex1);
#print "@associated_or_not_array\n";

$Fst_associated_perm=0;
$Fst_not_associated_perm=0;
$n_Fst_associated_perm=0;
$n_Fst_not_associated_perm=0;

# first check if the length of the permitation array is the same as the fst array
# for this analysis
if($#associated_or_not_array ne $#weighted_ave_fst_for_perms_complex1){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#weighted_ave_fst_for_perms_complex1,"\n";
}


@perm_diffs={};
for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$Fst_associated_perm=0;
	$Fst_not_associated_perm=0;
	$n_Fst_associated_perm=0;
	$n_Fst_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$Fst_associated_perm+= $weighted_ave_fst_for_perms_complex1[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $weighted_ave_fst_for_perms_complex1[$x];
			$n_Fst_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($Fst_associated_perm/$n_Fst_associated_perm) - ($Fst_not_associated_perm/$n_Fst_not_associated_perm));	
}

@perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
$switch=0;
$pval=0;
$counter=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for complex1 (if negative, then not significant):",$test_stat,"\n";
print "P_complex1 = ",1-($pval/$perms),"\n";



#################
# COMPLEX 3 perms
#################


# calculate test statistic
$test_stat = ($Fst_associated_complex3/$n_Fst_associated_complex3) - ($Fst_non_associated_complex3/$n_Fst_non_associated_complex3);

# do permutations for complex1_mt_interact vs all_complexes_no_mt_interact
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
@associated_or_not_array={};
@associated_or_not_array = (('1') x $n_Fst_associated_complex3, ('0') x $n_Fst_non_associated_complex3);
#print "@associated_or_not_array\n";

$Fst_associated_perm=0;
$Fst_not_associated_perm=0;
$n_Fst_associated_perm=0;
$n_Fst_not_associated_perm=0;

# first check if the length of the permitation array is the same as the fst array
# for this analysis
if($#associated_or_not_array ne $#weighted_ave_fst_for_perms_complex3){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#weighted_ave_fst_for_perms_complex3,"\n";
}


@perm_diffs={};
for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$Fst_associated_perm=0;
	$Fst_not_associated_perm=0;
	$n_Fst_associated_perm=0;
	$n_Fst_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$Fst_associated_perm+= $weighted_ave_fst_for_perms_complex3[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $weighted_ave_fst_for_perms_complex3[$x];
			$n_Fst_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($Fst_associated_perm/$n_Fst_associated_perm) - ($Fst_not_associated_perm/$n_Fst_not_associated_perm));	
}

@perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
$switch=0;
$pval=0;
$counter=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for complex3 (if negative, then not significant):",$test_stat,"\n";
print "P_complex3 = ",1-($pval/$perms),"\n";



#################
# COMPLEX 4 perms
#################


# calculate test statistic
$test_stat = ($Fst_associated_complex4/$n_Fst_associated_complex4) - ($Fst_non_associated_complex4/$n_Fst_non_associated_complex4);

# do permutations for complex1_mt_interact vs all_complexes_no_mt_interact
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
@associated_or_not_array={};
@associated_or_not_array = (('1') x $n_Fst_associated_complex4, ('0') x $n_Fst_non_associated_complex4);
#print "@associated_or_not_array\n";

$Fst_associated_perm=0;
$Fst_not_associated_perm=0;
$n_Fst_associated_perm=0;
$n_Fst_not_associated_perm=0;


# first check if the length of the permitation array is the same as the fst array
# for this analysis
if($#associated_or_not_array ne $#weighted_ave_fst_for_perms_complex4){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#weighted_ave_fst_for_perms_complex4,"\n";
}


@perm_diffs={};

for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$Fst_associated_perm=0;
	$Fst_not_associated_perm=0;
	$n_Fst_associated_perm=0;
	$n_Fst_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$Fst_associated_perm+= $weighted_ave_fst_for_perms_complex4[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $weighted_ave_fst_for_perms_complex4[$x];
			$n_Fst_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($Fst_associated_perm/$n_Fst_associated_perm) - ($Fst_not_associated_perm/$n_Fst_not_associated_perm));	
}

@perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
$switch=0;
$pval=0;
$counter=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for complex4 (if negative, then not significant):",$test_stat,"\n";
print "P_complex4 = ",1-($pval/$perms),"\n";


#################
# COMPLEX 5 perms
#################


# calculate test statistic
$test_stat = ($Fst_associated_complex5/$n_Fst_associated_complex5) - ($Fst_non_associated_complex5/$n_Fst_non_associated_complex5);

# do permutations for complex1_mt_interact vs all_complexes_no_mt_interact
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
@associated_or_not_array={};
@associated_or_not_array = (('1') x $n_Fst_associated_complex5, ('0') x $n_Fst_non_associated_complex5);
#print "@associated_or_not_array\n";

$Fst_associated_perm=0;
$Fst_not_associated_perm=0;
$n_Fst_associated_perm=0;
$n_Fst_not_associated_perm=0;

# first check if the length of the permitation array is the same as the fst array
# for this analysis
if($#associated_or_not_array ne $#weighted_ave_fst_for_perms_complex5){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length fst array ",$#weighted_ave_fst_for_perms_complex5,"\n";
}


@perm_diffs={};

for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$Fst_associated_perm=0;
	$Fst_not_associated_perm=0;
	$n_Fst_associated_perm=0;
	$n_Fst_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$Fst_associated_perm+= $weighted_ave_fst_for_perms_complex5[$x];
			$n_Fst_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$Fst_not_associated_perm+= $weighted_ave_fst_for_perms_complex5[$x];
			$n_Fst_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($Fst_associated_perm/$n_Fst_associated_perm) - ($Fst_not_associated_perm/$n_Fst_not_associated_perm));	
}

@perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
$switch=0;
$pval=0;
$counter=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for complex5 (if negative, then not significant):",$test_stat,"\n";
print "P_complex5 = ",1-($pval/$perms),"\n";


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
