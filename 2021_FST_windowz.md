Perl script for making input files for plotting
```perl
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
# ./All_N_mt_allinteract_fst_permutation.pl FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez.txt fst_bor/nem_mau_windowstats.concat

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

# this will print out a file that I may use for plotting later

my $outputfile = $inputfile2."_FST__density.txt"; # the name of the output file is from the commandline
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}

#print OUTFILE "FST\tN_interact\n";

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
					$gene_containing_window=1;
					$number_of_genes_in_this_window+=1; # only count the start of genes in windows
					$OXPHOS{$key}{"start_fst"} = $temp[8];
					$OXPHOS{$key}{"start_fst_sites"} = $temp[4];
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
			if(($temp1[0] eq $temp[0])&&($temp1[2] >= $temp[1])&&($temp1[2] <= $temp[2])) {  # end is in block
					$gene_containing_window=1;
					#$number_of_genes_in_this_window+=1; # if we count the end of genes too, we end up double counting
					$OXPHOS{$key}{"end_fst"} = $temp[8];
					$OXPHOS{$key}{"end_fst_sites"} = $temp[4];
					if($OXPHOS{$key}{"mt_interact"} == 1){
						$N_interact_window=1;
					#	$number_of_Ninteract_genes_in_this_window+=1; # if we count the end of genes too, we end up double counting
					}	
			}
		}
		if($temp[0] ne 'scaffold'){
			print OUTFILE $temp[0],"\t",$temp[1],"\t",$temp[8],"\t",$gene_containing_window,"\t",$N_interact_window,
			"\t",$number_of_genes_in_this_window,"\t",$number_of_Ninteract_genes_in_this_window,"\t",$Ninteract_acronym,"\n";
		}	
}

close OUTFILE;

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
	#else{
	#	print "Problem with fst value for ",$key,"\n";
	#	print "The problem is that there is no fst value for the gene\n";
	#}							
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

my $N_interact_counter=0;

# now calculate the average fst for associated and non-associated OXPHOS genes
foreach my $key (keys %OXPHOS){
	if(exists($OXPHOS{$key}{"weighted_ave_fst"})){
		if($OXPHOS{$key}{"mt_interact"} == 1){
			$N_interact_counter+=1;
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
print "Number of associated blocks ",$n_Fst_associated,"\n";
print "Number of N_interact genes ",$N_interact_counter,"\n";
print "Mean Fst all associated N_mt genes ",$Fst_associated/$n_Fst_associated,"\n";
print "Mean Fst complex 1 associated ",$Fst_associated_complex1/$n_Fst_associated_complex1,"\n";
print "Mean Fst complex 3 associated ",$Fst_associated_complex3/$n_Fst_associated_complex3,"\n";
print "Mean Fst complex 4 associated ",$Fst_associated_complex4/$n_Fst_associated_complex4,"\n";
print "Mean Fst complex 5 associated ",$Fst_associated_complex5/$n_Fst_associated_complex5,"\n";
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

my $perms=1;
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

R script for plotting:
```R
library(ggrepel)
library(ggplot2)
setwd('/Users/Shared/Previously Relocated Items/Security/projects/2017_SEAsian_macaque_genomz/Fst_windows/all_fst_outputz')

pairwise_vector <- c("nem_mau", "nem_ton", "nem_hec", "nem_nig", "mau_ton",
                     "mau_hec", "mau_nig", "ton_hec", "ton_nig", "hec_nig")
pairwise_vector_plot <- c("nem_mau_plot", "nem_ton_plot", "nem_hec_plot", "nem_nig_plot", "mau_ton_plot",
                     "mau_hec_plot", "mau_nig_plot", "ton_hec_plot", "ton_nig_plot", "hec_nig_plot")

for (pair in pairwise_vector){
  a <- read.table(paste(eval(pair), "_windowstats.concat_FST__density.txt", sep=""), header = T)
  assign(pair,a)
}

# now add faceting variables to each dataset 
nem_mau$pair <- "nem_mau"
nem_ton$pair <- "nem_ton"
nem_hec$pair <- "nem_hec"
nem_nig$pair <- "nem_nig"
mau_ton$pair <- "mau_ton"
mau_hec$pair <- "mau_hec"
mau_nig$pair <- "mau_nig"
ton_hec$pair <- "ton_hec"
ton_nig$pair <- "ton_nig"
hec_nig$pair <- "hec_nig"



#nem_mau
    my_data_only_genez <- nem_mau[nem_mau$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                      (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                              (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    nem_mau_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
                       ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
   
#nem_ton
    my_data_only_genez <- nem_ton[nem_ton$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    nem_ton_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
#nem_hec
    my_data_only_genez <- nem_hec[nem_hec$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    nem_hec_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
    #nem_nig
    my_data_only_genez <- nem_nig[nem_nig$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    nem_nig_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
    #mau_ton
    my_data_only_genez <- mau_ton[mau_ton$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    mau_ton_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
    #mau_hec
    my_data_only_genez <- mau_hec[mau_hec$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    mau_hec_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
    #mau_nig
    my_data_only_genez <- mau_nig[mau_nig$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    mau_nig_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")

    #ton_hec
    my_data_only_genez <- ton_hec[ton_hec$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    ton_hec_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")

    #ton_nig
    my_data_only_genez <- ton_nig[ton_nig$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    ton_nig_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
    #hec_nig
    my_data_only_genez <- hec_nig[hec_nig$containsgenes == 1,] 
    
    # explore relationship between Fst and number of genes in Ninteract windows
    my_data_only_nonNinteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 0,] 
    my_data_only_Ninteractgenez <- my_data_only_genez[my_data_only_genez$containsNinteractgenez == 1,] 
    my_data_only_nonNinteractgenez<- my_data_only_nonNinteractgenez[complete.cases(my_data_only_nonNinteractgenez),]
    
    my_data_only_genez<- my_data_only_genez[complete.cases(my_data_only_genez),]
    #dim(my_data_only_genez)
    #head(my_data_only_genez)
    
    # calculate a lm for all data (because there was not a significant interaction term)
    mod <- lm(Fst ~ number_of_genes, data=my_data_only_genez)
    # get the fitted values (y = mx+b)
    fitted <- mod$coefficients[2]*my_data_only_genez$number_of_genes + mod$coefficients[1]
    #cbind fitted to data
    my_data_only_genez <- cbind(my_data_only_genez,fitted)
    
    # calculate cooks d for all data
    cooksd <- cooks.distance(mod)
    #cbind cooksd to data
    my_data_only_genez <- cbind(my_data_only_genez,cooksd)
    # make a column to specify whether a gene is an Ninteract gene or not
    my_data_only_genez$color <- ifelse(my_data_only_genez$containsNinteractgenez == 1, "pink", "gray")
    my_data_only_genez$alpha <- ifelse(my_data_only_genez$color == "gray", 0.7, 1)
    # make a column that specifies whether cooksd suggests an outlier
    # but only for Ninteract genes
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "red"
    # now color the ones that are below the fitted line blue
    my_data_only_genez$color[(cooksd >=4*mean(cooksd, na.rm=T)) &
                               (my_data_only_genez$Fst < my_data_only_genez$fitted) &
                               (my_data_only_genez$containsNinteractgenez == 1)] <-  "blue"
    
    # make the color column into an ordered factor
    #my_data_only_genez$color <- factor(my_data_only_genez$color, levels = c("gray", "pink", "red"), 
    #                                   ordered = is.ordered(my_data_only_genez$color))
    # on now plot the data with the color representing outliers for N_interact genes
    hec_nig_plot <- ggplot(my_data_only_genez, aes(x = number_of_genes, y = Fst, color=color, fill=color)) +
      geom_smooth(data = my_data_only_genez, method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color == 'gray'), method=lm, se=T, fullrange=TRUE, colour="gray", fill = "gray") +
      #geom_smooth(data = subset(my_data_only_genez, color != 'gray'), method=lm, se=T, fullrange=TRUE, colour="pink", fill = "pink") +
      geom_point(data = subset(my_data_only_genez, color == "gray"),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'gray') +
      geom_point(data = subset(my_data_only_genez, color == 'pink'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha), color = 'pink') +
      geom_point(data = subset(my_data_only_genez, color == 'red'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_point(data = subset(my_data_only_genez, color == 'blue'),
                 aes(x = number_of_genes, y = Fst, alpha = alpha),color = 'red') +
      geom_text_repel( data = my_data_only_genez,
                       #mapping = aes(label = Ninteract_acronym),
                       mapping = aes(label = ifelse(color == "red",as.character(Ninteract_acronym),'')),
                       #force_pull = 0,
                       force = 1.2,
                       nudge_x = 15,
                       color = "black",
                       size = 2,
                       box.padding = 0.35, 
                       #point.padding = 0.5,
                       direction     = "y",
                       max.overlaps = Inf,
                       hjust = 0,
                       segment.size = 0.25,
                       segment.color = 'grey50'
      ) +
      labs(x = "Number of genes in window", y=expression(paste(italic(F[ST])))) +
      theme_bw() + theme(legend.position = "none")
    
    library(gridExtra)
    
    grid.arrange(nem_mau_plot, nem_ton_plot,
                 nem_hec_plot, nem_nig_plot,
                 mau_ton_plot, mau_hec_plot,
                 mau_nig_plot, ton_hec_plot,
                 ton_nig_plot, hec_nig_plot,ncol=2)
    
 ```
