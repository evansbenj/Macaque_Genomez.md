# Here are the perms for the introgression of N_interact loci in PF511 and PF626:
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# OXPHOS and non OXPHOS genes and their acronyms, and all other genes

# The other file is a file admix introgression blocks concatenated for all chrs.  

# First identify how many interacting OXPHOS genes are within introgression blocks
# then scramble them and check how many are expected by chance.

# run like this:
# ./Introgression_block_permutation.pl OXPHOS_ARP2_MRP_allinteractexceptC2_andallother_genez.txt ./TON_HEC_MAU/tonk_PF511_for_perm_.concat

my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];

my @windowsites;
my @Fst_values;
my $sumsites=0;
my @temp;
my $y;
my $x;
my %N_interact_hash;
my @interact_perm;

# first open up the N_interact_hash gene info (OXPHOS_ARP2_MRP_allinteractexceptC2_andallother_genez.txt)
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'gene'){ 
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0];
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5];
		push(@interact_perm,$temp[5]); # this will be used for the permutations later
	}	
}		
close DATAINPUT;

# now open up the introgression data
# consider an introgression window as
# any window with the homoz TON probability <0.5
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

my @temp1;
my $n_introgression_blocks_with_interacting_genez=0;
my $n_introgression_blocks_with_other_genez=0;
my $n_introgression_blocks_without_genez=0;
my $admixfrog_block_size=1000000; # yes, the admixfrog blocks are 1million bp!
my %introgression_blocks;
my %perm_blockz; 

while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split(',',$line);
	# ignore first line
	if($temp[0] ne 'chrom'){
		# check if this is an introgression block

		if($temp[7] < 0.5){ 
			$perm_blockz{$temp[0]."_".$temp[2]}{"introgression"}=1; # this is an introgression block
			# this is an introgression block; 
			# block size is $admixfrog_block_size bp
			# initially assign the block to have no genes 
			$introgression_blocks{$temp[0]."_".$temp[2]}{"genes"} = 0; # key has the lower limit of window
																	   # from admix frog.  the upper limit
																	   # is this plus $admixfrog_block_size
																	   # minus 1
			# also assume that it does not have any interacting genes
			$introgression_blocks{$temp[0]."_".$temp[2]}{"interacting"} = 0;

			# Now cycle through all the genes to see if any are in this block
			foreach my $key (keys %N_interact_hash){
				@temp1=split('_',$key);
				# now check if this block contains any genes
				if(
					($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[2])&&($temp1[1] <= ($temp[2]+$admixfrog_block_size-1))
					# beginning is in this block
					|
					($temp1[0] eq $temp[0])&&($temp1[2] >= $temp[2])&&($temp1[2] <= ($temp[2]+$admixfrog_block_size-1))
					# end is in this block	
					){
						# this block has a gene
						$introgression_blocks{$temp[0]."_".$temp[2]}{"genes"} = 1;
						# check if it is an interacting gene
						if($N_interact_hash{$key}{"mt_interact"} == 1){
							$introgression_blocks{$temp[0]."_".$temp[2]}{"interacting"} = 1;
						}
						# no else because we don't want to erase that assignment if we had an interacting
						# gene in this block already
				}
			}
		}
		else{
			$perm_blockz{$temp[0]."_".$temp[2]}{"introgression"}=0; # this is not an introgression block
		}
	}	
}

# ok now I have a hash that has information on whether or not a block has any genes
# and whether or not any of these genes have any N_interact genes.
# print these numbers

foreach my $key (keys %introgression_blocks){
	if($introgression_blocks{$key}{"genes"} == 1){
		if($introgression_blocks{$key}{"interacting"} == 1){
			$n_introgression_blocks_with_interacting_genez+=1;
			print "Introgression_with_interacting ",$key,"\n";
		}
		else{
			$n_introgression_blocks_with_other_genez+=1;
		}	
	}
	else{
		$n_introgression_blocks_without_genez+=1;
	}	
}	

print "Number of introgression blocks with N_interact genes: ",$n_introgression_blocks_with_interacting_genez,"\n";
print "Number of introgression blocks with other genes: ",$n_introgression_blocks_with_other_genez,"\n";
print "Number of introgression blocks without genes: ",$n_introgression_blocks_without_genez,"\n";
print "Proportion of blocks with genes that have N_interact genes",$n_introgression_blocks_with_interacting_genez/
($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez),"\n";

my $size = keys %N_interact_hash;
print "hello ",$#interact_perm,"\t", $size,"\n";

######################
# Permutations
######################

my $perms=1000;
my $counter=0;
my $introg_interacter=0;
my $introg_withgenez=0;
my $introg_without_genez=0;
my @permed_interactorz;
my $introg_withgenez_switch=0;
my $introg_interacter_switch=0;

for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@interact_perm );    # permutes the N_interact assignment for each gene
	# reset the variables used in the tabulation for each perm
	$introg_withgenez=0;
	$introg_interacter=0;
	$introg_without_genez=0;
	foreach my $permkey (keys %perm_blockz){ # go through all the introgression blocks
		if($perm_blockz{$permkey}{"introgression"} == 1){ # this is an introgression block
			# reset the variables used in the tabulation for each block
			$introg_withgenez_switch=0;
			$introg_interacter_switch=0;
			$counter=0;
			# the key is chr_start for blocks of $admixfrog_block_size base pairs
			@temp=split('_',$permkey);
			# for each introgression block go through all genes to see if ther is an interaction gene
			foreach my $key (keys %N_interact_hash){ # go through each gene where the key is chr_start_end
				@temp1=split('_',$key);
				# now check if this block contains any genes
					if(
						($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= ($temp[1]+$admixfrog_block_size-1))
						# beginning is in this block
						|
						($temp1[0] eq $temp[0])&&($temp1[2] >= $temp[1])&&($temp1[2] <= ($temp[1]+$admixfrog_block_size-1))
						# end is in this block	
						){
							# this block has a gene
							$introg_withgenez_switch = 1;
							#print "counter ",$counter,"\n";
							# check if it is a (permuted) interacting gene
							if($interact_perm[$counter] == 1){
								$introg_interacter_switch = 1;
							}
							# no else because we don't want to erase that assignment if we had an interacting
							# gene in this block already
					}
					$counter+=1; # increment for every gene
			} # ok for this block we checked all genes to see if there are any in this block 
			  # and if so if any are interactors
			if($introg_withgenez_switch == 1){
				$introg_withgenez+=1; # this keeps track of the number of introgression blocks with genes
									# for each perm - should be the same as the observed
				if($introg_interacter_switch == 1){
					$introg_interacter += 1; # this is the one that should differ among the perms and the observed
											 # this is the number of interactors by chance
				}
			}
			else{
				# no gene here!
				$introg_without_genez+=1; # this should be the same as the observed
			}
		}
	}	
	# ok now for that perm, push the proportion of introgression blocks with interactors into an array
	push(@permed_interactorz,($introg_interacter/$introg_withgenez));
	print $introg_interacter,"\t",$introg_withgenez,"\t",$introg_without_genez,"\n";
}

my $test_stat = $n_introgression_blocks_with_interacting_genez/
($n_introgression_blocks_with_interacting_genez+$n_introgression_blocks_with_other_genez);
my @permz_sorted = sort { $a <=> $b } @permed_interactorz;
my $switch=0;
my $pval=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#permz_sorted; $y++ ) {
	if(($test_stat <= $permz_sorted[$y])&&($switch==0)){
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
