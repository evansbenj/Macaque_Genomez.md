# MtDNA and yDNA

On redfin, goblin, etc, I am working in this directory: `/work/ben/2017_SEAsian_macaques/ben_scripts`

Because these loci are haploid, we need to treat them differently from autosomal DNA (and maybe xDNA). One option is to just call bases based on the highest depth of coverage.  I have a script to do this that starts from a combined vcf file.  First I generate a combined vcf file like this, on redfin (9_gatk_GenotypeGVCFs_bychr_sqsub_females_and_males.pl):

``` perl
#!/usr/bin/perl
# This script will run quake on trimmed fq files

my $gatkpath = "/work/ben/2017_SEAsian_macaques/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/";
#my $referencegenome="/scratch/ben/MacaM/MacaM_mt_y.fa";
my $referencegenome="/work/ben/2017_SEAsian_macaques/MacaM/MacaM_mt_female.fa";
my $majorpathfemales = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females/";
my $majorpathmales = "/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/males/";
my @chromosomes =("chr01","chr02a","chr02b","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","ch
r19","chrX","chrM","chrY");

    foreach my $chromosome (@chromosomes){
	if($chromosome =~ /chrM/){
	    my $commandline = "sqsub -r 2d --mpp 16G -o catvar_".$chromosome."\.log ";
	    $commandline = $commandline."/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx8G -jar ".$gatkpath."GenomeAnalysisTK.jar -T GenotypeGVCFs -R ".$referencegenome;
		# the combined g.vcf files are in the female directory
		@files = glob($majorpathfemales."all_".$chromosome."_noBSQR.g.vcf.gz");
		foreach $file (@files){
		    $commandline = $commandline." -V ".$file;
		}
		$commandline = $commandline." --includeNonVariantSites -o ".$majorpathfemales."all_".$chromosome."_noBSQR_allsites.vcf.gz";
	    print $commandline,"\n";
#    $status = system($commandline);
	 }
 }

```

# Genotype by depth

I genotyped mtDNA by depth by setting all individuals as zero (males) using the script below.

10_Genotypes_only_male_chrX_based_on_allelic_depth.pl:

``` perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::Util 'max';
use List::Util qw(shuffle);


# This program reads in a vcf file then genotypes chrX sequences
# based on the AD (allelic depth) annotation for males
# Females are left as is

# To run type this:
# Genotypes_only_male_chrX_based_on_allelic_depth copy input.vcf 0101010 output.tab

# where 0101010 indicates for each ingroup 
# sample whether the individual is not (0) or is (1)
# a female

# It takes as input a vcf file and outputs a tab delimited file


my $inputfile = $ARGV[0];
my $outputfile = $ARGV[2];

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

if($inputfile =~ /.gz/){
    print "opening compressed file\n";
#    unless (open DATAINPUT, '/work/ben/2017_SEAsian_macaques/bin/htslib-1.6/bin/bgzip -c $inputfile | ') {
#    open DATAINPUT, '<:gzip' $inputfile or die "Can not find the input file, jackass.\n";
    open DATAINPUT, '-|','gunzip','-c',$inputfile;
}
else{
    unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file, jackass.\n";
	exit;
    }
}

my @sexes = split("",$ARGV[1]);

my $y;
my $x;
my @columns=();
my @fields;
my $AD;
my $GT;
my $counter=0;
my @genotypes;
my $genotypez;
my @alleledepth;
my $max;
my @maxcounter=();
my $counter2=0;
my @altalleles=();
my @allelieos;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@columns=split("\t",$line);
		if(substr($columns[0],0,1) ne '#'){ # this is not a comment
		    @fields=split(":",$columns[8]);
			$counter=0;
			$AD=0;
			$GT=0;
			@altalleles=split(",",$columns[4]);
			# first find out where the AD and GT columns are
			foreach(@fields){
				if($_ eq 'AD'){
					$AD=$counter;
				}
				elsif($_ eq 'GT'){
					$GT=$counter;
				}
				$counter+=1;
			}
			# now print out genotypes
			# first check if we have no data for all individuals
			$genotypez=();
			for ($y = 9 ; $y <= $#columns; $y++ ) {
				@genotypes=split(":",$columns[$y]);
				$genotypez=$genotypez.$genotypes[$GT];
			}
			#print "genotypez ",$genotypez,"\n";
			if(
				(index($genotypez,'0') != -1)||
				(index($genotypez,'1') != -1)||
				(index($genotypez,'2') != -1)){
				# there is at least one genotype in the ingroup
				# if $AD==0 then all individuals are ref
				if($AD==0){ # this probably never happens
					print OUTFILE $columns[0],"\t",$columns[1],"\t",$columns[3];
					for ($y = 9 ; $y <= $#columns; $y++ ) {
						@genotypes=split(":",$columns[$y]);
						if($genotypes[$GT] eq '.\/.'){
							# check if this is a male or female
							if($sexes[$y-9] == 0){ # it is a male
								print OUTFILE "\t\.\/";
							}
							else{ # this is a female
								print OUTFILE "\t\.\/.";
							}	
						}
						elsif($genotypes[$GT] eq '0/0'){
							if($sexes[$y-9] == 0){ # it is a male
								print OUTFILE "\t".$columns[3]."\/";
							}
							else{ # this is a female
								print OUTFILE "\t$columns[3]\/$columns[3]";
							}	
								
						}
						else{
							print "Something is weird with the invariant genotypes\n";
						}
					}	
					print OUTFILE "\n";
				}
				else{ # This is probably what happens all the time
					#print "AD $AD\n";
					print OUTFILE $columns[0],"\t",$columns[1],"\t",$columns[3];
					for ($y = 9 ; $y <= $#columns; $y++ ) {
						if($sexes[$y-9] == 0){ # this is a male
							@alleledepth=();
							@genotypes=();
							@genotypes=split(":",$columns[$y]);
							@alleledepth=split(",",$genotypes[$AD]);
							@allelieos = split("/",$genotypes[$GT]); # these are the alleles with numbers 0, 1, 2, etc
							@maxcounter=();
							$counter2=0;
							$max=0;
							$max=max @alleledepth;
							# now cycle through each allele depth to find highest and see if there is a tie
							foreach my $alleledepth (@alleledepth){
								if($alleledepth == $max){
									push(@maxcounter,$counter2);
								}
								$counter2+=1;
							}	
							@maxcounter = shuffle @maxcounter;
							if($genotypes[$GT] eq './.'){
								print OUTFILE "\t\.\/";
							}
							elsif($maxcounter[0] == 0){ # the highest depth is the REF allele
								if($columns[3] ne '*'){
									print OUTFILE "\t".$columns[3]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							elsif($maxcounter[0] == 1){ # the highest depth is the first alt allele
								if($altalleles[0] ne '*'){
									print OUTFILE "\t".$altalleles[0]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							elsif($maxcounter[0] == 2){ # the highest depth is the second alt allele
								if($altalleles[1] ne '*'){
									print OUTFILE "\t".$altalleles[1]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							elsif($maxcounter[0] == 3){ # the highest depth is the third alt allele
								if($altalleles[2] ne '*'){
									print OUTFILE "\t".$altalleles[2]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							elsif($maxcounter[0] == 4){ # the highest depth is the fourth alt allele
								if($altalleles[3] ne '*'){
									print OUTFILE "\t".$altalleles[3]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							elsif($maxcounter[0] == 5){ # the highest depth is the fifth alt allele
								if($altalleles[4] ne '*'){
									print OUTFILE "\t".$altalleles[4]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							elsif($maxcounter[0] == 6){ # the highest depth is the fifth alt allele
								if($altalleles[5] ne '*'){
									print OUTFILE "\t".$altalleles[5]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							elsif($maxcounter[0] == 7){ # the highest depth is the fifth alt allele
								if($altalleles[6] ne '*'){
									print OUTFILE "\t".$altalleles[6]."\/";
								}
								else{
									print OUTFILE "\t\.\/";
								}	
							}
							else{
								print OUTFILE "\t\.\/";
							}
						}
						else{ # this is a female
							@genotypes=();
							@genotypes=split(":",$columns[$y]);
							if($genotypes[$GT] eq './.'){
								print OUTFILE "\t\.\/\.";
							}
							elsif($genotypes[$GT] eq '0/0'){
								print OUTFILE "\t".$columns[3]."\/".$columns[3];
							}
							else{ # this female is heterozygous
								@altalleles = split(",",$columns[4]);
								@allelieos = split("/",$genotypes[$GT]);
								# first print first allele for this female
								if($altalleles[$allelieos[0]-1] ne '*'){
									if($allelieos[0] eq '0'){
										print OUTFILE "\t".$columns[3]."\/";
									}
									else{
										print OUTFILE "\t".$altalleles[$allelieos[0]-1]."\/";
									}
								}
								else{
									print OUTFILE "\t\.\/";
								}
								# now print the second allele	
								if($altalleles[$allelieos[1]-1] ne '*'){
									if($allelieos[1] eq '0'){
										print OUTFILE $columns[3];
									}
									else{
										print OUTFILE $altalleles[$allelieos[1]-1];
									}
								}
								else{
									print OUTFILE "\.";
								}
							}

						}	
					}
					print OUTFILE "\n";	
				}
			}
		}# endif
		elsif(substr($columns[0],0,6) eq '#CHROM'){ # print the first line
			print OUTFILE "#CHROM	POS	REF";
				for ($y = 9 ; $y <= $#columns; $y++ ) {
					print OUTFILE "\t",$columns[$y];
				}
				print OUTFILE "\n";			
		}
}# end while
close DATAINPUT;
close OUTFILE;

```

# Convert tab to nexus

This script should be able to handle tab files that have haploid genotypes.  I modified it from (21_tab_to_interleave_nexus.pl) to take as input tab files with a variable number of reference genomes (11_tab_to_interleave_nexus.pl). Note that this script removes indels that are present in the ref relative to the ingroup or vice versa.  This is done because it is coded to enforce only one REF base for each ingroup genotype and then when longer genotypes are encountered in an ingroup taxa, a gap/missing site is used. This way the alignmnet in the nexus file is preserved.

These files are here on goblin:
```
/work/ben/2017_SEAsian_macaques/SEAsian_macaques_bam
```

commandline:
``` bash
./11_tab_to_interleave_nexus.pl ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf.gz.bydepth.tab ../SEAsian_macaques_bam/females/all_chrM_noBSQR_allsites.vcf.gz.bydepth.nxs 1
```
``` perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program '17_adds_outgroup_to_lots_of_tab_files.pl'
#  and generates an interleaved nexus file that includes degenerate bases and gaps

# it is hardcoded to expect one base from two outgroup sequences after the reference (rhesus) seq


# run it like this
# 21_tab_to_interleave_nexus.pl input.tab output_interleave.nxs 1


my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];
my $number_of_reference_genomez = $ARGV[2]; # this indicates how many reference genomes preceed the data

unless (open DATAINPUT, $inputfile) {
	print 'Can not find the input file.\n';
	exit;
}


my @temp;
my @temp1;
my @names;
my %datahash;
my $y;
my $x;
my $watisitnow;
my $count=0;
my $interleave=0;
my $ref_length_check=0;

# Read in datainput file

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		@names=@temp;
		for ($y = 2; $y <= $#names; $y++ ) {
			$names[$y] =~ s/-//gi; # get rid of dashes in names
			$datahash{$names[$y]}=''; # initialize the hash
		}	
		# print preamble to output file
		unless (open(OUTFILE, ">$outputfile"))  {
			print "I can\'t write to $outputfile\n";
			exit;
		}
		print "Creating output file: $outputfile\n";
		print OUTFILE "#NEXUS\n\n";
		print OUTFILE "BEGIN DATA\;\nDIMENSIONS NTAX=",$#names-1," NCHAR= XXXX\;\n";
		print OUTFILE "FORMAT DATATYPE=DNA  MISSING=? INTERLEAVE GAP=- \;\n";
		print OUTFILE "MATRIX\n";
		print OUTFILE "\n";
	}
	else{	
		# only print ones that are not microsats or indels in the outgroups
		if($interleave>79){  # print this section of the data
			for ($y = 2; $y <= $#names; $y++ ) {
				print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
			}
			print OUTFILE "\n";
			$interleave=0;
			# clear the hash
			for ($y = 2 ; $y <= $#names; $y++ ) {
				$datahash{$names[$y]}='';
			}	
		}
		$ref_length_check=0;
		for ($y = 2 ; $y < 2+$number_of_reference_genomez; $y++ ) {
		    if(length($temp[$y]) != 1){
			$ref_length_check=1;
		       }
		}
	        if($ref_length_check == 0){ # all references have length of 1, so process the site as usual
			$count=$count+1; # this is the count of all positions
			$interleave=$interleave+1; # this is the count of the interleave length
			
			
			# now add data to the hash
			for ($y = 2 ; $y < 2+$number_of_reference_genomez; $y++ ) {	# first the three outgroups which are haploid
				if($temp[$y] ne '*'){
					$datahash{$names[$y]} = $datahash{$names[$y]}.uc($temp[$y]);
				}
				else{
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}	
			}
			for ($y = (2+$number_of_reference_genomez) ; $y <= $#temp; $y++ ) { # now the ingroups, which are diploid, usually (except chrX and chrY)
				# for these, we need to use IUPAC codes
				if(($temp[$y] eq 'G/G')||($temp[$y] eq 'C/C')||($temp[$y] eq 'T/T')||($temp[$y] eq 'A/A')||($temp[$y] eq 'G/')||($temp[$y] eq 'C/')||($temp[$y] eq 'T/')||($temp[$y] eq 'A/')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.substr($temp[$y],0,1);
				}
				elsif(($temp[$y] eq './.')||($temp[$y] eq './')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}
				elsif(($temp[$y] eq 'C/T')||($temp[$y] eq 'T/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'Y';
				}
				elsif(($temp[$y] eq 'A/G')||($temp[$y] eq 'G/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'R';
				}
				elsif(($temp[$y] eq 'A/C')||($temp[$y] eq 'C/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'M';
				}
				elsif(($temp[$y] eq 'A/T')||($temp[$y] eq 'T/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'W';
				}
				elsif(($temp[$y] eq 'C/G')||($temp[$y] eq 'G/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'S';
				}
				elsif(($temp[$y] eq 'G/T')||($temp[$y] eq 'T/G')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'K';
				}
				else{ # this is a microsat, so substitute missing data
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}
			}	
		}
	}
}		

# print the last line
for ($y = 2; $y <= $#names; $y++ ) {
	print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
}
#print OUTFILE "\n";



print OUTFILE "\;\nEND\;";
print OUTFILE "\n";
print OUTFILE "\n";


close OUTFILE;
print "The number of sites is $count\n";

# now update the number of bases

my $status;
$status = system("perl -p -i -e 's/XXXX/$count/g' $outputfile");

```

# Another tab to nexus script

This one assumes haploid genotypes and includes indels in the ref and the ingroup:
```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );


#  This program reads in a tab delimited genotype file generated
#  by the perl program '17_adds_outgroup_to_lots_of_tab_files.pl'
#  and generates an interleaved nexus file that includes degenerate bases and gaps

# it can accommodate any number of haploid reference (rhesus) genomes before the ingroup

# it will include haploid indels in the ref or ingroup

# so this should be used for mtDNA or yDNA generated either by gatk calling haploid genotypes
# or by my genotype_by_depth script where all individuals are called as haploid

# I am restricting it to haploids (at least for now). 

# probably a good idea to check alignments manually (or convert to fasta and do it with a program)
# because the indels will not be properly aligned

# run it like this
# 21_tab_to_interleave_nexus.pl input.tab output_interleave.nxs 1


my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];
my $number_of_reference_genomez = $ARGV[2]; # this indicates how many reference genomes preceed the data

unless (open DATAINPUT, $inputfile) {
	print 'Can not find the input file.\n';
	exit;
}


my @temp;
my @temp1;
my @names;
my %datahash;
my $z;
my $y;
my $x;
my $watisitnow;
my $count=0;
my $interleave=0;
my @largest_length=1;
my $max;

# Read in datainput file

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		@names=@temp;
		for ($y = 2; $y <= $#names; $y++ ) {
			$names[$y] =~ s/-//gi; # get rid of dashes in names
			$datahash{$names[$y]}=''; # initialize the hash
		}	
		# print preamble to output file
		unless (open(OUTFILE, ">$outputfile"))  {
			print "I can\'t write to $outputfile\n";
			exit;
		}
		print "Creating output file: $outputfile\n";
		print OUTFILE "#NEXUS\n\n";
		print OUTFILE "BEGIN DATA\;\nDIMENSIONS NTAX=",$#names-1," NCHAR= XXXX\;\n";
		print OUTFILE "FORMAT DATATYPE=DNA  MISSING=? INTERLEAVE GAP=- \;\n";
		print OUTFILE "MATRIX\n";
		print OUTFILE "\n";
	}
	else{	
		if($interleave>=80){  # print this section of the data
			# print out the first 80bp, trim this out of the hash, and move on
			# we assume here that indels never are longer than 80 bp
			# if they were, I think the last interleave will be longer than the others
			for ($y = 2; $y <= $#names; $y++ ) {
				print OUTFILE $names[$y],"\t\t\t",substr($datahash{$names[$y]},0,79),"\n";
			}
			print OUTFILE "\n";
			$interleave=0;
			# clear the hash
			for ($y = 2 ; $y <= $#names; $y++ ) {
				#$datahash{$names[$y]}='';
				substr($datahash{$names[$y]},0,79,"")
			}
			$interleave=length($datahash{$names[2]});
		}
		@largest_length=();
		# store length each genotype in ref
		for ($y = 2 ; $y < 2+$number_of_reference_genomez; $y++ ) {
			$largest_length[$y-2] = length($temp[$y]);
		}
		# and in the ingrpup
		for ($y = (2+$number_of_reference_genomez) ; $y <= $#temp; $y++ ) {
			@temp1=split('\/',$temp[$y]);
			$largest_length[$y-2] = length($temp1[0]);
		}	
		# now add data to the hash
		for ($y = 2 ; $y < 2+$number_of_reference_genomez; $y++ ) {	# first the three outgroups which are haploid
			if($temp[$y] ne '*'){
				$datahash{$names[$y]} = $datahash{$names[$y]}.uc($temp[$y]);
			}
			else{
				$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
			}	
		}
		for ($y = (2+$number_of_reference_genomez) ; $y <= $#temp; $y++ ) { 
			# now the ingroups, which are haploid
			@temp1=split('\/',$temp[$y]);
			if(($temp1[0] ne '*')&&($temp1[0] ne '.')&&($temp1[0] !~ /NON_REF/)){
				$datahash{$names[$y]} = $datahash{$names[$y]}.$temp1[0];
			}	
			else{ # this is a microsat, so substitute missing data
				$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
			}
		}	
		# add gaps if necessary
		$max = max @largest_length;
		if($max >1){
			print $max,"\n";
		}
		for ($y = 2 ; $y <= $#temp; $y++ ) {
			for ($z = $largest_length[$y-2]; $z < $max; $z++ ) {
				#add gaps to hash
				$datahash{$names[$y]}=$datahash{$names[$y]}.'-';
			}
		}
		$count=$count+$max; # this is the count of all positions
		$interleave=$interleave+$max; # this is the count of the interleave length

	}
}		

# print the last line
for ($y = 2; $y <= $#names; $y++ ) {
	print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
}
#print OUTFILE "\n";



print OUTFILE "\;\nEND\;";
print OUTFILE "\n";
print OUTFILE "\n";


close OUTFILE;
print "The number of sites is $count\n";

# now update the number of bases

my $status;
$status = system("perl -p -i -e 's/XXXX/$count/g' $outputfile");
```



