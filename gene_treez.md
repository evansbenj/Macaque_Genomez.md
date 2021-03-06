# Gene treez
directory with full genotype files (not just SNPS):
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio
```
Get the coordinates of some interesting genes:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'MRPL1'
```
Here they are:
```
NDUFAF3  chr03	104846075	104847927
NDUFA2 chr05	138178553	138180781
TACO1 chr17	56920802	56928897
HARS2 chr05	138224999	138233197
COX16 chr14	132330617	132354632
C11orf83   chr11	11484179	11486231
MRPS21  chr01	124258900	124272013
UQCRC1   chr03	104359835	104371673
ATP5B    chr12	54910245	54920096
MRPL48  chr11	65083808	65165565
NDUFB11   chrX	47194137	47196604
NDUFB5  chr03	84594254	84613076
MRPL1  chr04	56478980	56572561
MRPL47 chr03	84577965	84594079
```
Extract section of gene with vcftools:
```
module load StdEnv/2020 vcftools/0.1.16
vcftools --gzvcf FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered1.vcf.gz --chr chr01 --from-bp 124258900 --to-bp 124272013 --out MRPS21.vcf.gz --recode
```
compress the vcf filez
```
module load nixpkgs/16.09  intel/2016.4 tabix/0.2.6
bgzip -c TACO1.vcf.gz.recode.vcf > TACO1.vcf.gz.recode.vcf.gz
```
convert to tab
```
zcat MRPL47.vcf.gz.recode.vcf.gz | vcf-to-tab > MRPL47.tab
```
Convert to nexus
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program '17_adds_outgroup_to_lots_of_tab_files.pl'
#  and generates an interleaved nexus file that includes degenerate bases and gaps

# it is hardcoded to expect one base from two outgroup sequences after the reference (rhesus) seq


# run it like this
# 21_tab_to_interleave_nexus.pl input.tab output_interleave.nxs


my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

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
		if((length($temp[2]) == 1)){ # the ref seq is single bp
			$count=$count+1; # this is the count of all positions
			$interleave=$interleave+1; # this is the count of the interleave length
			
			
			# now add data to the hash
			for ($y = 2 ; $y <= 2; $y++ ) {	# first refseq which is haploid
				if($temp[$y] ne '*'){
					$datahash{$names[$y]} = $datahash{$names[$y]}.uc($temp[$y]);
				}
				else{
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}	
			}
			for ($y = 3 ; $y <= $#temp; $y++ ) { # now the ingroups, which are diploid, usually (except chrX and chrY)
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


print OUTFILE "BEGIN PAUP\;\n";
print OUTFILE "OutGroup download\;\n";
print OUTFILE "END\;\n";
#print OUTFILE "Lset  nst=6  rates=invgamma\;\n";
#print OUTFILE "mcmc ngen=2000000 savebrlens=yes\;\n";
#print OUTFILE "sumt burnin=10000\;\n";
#print OUTFILE "quit\;\n";

close OUTFILE;
print "The number of sites is $count\n";

# now update the number of bases
my $status;
$status = system("perl -p -i -e 's/XXXX/$count/g' $outputfile");

```
