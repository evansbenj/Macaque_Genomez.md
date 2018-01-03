# Polymorphism and divergence

I have a script called "Boot_from_tab_diverge_poly_2015.pl" that calculates polymorphism and divergence statistics from a tab delimited file. For each chr, this is the commandline:

For autosomes in this format: #CHROM	POS	REF	SAMN03083651	SAMN03264597	SAMN03264605	SAMN03264609	SAMN03264610	SAMN03264614	SAMN03264646	SAMN03264649	SAMN03264650	SAMN
03264651	SAMN03264652	SAMN03264653	SAMN03264658	SAMN03264659	SAMN03264660	SAMN03264661	SAMN03264662	SAMN03264666	SAMN03264667	SAMN03264673
	SAMN03264674	SAMN03264677	SAMN03264678	SAMN03264683	SAMN03264685	SAMN03264696	SAMN03264703	SAMN03264706	SAMN03264707	SAMN03264708	SAMN
03264709	SAMN03264710	SAMN03264712	SAMN03264713	SAMN03264715	SAMN03264718	SAMN03264719	SAMN03264726	SAMN03264733	SAMN03264734 use this:

```
Boot_from_tab_diverge_poly_2015.pl /work/ben/2017_rhesus_genomez/F_and_M/FandM_chr01_BSQR_jointgeno_allsites_filtered.vcf.gz.tab 1100000100011010101001111110110101001011 3_4_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39_40 rhesus_chr01_poly_and_diverge.txt
```

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
		$commandline = "Boot_from_tab_diverge_poly_2015.pl ".$_." ".$bin_sex." 3_4_".$numberz[$y]." ".$species[$y]."_".$what_what[0]."_boot.poly";
		$status = system($commandline);
	}
}
```
