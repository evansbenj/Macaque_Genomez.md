# It would be neat to look at the OXPHOS genes
Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general
```

Here's how I got the coordinates of all CDS in exons in the MacaM genome:
```
grep 'CDS' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | cut -f1,4,5 > coordinates_all_CDS.txt
```
Or only for chrX:
```
grep 'CDS' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | grep 'chrX' | cut -f1,4,5 > coordinates_chrX_all_CDS.txt
```

Here's how I got the coordinates of all CDS in OXPHOS genes:
```
grep 'CDS' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'COX|NDUF|UQCR|ATP5|CYC1|SDHB|SDHA|SDHC|SDHD' | cut -f1,4,5 > coordinates_OXPHOS_all_CDS.txt
```
Or only CDS of OXPHOS genes on only chrX:
```
grep 'CDS' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'COX|NDUF|UQCR|ATP5|CYC1|SDHB|SDHA|SDHC|SDHD' | grep 'chrX'| cut -f1,4,5 > coordinates_chrX_OXPHOS_all_CDS.txt
```
Heres how I got the final list of genes in complex 1 based on this paper in Essays in Biochemistry (2018) 62 255-270:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'NDUFA2|NDUFAF6|NDUFS1|NDUFV1|NDUFV2|NUBPL' | cut -f1,4,5 > coordinates_OXPHOS_Nmodule_complex1.txt

grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'NDUFA3|NDUFA8|NDUFA13|NDUFAF5|NDUFAF7|NDUFS2|NDUFS3|NDUFS7|NDUFS8'| cut -f1,4,5 > coordinates_OXPHOS_Qmodule_complex1.txt

grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'NDUFAF3|NDUFAF4|TIMMDC1'| cut -f1,4,5 > coordinates_OXPHOS_QNmodule_complex1.txt

grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'ATP5SL|NDUFB6|FOXRED1|NDUFB10|NDUFB11|NDUFB4|NDUFB5|TMEM70' | cut -f1,4,5 > coordinates_OXPHOS_ND4module_complex1.txt

grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'ACAD9|=COA1|ECSIT|NDUFAF1|NDUFC1|NDUFC2|TMEM126B|TMEM186' | cut -f1,4,5 > coordinates_OXPHOS_ND2module_complex1.txt

grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'C9orf123|NDUFAB1|NDUFB2|NDUFB3|NDUFB7|NDUFB8|NDUFB9' | cut -f1,4,5 > coordinates_OXPHOS_ND5module_complex1.txt

grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'NDUFA6|NDUFA7|NDUFA12|NDUFS4|NDUFS6|NDUFV3' | cut -f1,4,5 > coordinates_OXPHOS_other_complex1.txt
```
And here are the coordinates of nuclear genes in complex 1 (with names; I could not find three of them in the gff file):
```
NDUFA2	1	chr05	138178553	138180781
NDUFAF6	1	chr08	93359100	93394552
NDUFS1	1	chr02b	93979633	94015898
NDUFV1	1	chr11	6690185	6696071
NDUFV2	1	* not in gff file		
NUBPL	1	chr14	92789615	93049799
NDUFA3	1	chr19	9634493	9638249
NDUFA5	1	* not in gff file		
NDUFA8	1	chr09	94025131	94040826
NDUFA13	1	chr19	19317399	19330902
NDUFAF5	1	chr15	46796983	46833911
NDUFAF7	1	chr02a	37732212	37749061
NDUFS2	1	chr01	135484416	135497548
NDUFS3	1	chr11	18260513	18265960
NDUFS7	1	chr19	1153942	1165810
NDUFS8	1	chr11	6532924	6539484
NDUFAF3	1	chr03	104846075	104847927
NDUFAF4	1	chr06	95211059	95218436
TIMMDC1	1	chr03	157411293	157436295
ATP5SL	1	chr19	36878373	36887017
NDUFB6	1	chr09	59806087	59825964
FOXRED1	1	chr11	118312703	118322904
NDUFB10	1	chr16	1930096	1932930
NDUFB11	1	chrX	47194137	47196604
NDUFB1	1	* not in gff file		
NDUFB4	1	chr03	156309713	156315708
NDUFB5	1	chr03	84594254	84613076
TMEM70	1	chr08	72064967	72071542
ACAD9	1	chr03	148084830	148114936
COA1	1	chr07	64839692	64850427
ECSIT	1	chr19	11316389	11345234
NDUFAF1	1	chr14	17717208	17731422
NDUFC1	1	chr04	139226852	139232639
NDUFC2	1	chr11	69474342	69546903
TMEM126B	1	chr11	77465139	77472013
TMEM186	1	chr16	8788484	8790936
C9orf123 (DMAC1)	1	chr09	34632092	34635355
NDUFAB1	1	chr16	22236971	22249112
NDUFB2	1	chr07	166649885	166659561
NDUFB3	1	chr02b	88692897	88705887
NDUFB7	1	chr19	14340049	14345977
NDUFB8	1	chr10	95956681	95962600
NDUFB9	1	chr08	123481402	123491824
NDUFA6	1	chr15	84086730	84096073
NDUFA7	1	chr19	8302280	8312129
NDUFA12	1	chr12	94150466	94193941
NDUFS4	1	chr05	54356765	54473046
NDUFS6	1	chr05	1601204	1616216
NDUFV3	1	chr07	3770437	3785385
```
Here's how I got the final list of genes in complex 2:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'SDHA|SDHAF1|SDHAF2|SDHB|SDHC|SDHD|ACN9|C6orf57' | cut -f1,4,5 > coordinates_OXPHOS_complex2.txt
```
here are the coordinates for complex2:
```
	SDHA	2	chr05	191800	225299
	SDHAF1	2	chr19	31551177	31552175
	SDHAF2	2	chr11	12844602	12863522
	ACN9 (SDHAF3)	2	chr07	122255513	122331494
	C6orf57 (SDHAF4)	2	chr06	68582810	68605348
	SDHB	2	chr01	15628447	15662313
	SDHC	2	chr01	135601364	135647524
	SDHD	2	chr11	104210163	104219988
```

Here's how I got the final list of genes in complex 3:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'UQCRB|UQCRQ|UQCRC1|UQCRC2|CYC1|UQCRH_|UQCR10|UQCR11|BCS1L|TTC19|UQCC|MNF1|C11orf83' | cut -f1,4,5 > coordinates_OXPHOS_complex3.txt
```
genes in comlex3:
```
Complex3					
CYB (32 aas)	UQCC (UQCC1, Cbp3)	3	chr15	28688308	28798664
	MNF1 (UQCC2,Cbp6)	3	chr06	34463088	34477872
	C11orf83 (UQCC3,Cbp4)	3	chr11	11484179	11486231
	UQCRB	3	chr08	94674995	94684605
	UQCRQ	3	chr05	130350642	130352620
	UQCRC1	3	chr03	104359835	104371673
	UQCRC2	3	chr16	20500148	20530512
	CYC1	3	chr08	143057509	143059999
	UQCRH	3	chr01	45636133	45647963
	UQCR10	3	chr15	71526832	71530391
	UQCR11	3	chr19	1385443	1392486
	BCS1L	3	chr02b	106853441	106857546
	MZM1L (LYRM7)	3	* not in gff file		
	UQCRFS1	3	* not in gff file		
	TTC19	3	chr17	16554117	16582008
 ```

here's the genes in complex 4:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'TACO1|COX14|COA3|CMC1|COA1|SURF1|C4orf52|COX17|COX11|COX19|COX10|COX15|COX4I1|COX5A|HIGD1A|COX17|COX16|COA6|SCO2|=SCO1|COX20|COX18|TMEM177|COX5B|COX6C|COX7B|COX7C|COX8A|=MR1_|PET100|PET117|COX6A1|COX6B1|COX7A1|NDUFA4' | cut -f1,4,5 > coordinates_OXPHOS_complex4.txt
```
coordinates:
```
TACO1	4	chr17	56920802	56928897	
COX14	4	chr12	48300943	48310011	
COA3	4	chr17	50660415	50661510	
CMC1	4	chr03	19605893	19685329	
COA1	4	chr07	64839692	64850427	
SURF1	4	chr09	106058130	106063340	
C4orf52 (MITRAC7)	4	chr04	25279951	25295066	
COX17	4	chr03	157252481	157260548	* possibly with CO1 and CO2
COX11	4	chr17	38276744	38283435	
COX19	4	chr07	34825571	34833611	
COX10	4	chr17	14089633	14230045	
COX15	4	chr10	95151638	95177391	
COX4I1	4	chr16	70395661	70403728	
COX5A	4	chr14	52057142	52077913	
HIGD1A	4	chr03	98375417	98385870	
COX16	4	chr14	132330617	132354632	
COA6	4	chr01	210102822	210113992	
SCO2	4	chr15	92380326	92383164	
SCO1	4	chr17	10622390	10637333	
COX20	4	chr01	220779458	220788687	
COX18	4	chr04	61522247	61536172	
TMEM177	4	chr02b	11434587	11435874	
COX5B	4	chr02a	92746891	92749017	
COX6C	4	chr08	98367530	98383414	
COX7B	4	chrX	71633448	71639804	
COX7C	4	chr05	56278755	56279160	
COX8A	4	chr11	8907900	8910660	
MR1 (MR1S)	4	chr01	185833218	185860646	
PET100	4	chr19	7673574	7675433	
PET117	4	chr15	51279255	51284434	
COX6A1 (COX6A)	4	chr12	120032059	120034809	
COX6B1	4	chr19	31207304	31214924	
COX7A	4	chr19	31714988	31716937	
NDUFA4	4	chr07	97408775	97416979	
```

Here's how I got the final list of genes in complex 5 based on this paper Essays in Biochemistry (2018) 62 255–270:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'ATPA|TMEM70|ATPIF1|USMG5|C14orf2_|ATP5O|ATP5F1|ATP5I|ATP5L|ATP5J2|ATP5C1|ATP5D|ATP5E|ATP5B|ATP5A1|ATP5G|ATP5J|ATP5H' | cut -f1,4,5 > coordinates_OXPHOS_complex5.txt
```
Here the coordinates are (with names):
```
ATPAF1	5	chr01 45953218	45987492
TMEM70	5	chr08 72064967	72071542
ATPAF2	5	chr17	18045894	18066877
ATPIF1 (IF1)	5	chr01	26938483	26940445
USMG5 (ATP5MD, DAPIT)	5	chr10	99032999	99036509
C14orf2 (ATP5MPL, MP68, PLPM, MLQ,6,8PL)	5	chr14	166432920	166436127
ATP5O (OSCP)	5	chr07	12722004	12734438
ATP5C1	5	chr10	7723482	7740495
ATP5D	5	chr19	1002804	1005971
ATP5E	5	chr15	5385192	5389178
ATP5B	5	chr12	54910245	54920096
ATP5A1	5	chr18	35035298	35047335
ATPG1	5	chr17	44211923	44214912
ATPG2	5	chr12	51828101	51839397
ATPG3	5	chr02b	62111164	62116182
ATP5F1	5	chr01	112030533	112043548
ATP5H (ATPH)	5	chr17	68443001	68447706
ATP5J	5	chr07	21199179	21209971
ATP5I	5	chr04	634919	636858
ATP5L	5	chr11	110459798	110468428
ATPJ2	5	chr07	124582469	124590730
```

The cool thing is that the two mt genes in complex 5 (ATP6, ATP8) interact with others in the peripheral stalk, which includes 3 genes (alternative acronyms are in parentheses): ATPIF1 (IF1), USMG5 (ATP5MD, DAPIT), C14orf2 (ATP5MPL, MP68, PLPM, MLQ,6,8PL) so we have a clear hypothesis here.

Here's a preliminary list of all the autosomally encoded proteins in in the OXPHOS complex based on this paper: Genome Res. 2018. 28: 952-967:

This gets all the OXPHOS genes:
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'NDUFA2|NDUFAF6|NDUFS1|NDUFV1|NDUFV2|NUBPL|NDUFA3|NDUFA8|NDUFA13|NDUFAF5|NDUFAF7|NDUFS2|NDUFS3|NDUFS7|NDUFS8|NDUFAF3|NDUFAF4|TIMMDC1|ATP5SL|NDUFB6|FOXRED1|NDUFB10|NDUFB11|NDUFB4|NDUFB5|TMEM70|ACAD9|=COA1|ECSIT|NDUFAF1|NDUFC1|NDUFC2|TMEM126B|TMEM186|C9orf123|NDUFAB1|NDUFB2|NDUFB3|NDUFB7|NDUFB8|NDUFB9|NDUFA6|NDUFA7|NDUFA12|NDUFS4|NDUFS6|NDUFV3|SDHA|SDHAF1|SDHAF2|SDHB|SDHC|SDHD|ACN9|C6orf57|UQCRB|UQCRQ|UQCRC1|UQCRC2|CYC1|UQCRH_|UQCR10|UQCR11|BCS1L|TTC19|UQCC|MNF1|C11orf83|TACO1|COX14|COA3|CMC1|COA1|SURF1|C4orf52|COX17|COX11|COX19|COX10|COX15|COX4I1|COX5A|HIGD1A|COX17|COX16|COA6|SCO2|=SCO1|COX20|COX18|TMEM177|COX5B|COX6C|COX7B|COX7C|COX8A|=MR1_|PET100|PET117|COX6A1|COX6B1|COX7A1|NDUFA4|ATPA|TMEM70|ATPIF1|USMG5|C14orf2_|ATP5O|ATP5F1|ATP5I|ATP5L|ATP5J2|ATP5C1|ATP5D|ATP5E|ATP5B|ATP5A1|ATP5G|ATP5J|ATP5H' | egrep 'transcript_01' > OXPHOS.txt
```
This gets all the non-OXPHOS genes (using the -v after egrep)
```
grep 'mRNA' ~/projects/rrg-ben/ben/2017_SEAsian_macaques/MacaM/MacaM_Rhesus_Genome_Annotation_v7.6.8.gff | egrep 'NDUFA2|NDUFAF6|NDUFS1|NDUFV1|NDUFV2|NUBPL|NDUFA3|NDUFA8|NDUFA13|NDUFAF5|NDUFAF7|NDUFS2|NDUFS3|NDUFS7|NDUFS8|NDUFAF3|NDUFAF4|TIMMDC1|ATP5SL|NDUFB6|FOXRED1|NDUFB10|NDUFB11|NDUFB4|NDUFB5|TMEM70|ACAD9|=COA1|ECSIT|NDUFAF1|NDUFC1|NDUFC2|TMEM126B|TMEM186|C9orf123|NDUFAB1|NDUFB2|NDUFB3|NDUFB7|NDUFB8|NDUFB9|NDUFA6|NDUFA7|NDUFA12|NDUFS4|NDUFS6|NDUFV3|SDHA|SDHAF1|SDHAF2|SDHB|SDHC|SDHD|ACN9|C6orf57|UQCRB|UQCRQ|UQCRC1|UQCRC2|CYC1|UQCRH_|UQCR10|UQCR11|BCS1L|TTC19|UQCC|MNF1|C11orf83|TACO1|COX14|COA3|CMC1|COA1|SURF1|C4orf52|COX17|COX11|COX19|COX10|COX15|COX4I1|COX5A|HIGD1A|COX17|COX16|COA6|SCO2|=SCO1|COX20|COX18|TMEM177|COX5B|COX6C|COX7B|COX7C|COX8A|=MR1_|PET100|PET117|COX6A1|COX6B1|COX7A1|NDUFA4|ATPA|TMEM70|ATPIF1|USMG5|C14orf2_|ATP5O|ATP5F1|ATP5I|ATP5L|ATP5J2|ATP5C1|ATP5D|ATP5E|ATP5B|ATP5A1|ATP5G|ATP5J|ATP5H' | egrep -v 'transcript_01' > nonOXPHOS.txt
```

For chrX PopGenome needed a "pseudodiploid" vcf file.  In this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data
```
I made one like this:
```
zcat all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode.vcf.gz > temp.vcf
cat temp.vcf | sed 's/   0:/     0\/0:/g' > temp2.vcf
cat temp2.vcf | sed 's/  1:/     1\/1:/g' > temp3.vcf
cat temp3.vcf | sed 's/  .:/     .\/.:/g' > temp4.vcf
module load nixpkgs/16.09  intel/2016.4 tabix/0.2.6
bgzip -c temp4.vcf > temp4.vcf.gz
tabix -p vcf temp4.vcf.gz
mv temp4.vcf.gz chrX_pseudo_diploid_forPopGenome.vcf.gz
mv temp4.vcf.gz.tbi chrX_pseudo_diploid_forPopGenome.vcf.gz.tbi
```
And here is the R code that gives me S and NS polymorphism and divergence by population (!):
```
## Working directory
# https://cran.r-project.org/web/packages/PopGenome/vignettes/Whole_genome_analyses_using_VCF_files.pdf
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2017_SEAsian_macaque_genomz/PopGenome_OXPHOS")
library(ggplot2)
library(PopGenome)

#bru_PF707
#hecki_PF505
#hecki_PF643
#hecki_PF644
#hecki_PF647
#hecki_PF648
#maura_PF615
#maura_PF713
#maura_PM613
#maura_PM614
#maura_PM616
#nem_GumGum_female
#nem_Ngsang_sumatra_female
#nem_PM1206
#nem_PM664
#nem_PM665
#nem_Sukai_male
#nigra_PF1001
#nigra_PF660
#nigra_PM1003
#nigrescens_PM1011
#nigrescens_PM654
#tog_PF549
#tonk_PF511
#tonk_PF559
#tonk_PF563
#tonk_PF597
#tonk_PF626
#tonk_PM592

#chr01	225002135
#chr02a	108967917
#chr02b	131519175
#chr03	198060209
#chr04	190369981
#chr05	179725205
#chr06	171866349
#chr07	185267708
#chr08	144034664
#chr09	111027318
#chr10	129655328
#chr11	127102482
#chr12	133317794
#chr13	95368959
#chr14	169736342
#chr15	92674614
#chr16	74750809
#chr17	76917410
#chr18	70128972
#chr19	53113586
#chrX	148935249

#OXPHOS coordinates (all chr except chr13)
#chr01	15628447	15662313
#chr01	37992642	37997542
#chr01	45636133	45647963
#chr01	112030533	112043548
#chr01	135484416	135497548
#chr01	135601364	135647524
#chr01	220779458	220788687
#chr02a	37732212	37749061
#chr02a	42997182	43010254
#chr02a	92746891	92749017
#chr02b	62111164	62116182
#chr02b	88692897	88705887
#chr02b	93979633	94015898
#chr02b	128742129	128807570
#chr03	84594254	84613076
#chr03	104359835	104371673
#chr03	104846075	104847927
#chr03	114610662	114642316
#chr03	156309713	156315708
#chr03	156309713	156315708
#chr03	157252481	157260548
#chr04	634919	636858
#chr04	8603394	8653763
#chr04	46628495	46628917
#chr04	61522247	61536172
#chr04	97928438	97928963
#chr04	139226852	139232639
#chr05	191800	225299
#chr05	1601204	1616216
#chr05	54356765	54473046
#chr05	56278755	56279160
#chr05	58728795	58927343
#chr05	130350642	130352620
#chr05	138178553	138180781
#chr06	73250180	73256003
#chr06	95211059	95218436
#chr07	3770437	3785385
#chr07	12722004	12734438
#chr07	21199179	21209971
#chr07	34825571	34833611
#chr07	97408775	97416979
#chr07	124582469	124590730
#chr07	166649885	166659561
#chr08	93359100	93394552
#chr08	94674995	94684605
#chr08	94674995	94684605
#chr08	98367530	98383414
#chr08	123481402	123491824
#chr08	143057509	143059999
#chr09	59806087	59825964
#chr09	59806087	59825964
#chr09	94025131	94040826
#chr10	7723482	7740495
#chr10	95151638	95177391
#chr10	95151638	95177391
#chr10	95956681	95962600
#chr11	6532924	6539484
#chr11	6690185	6696071
#chr11	8907900	8910660
#chr11	12844602	12863522
#chr11	18260513	18265960
#chr11	69474342	69546903
#chr11	69534557	69546903
#chr11	69534557	69546903
#chr11	69534557	69546903
#chr11	104210163	104219988
#chr11	110459798	110468428
#chr12	4801409	4840545
#chr12	48300943	48310011
#chr12	51828101	51839397
#chr12	54910245	54920096
#chr12	55497421	55499859
#chr12	94150466	94193941
#chr12	94150482	94193941
#chr12	120032059	120034809
#chr14	17717208	17731422
#chr14	52057142	52077913
#chr14	112154408	112167663
#chr14	132329422	132413355
#chr14	132330617	132354632
#chr15	5385192	5389178
#chr15	32433903	32443159
#chr15	46796983	46833911
#chr15	71526832	71530391
#chr15	71526832	71530391
#chr15	84086730	84096073
#chr16	1930096	1932930
#chr16	20500148	20530512
#chr16	22236971	22249112
#chr16	29233433	29234114
#chr16	70395661	70403728
#chr17	14089633	14230045
#chr17	38276744	38283435
#chr17	44211923	44214912
#chr17	68443001	68447706
#chr17	69386150	69427976
#chr18	35035298	35047335
#chr19	1002804	1005971
#chr19	1153942	1165810
#chr19	1385443	1392486
#chr19	5934304	5942627
#chr19	8302280	8312129
#chr19	9634493	9638249
#chr19	14340049	14345977
#chr19	19317399	19330902
#chr19	31207304	31214924
#chr19	31551177	31552175
#chr19	31714988	31716937
#chr19	36878373	36887017
#chr19	36878373	36887017
#chr19	36878373	36886666
#chr19	50026620	50032055
#chrX	47194137	47196604
#chrX	71633448	71639804
#chrX	113497709	113502918

# read the vcf file
# GENOME.class <- readVCF("FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz",
GENOME.class <- readVCF("temp4.vcf.gz",
                        numcols=10000,
                        ##########
                        tid="chrX",
                        from=1, 
                        ##########
                        #to= 225002135, # for chr01
                        #to= 108967917, # for chr02a
                        #to= 131519175, # for chr02b
                        #to= 198060209, # for chr03
                        #to= 190369981, # for chr04
                        #to= 179725205, # for chr05
                        #to= 171866349, # for chr06
                        #to= 185267708, # for chr07
                        to= 148935249, # for chrX
                        approx=FALSE, 
                        out="", 
                        parallel=FALSE, 
                        gffpath="MacaM_Rhesus_Genome_Annotation_v7.6.8.gff")
# define the populations
GENOME.class <- set.populations(GENOME.class,
                                list(
                                  c("bru_PF707"),
                                  c("hecki_PF505","hecki_PF643","hecki_PF644","hecki_PF647","hecki_PF648"),
                                  c("maura_PF615","maura_PF713","maura_PM613","maura_PM614","maura_PM616"),
                                  c("nem_Ngsang_sumatra_female"),
                                  c("nem_GumGum_female","nem_PM1206","nem_PM664","nem_PM665","nem_Sukai_male"),
                                  c("nigra_PF1001","nigra_PF660","nigra_PM1003"),
                                  c("nigrescens_PM1011","nigrescens_PM654"),
                                  c("tog_PF549"),
                                  c("tonk_PF511","tonk_PF559","tonk_PF563","tonk_PF597","tonk_PF626","tonk_PM592")), 
                                diploid=TRUE)

# verify the SNPs (there is a typo on pg 8 of the manual about this command)
# Set syn & nonsyn SNPs: The results are stored
# in the slot GENOME.class@region.data@synonymous
# The input of the set.synnonsyn function is an object of
# class GENOME and a reference chromosome in FASTA format.
GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="MacaM_mt_female.fa",save.codons=T)

# this is a list of whether each SNP is synonymous (1) or not (0)
# GENOME.class@region.data@synonymous

# this is the number of synonymous SNPs:
sum(GENOME.class@region.data@synonymous[[1]]==1, na.rm=TRUE)
# this is the number of nonsynonymous SNPs:
sum(GENOME.class@region.data@synonymous[[1]]==0, na.rm=TRUE)
# Here, we have to define the parameter na.rm=TRUE because NaN values in this slot
# indicate that the observed SNP is in a non-coding region
# We now could split the data into gene regions



##########
genePos <- get_gff_info(gff.file="MacaM_Rhesus_Genome_Annotation_v7.6.8.gff",chr="chr07", 
                        feature="mRNA")
genes <- splitting.data(GENOME.class, positions=genePos, type=2)

# Now we perform The Tajima’s D statistic on the whole data set and
# consider only nonsyn SNPs in each gene/region.
#genes <- neutrality.stats(genes, subsites="nonsyn", FAST=TRUE)
#nonsynTaj <- genes@Tajima.D

# The same now for synonymous SNPs
#genes <- neutrality.stats(genes, subsites="syn", FAST=TRUE)
#synTaj <- genes@Tajima.D

# To have a look at the differences of syn and nonsyn Tajima D values in each gene we
# could do the following plot:
# plot(nonsynTaj, synTaj, main="2L: Genes : Tajima’s D ")

# McDonald Kreitman test
# Verify the syn/non-syn SNPs
# See # verify the SNPs (there is a typo on pg 8 of the manual about this command) ABOVE
# Set the populations
# See # define the populations ABOVE 
# Splitting the data into genes
# See # We now could split the data into gene regions ABOVE
# Peform the MKT
genes <- MKT(genes)
# lengths(genes@MKT,use.names = T)
# To look at gene 9 we can do the following:
# genes@MKT[[9]]


# OK let's figure out which genes are OXPHOS genes
# gff_info <- get_gff_info(genes, position=3, chr="chr19", gff.file="MacaM_Rhesus_Genome_Annotation_v7.6.8_chr19.gff")[[1]]
# using bash there are:
# grep 'mRNA' MacaM_Rhesus_Genome_Annotation_v7.6.8_chr19.gff | wc -l
# 1182 genes total, not all are polymorphic
# I got the coordinates of the OXPHOS genes on chr19 like this:
# grep 'mRNA' MacaM_Rhesus_Genome_Annotation_v7.6.8_chr19.gff | egrep 'COX|NDUF|UQCR|ATP5|CYC1|SDHB|SDHA|SDHC|SDHD' | cut -f1,4,5 > coordinates_OXPHOS_mRNAs.txt
# there are 15 OXPHOS genes
# grep 'CDS' MacaM_Rhesus_Genome_Annotation_v7.6.8_chr19.gff | egrep 'COX|NDUF|UQCR|ATP5|CYC1|SDHB|SDHA|SDHC|SDHD' | cut -f1,4,5 | wc -l
# the 15 OXPHOS genes have a total of 59 exons



# Extracts the information stored in the slot region.names
# and converts the strings to numeric values (start position of the
# region and end position)
# this makes a vector with only the start position for all genes
from.pos <- sapply(genes@region.names,function(x)
{return(as.numeric(strsplit(x," ")[[1]][1]))})
# this makes a vector with only the end position for all genes
to.pos <- sapply(genes@region.names,function(x)
{return(as.numeric(strsplit(x," ")[[1]][3]))})
# Lets concatenate the values into a matrix
DATA <- cbind(from.pos, to.pos, genes@MKT)
# now check out the values for the 7 OXPHOS genes on chr01
DATA[from.pos == '15628447']
DATA[from.pos == '37992642'] # NULL
DATA[from.pos == '45636133']
DATA[from.pos == '112030533']
DATA[from.pos == '135484416']
DATA[from.pos == '135601364']
DATA[from.pos == '220779458']
# now check out the values for the 3 OXPHOS genes on chr02a
DATA[from.pos == '37732212']
DATA[from.pos == '42997182']
DATA[from.pos == '92746891']
# now check out the values for the 4 OXPHOS genes on chr02b
DATA[from.pos == '62111164']
DATA[from.pos == '88692897']
DATA[from.pos == '93979633']
DATA[from.pos == '128742129']
# now check out the values for the 3 OXPHOS genes on chr03
DATA[from.pos == '84594254']
DATA[from.pos == '104359835']
DATA[from.pos == '104846075']
DATA[from.pos == '114610662']
DATA[from.pos == '156309713']
DATA[from.pos == '156309713']
DATA[from.pos == '157252481']
# now check out the values for the 6 OXPHOS genes on chr04
DATA[from.pos == '634919']
DATA[from.pos == '8603394']
DATA[from.pos == '46628495']
DATA[from.pos == '61522247']
DATA[from.pos == '97928438']
DATA[from.pos == '139226852']
# now check out the values for the 7 OXPHOS genes on chr05
DATA[from.pos == '191800']
DATA[from.pos == '1601204']
DATA[from.pos == '54356765']
DATA[from.pos == '56278755']
DATA[from.pos == '58728795']
DATA[from.pos == '130350642']
DATA[from.pos == '138178553']
# now check out the values for the 2 OXPHOS genes on chr06
DATA[from.pos == '73250180']
DATA[from.pos == '95211059']
# now check out the values for the 7 OXPHOS genes on chr07
DATA[from.pos == '3770437']
DATA[from.pos == '12722004']
DATA[from.pos == '21199179']
DATA[from.pos == '34825571']
DATA[from.pos == '97408775']
DATA[from.pos == '124582469']
DATA[from.pos == '166649885']
```

# Permutations

One expectation, particularly for PF511 and PM626 is that they should have mt-interacting nuclear genes from the species whose mtDNA they carry.  Most of the introgressed regions are heterozygous, but this is worth testing any how. This really should be compared among windows that carry genes, not across all genes. I will do this as follows: read in the coordinates of the interacting OXPHOS (or any) genes plus all the other genes. And then read in coordinates of heterozygous and homoz introgressed bits.  Then see how many interactors are in the introgressed bits, shuffle the interaction designation many times, and see if it is more (PF511) or less (PM626) than expected by chance.

First step is to generate the introgression file from the admixfrog output.  This has to be done in a tedious way for each chr because the chr name needs to be substituted:

```
xzcat tonk_PF511_chr01_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr01,/g' > tonk_PF511_for_perm_chr01_admix.out
xzcat tonk_PF511_chr02a_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr02a,/g' > tonk_PF511_for_perm_chr02a_admix.out
xzcat tonk_PF511_chr02b_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr02b,/g' > tonk_PF511_for_perm_chr02b_admix.out
xzcat tonk_PF511_chr03_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr03,/g' > tonk_PF511_for_perm_chr03_admix.out
xzcat tonk_PF511_chr04_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr04,/g' > tonk_PF511_for_perm_chr04_admix.out
xzcat tonk_PF511_chr05_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr05,/g' > tonk_PF511_for_perm_chr05_admix.out
xzcat tonk_PF511_chr06_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr06,/g' > tonk_PF511_for_perm_chr06_admix.out
xzcat tonk_PF511_chr07_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr07,/g' > tonk_PF511_for_perm_chr07_admix.out
xzcat tonk_PF511_chr08_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr08,/g' > tonk_PF511_for_perm_chr08_admix.out
xzcat tonk_PF511_chr09_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr09,/g' > tonk_PF511_for_perm_chr09_admix.out
xzcat tonk_PF511_chr10_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr10,/g' > tonk_PF511_for_perm_chr10_admix.out
xzcat tonk_PF511_chr11_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr11,/g' > tonk_PF511_for_perm_chr11_admix.out
xzcat tonk_PF511_chr12_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr12,/g' > tonk_PF511_for_perm_chr12_admix.out
xzcat tonk_PF511_chr13_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr13,/g' > tonk_PF511_for_perm_chr13_admix.out
xzcat tonk_PF511_chr14_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr14,/g' > tonk_PF511_for_perm_chr14_admix.out
xzcat tonk_PF511_chr15_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr15,/g' > tonk_PF511_for_perm_chr15_admix.out
xzcat tonk_PF511_chr16_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr16,/g' > tonk_PF511_for_perm_chr16_admix.out
xzcat tonk_PF511_chr17_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr17,/g' > tonk_PF511_for_perm_chr17_admix.out
xzcat tonk_PF511_chr18_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr18,/g' > tonk_PF511_for_perm_chr18_admix.out
xzcat tonk_PF511_chr19_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chr19,/g' > tonk_PF511_for_perm_chr19_admix.out
xzcat tonk_PF511_chrX_TON_HEC_MAU.out.bin.xz | sed 's/ch,/chrX,/g' > tonk_PF511_for_perm_chrX_admix.out

```

And the same for PF626. And then use the famous "Makes_inputfile_for_jackknife.pl" script to concatenate them.

