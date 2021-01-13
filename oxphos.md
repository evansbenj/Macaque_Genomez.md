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

Here's a preliminary list of all the autosomally encoded proteins in in the OXPHOS complex based on this paper: Genome Res. 2018. 28: 952-967:
```
chr01	sim4cc	mRNA	220779458	220788687	.	+	.	ID=COX20_transcript_01;GID=116228
chr02a	sim4cc	mRNA	42997182	43010254	.	-	.	ID=COX7A2L_transcript_01;GID=9167
chr02a	sim4cc	mRNA	92746891	92749017	.	+	.	ID=COX5B_transcript_01;GID=1329
chr03	sim4cc	mRNA	114610662	114642316	.	-	.	ID=ACOX2_transcript_01;GID=8309
chr03	sim4cc	mRNA	157252481	157260548	.	+	.	ID=COX17_transcript_01;GID=10063
chr04	sim4cc	mRNA	8603394	8653763	.	-	.	ID=ACOX3_transcript_01;GID=8310
chr04	sim4cc	mRNA	46628495	46628917	.	-	.	ID=COX7B2_transcript_01;GID=170712
chr04	sim4cc	mRNA	61522247	61536172	.	+	.	ID=COX18_transcript_01;GID=285521
chr05	sim4cc	mRNA	56278755	56279160	.	+	.	ID=COX7C_transcript_01;GID=1350
chr06	sim4cc	mRNA	73250180	73256003	.	-	.	ID=COX7A2_transcript_01;GID=1347
chr07	sim4cc	mRNA	34825571	34833611	.	-	.	ID=COX19_transcript_01;GID=90639
chr08	sim4cc	mRNA	98367530	98383414	.	-	.	ID=COX6C_transcript_01;GID=1345
chr10	sim4cc	mRNA	95151638	95177391	.	-	.	ID=COX15_transcript_01;GID=1355
chr10	sim4cc	mRNA	95151638	95177391	.	-	.	ID=COX15_transcript_02;GID=1355
chr11	sim4cc	mRNA	8907900	8910660	.	+	.	ID=COX8A_transcript_01;GID=1351
chr12	sim4cc	mRNA	48300943	48310011	.	+	.	ID=COX14_transcript_01;GID=84987
chr12	sim4cc	mRNA	120032059	120034809	.	+	.	ID=COX6A1_transcript_01;GID=1337
chr14	Cufflinks	mRNA	52057142	52077913	.	-	.	ID=COX5A_transcript_01;GID=9377
chr14	sim4cc	mRNA	132329422	132413355	.	-	.	ID=SYNJ2BP-COX16_transcript_01;GID=100529257
chr14	sim4cc	mRNA	132330617	132354632	.	-	.	ID=COX16_transcript_01;GID=51241
chr15	sim4cc	mRNA	32433903	32443159	.	-	.	ID=COX4I2_transcript_01;GID=84701
chr16	sim4cc	mRNA	29233433	29234114	.	-	.	ID=COX6A2_transcript_01;GID=1339
chr16	sim4cc	mRNA	70395661	70403728	.	+	.	ID=COX4I1_transcript_01;GID=1327
chr17	sim4cc	mRNA	14089633	14230045	.	+	.	ID=COX10_transcript_01;GID=1352
chr17	sim4cc	mRNA	38276744	38283435	.	+	.	ID=COX11_transcript_01;GID=1353
chr17	sim4cc	mRNA	69386150	69427976	.	-	.	ID=ACOX1_transcript_01;GID=51
chr19	sim4cc	mRNA	31207304	31214924	.	+	.	ID=COX6B1_transcript_01;GID=1340
chr19	sim4cc	mRNA	31714988	31716937	.	-	.	ID=COX7A1_transcript_01;GID=1346
chr19	sim4cc	mRNA	50026620	50032055	.	-	.	ID=COX6B2_transcript_01;GID=125965
chrX	sim4cc	mRNA	71633448	71639804	.	+	.	ID=COX7B_transcript_01;GID=1349
chr01	sim4cc	mRNA	37992642	37997542	.	+	.	ID=NDUFS5_transcript_01;GID=4725
chr01	sim4cc	mRNA	135484416	135497548	.	+	.	ID=NDUFS2_transcript_01;GID=4720
chr02a	sim4cc	mRNA	37732212	37749061	.	+	.	ID=NDUFAF7_transcript_01;GID=55471
chr02b	sim4cc	mRNA	88692897	88705887	.	+	.	ID=NDUFB3_transcript_01;GID=4709
chr02b	sim4cc	mRNA	93979633	94015898	.	-	.	ID=NDUFS1_transcript_01;GID=4719
chr02b	sim4cc	mRNA	128742129	128807570	.	-	.	ID=NDUFA10_transcript_01;GID=4705
chr03	sim4cc	mRNA	84594254	84613076	.	+	.	ID=NDUFB5_transcript_01;GID=4711
chr03	sim4cc	mRNA	104846075	104847927	.	+	.	ID=NDUFAF3_transcript_01;GID=25915
chr03	sim4cc	mRNA	156309713	156315708	.	-	.	ID=NDUFB4_transcript_01;GID=4710
chr03	sim4cc	mRNA	156309713	156315708	.	-	.	ID=NDUFB4_transcript_02;GID=4710
chr04	sim4cc	mRNA	139226852	139232639	.	-	.	ID=NDUFC1_transcript_01;GID=4717
chr05	sim4cc	mRNA	1601204	1616216	.	+	.	ID=NDUFS6_transcript_01;GID=4726
chr05	sim4cc	mRNA	54356765	54473046	.	+	.	ID=NDUFS4_transcript_01;GID=4724
chr05	sim4cc	mRNA	58728795	58927343	.	+	.	ID=NDUFAF2_transcript_01;GID=91942
chr05	sim4cc	mRNA	138178553	138180781	.	-	.	ID=NDUFA2_transcript_01;GID=4695
chr06	sim4cc	mRNA	95211059	95218436	.	-	.	ID=NDUFAF4_transcript_01;GID=29078
chr07	sim4cc	mRNA	3770437	3785385	.	-	.	ID=NDUFV3_transcript_01;GID=4731
chr07	sim4cc	mRNA	97408775	97416979	.	+	.	ID=NDUFA4_transcript_01;GID=4697
chr07	Cufflinks	mRNA	166649885	166659561	.	+	.	ID=NDUFB2_transcript_01;GID=4708
chr08	sim4cc	mRNA	93359100	93394552	.	+	.	ID=NDUFAF6_transcript_01;GID=137682
chr08	sim4cc	mRNA	123481402	123491824	.	+	.	ID=NDUFB9_transcript_01;GID=4715
chr09	sim4cc	mRNA	59806087	59825964	.	-	.	ID=NDUFB6_transcript_01;GID=4712
chr09	sim4cc	mRNA	59806087	59825964	.	-	.	ID=NDUFB6_transcript_02;GID=4712
chr09	sim4cc	mRNA	94025131	94040826	.	-	.	ID=NDUFA8_transcript_01;GID=4702
chr10	sim4cc	mRNA	95956681	95962600	.	-	.	ID=NDUFB8_transcript_01;GID=4714
chr11	sim4cc	mRNA	6532924	6539484	.	-	.	ID=NDUFS8_transcript_01;GID=4728
chr11	sim4cc	mRNA	6690185	6696071	.	-	.	ID=NDUFV1_transcript_01;GID=4723
chr11	sim4cc	mRNA	18260513	18265960	.	-	.	ID=NDUFS3_transcript_01;GID=4722
chr11	sim4cc	mRNA	69474342	69546903	.	-	.	ID=NDUFC2-KCTD14_transcript_01;GID=100532726
chr11	sim4cc	mRNA	69534557	69546903	.	-	.	ID=NDUFC2_transcript_03;GID=4718
chr11	sim4cc	mRNA	69534557	69546903	.	-	.	ID=NDUFC2_transcript_01;GID=4718
chr11	sim4cc	mRNA	69534557	69546903	.	-	.	ID=NDUFC2_transcript_02;GID=4718
chr12	sim4cc	mRNA	4801409	4840545	.	+	.	ID=NDUFA9_transcript_01;GID=4704
chr12	sim4cc	mRNA	55497421	55499859	.	-	.	ID=NDUFA4L2_transcript_01;GID=56901
chr12	sim4cc	mRNA	94150466	94193941	.	-	.	ID=NDUFA12_transcript_02;GID=55967
chr12	sim4cc	mRNA	94150482	94193941	.	-	.	ID=NDUFA12_transcript_01;GID=55967
chr14	sim4cc	mRNA	17717208	17731422	.	-	.	ID=NDUFAF1_transcript_01;GID=51103
chr15	sim4cc	mRNA	46796983	46833911	.	+	.	ID=NDUFAF5_transcript_01;GID=79133
chr15	sim4cc	mRNA	84086730	84096073	.	-	.	ID=NDUFA6_transcript_01;GID=4700
chr16	sim4cc	mRNA	1930096	1932930	.	+	.	ID=NDUFB10_transcript_01;GID=4716
chr16	sim4cc	mRNA	22236971	22249112	.	-	.	ID=NDUFAB1_transcript_01;GID=4706
chr19	sim4cc	mRNA	1153942	1165810	.	+	.	ID=NDUFS7_transcript_01;GID=374291
chr19	sim4cc	mRNA	5934304	5942627	.	-	.	ID=NDUFA11_transcript_01;GID=126328
chr19	sim4cc	mRNA	8302280	8312129	.	-	.	ID=NDUFA7_transcript_01;GID=4701
chr19	sim4cc	mRNA	9634493	9638249	.	-	.	ID=NDUFA3_transcript_01;GID=4696
chr19	sim4cc	mRNA	14340049	14345977	.	-	.	ID=NDUFB7_transcript_01;GID=4713
chr19	sim4cc	mRNA	19317399	19330902	.	+	.	ID=NDUFA13_transcript_01;GID=51079
chrX	sim4cc	mRNA	47194137	47196604	.	-	.	ID=NDUFB11_transcript_01;GID=54539
chrX	sim4cc	mRNA	113497709	113502918	.	+	.	ID=NDUFA1_transcript_01;GID=4694
chr01	sim4cc	mRNA	45636133	45647963	.	+	.	ID=UQCRH_transcript_01;GID=7388
chr03	sim4cc	mRNA	104359835	104371673	.	-	.	ID=UQCRC1_transcript_01;GID=7384
chr04	blast	mRNA	97928438	97928963	.	-	.	ID=UQCRHL_transcript_01;GID=440567
chr05	sim4cc	mRNA	130350642	130352620	.	+	.	ID=UQCRQ_transcript_01;GID=27089
chr08	sim4cc	mRNA	94674995	94684605	.	-	.	ID=UQCRB_transcript_02;GID=7381
chr08	sim4cc	mRNA	94674995	94684605	.	-	.	ID=UQCRB_transcript_01;GID=7381
chr15	sim4cc	mRNA	71526832	71530391	.	+	.	ID=UQCR10_transcript_02;GID=29796
chr15	sim4cc	mRNA	71526832	71530391	.	+	.	ID=UQCR10_transcript_01;GID=29796
chr16	sim4cc	mRNA	20500148	20530512	.	+	.	ID=UQCRC2_transcript_01;GID=7385
chr19	sim4cc	mRNA	1385443	1392486	.	-	.	ID=UQCR11_transcript_01;GID=10975
chr01	sim4cc	mRNA	112030533	112043548	.	+	.	ID=ATP5F1_transcript_01;GID=515
chr02b	sim4cc	mRNA	62111164	62116182	.	-	.	ID=ATP5G3_transcript_01;GID=518
chr04	sim4cc	mRNA	634919	636858	.	-	.	ID=ATP5I_transcript_01;GID=521
chr07	sim4cc	mRNA	12722004	12734438	.	+	.	ID=ATP5O_transcript_01;GID=539
chr07	sim4cc	mRNA	21199179	21209971	.	+	.	ID=ATP5J_transcript_01;GID=522
chr07	sim4cc	mRNA	124582469	124590730	.	-	.	ID=ATP5J2_transcript_01;GID=9551
chr10	sim4cc	mRNA	7723482	7740495	.	+	.	ID=ATP5C1_transcript_01;GID=509
chr11	sim4cc	mRNA	110459798	110468428	.	+	.	ID=ATP5L_transcript_01;GID=10632
chr12	sim4cc	mRNA	51828101	51839397	.	-	.	ID=ATP5G2_transcript_01;GID=517
chr12	sim4cc	mRNA	54910245	54920096	.	-	.	ID=ATP5B_transcript_01;GID=506
chr14	sim4cc	mRNA	112154408	112167663	.	+	.	ID=ATP5S_transcript_01;GID=27109
chr15	sim4cc	mRNA	5385192	5389178	.	+	.	ID=ATP5E_transcript_01;GID=514
chr17	sim4cc	mRNA	44211923	44214912	.	-	.	ID=ATP5G1_transcript_01;GID=516
chr17	sim4cc	mRNA	68443001	68447706	.	-	.	ID=ATP5H_transcript_01;GID=10476
chr18	sim4cc	mRNA	35035298	35047335	.	-	.	ID=ATP5A1_transcript_01;GID=498
chr19	sim4cc	mRNA	1002804	1005971	.	+	.	ID=ATP5D_transcript_01;GID=513
chr19	sim4cc	mRNA	36878373	36887017	.	-	.	ID=ATP5SL_transcript_03;GID=55101
chr19	sim4cc	mRNA	36878373	36887017	.	-	.	ID=ATP5SL_transcript_02;GID=55101
chr19	sim4cc	mRNA	36878373	36886666	.	-	.	ID=ATP5SL_transcript_01;GID=55101
chr08	sim4cc	mRNA	143057509	143059999	.	+	.	ID=CYC1_transcript_01;GID=1537
chr01	sim4cc	mRNA	15628447	15662313	.	-	.	ID=SDHB_transcript_01;GID=6390
chr05	sim4cc	mRNA	191800	225299	.	+	.	ID=SDHA_transcript_01;GID=6389
chr11	sim4cc	mRNA	12844602	12863522	.	+	.	ID=SDHAF2_transcript_01;GID=54949
chr19	sim4cc	mRNA	31551177	31552175	.	+	.	ID=SDHAF1_transcript_01;GID=644096
chr01	sim4cc	mRNA	135601364	135647524	.	+	.	ID=SDHC_transcript_01;GID=6391
chr11	sim4cc	mRNA	104210163	104219988	.	+	.	ID=SDHD_transcript_01;GID=6392
```
