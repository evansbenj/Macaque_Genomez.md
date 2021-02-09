# Testing for LD between nonsynonymous SNPs in the mtDNA and variation in the nuclear genome

Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data
```

First I added a Papio anubis mtDNA genome (GenBank accessionm NC 020006) to the alignment and separated the mtDNA genomes by the 13 protein coding genes.  I exported the protein coding sequences, opened them each in paup and excluded gapped, constant and uninf (parsimony non-informative) positions.  The last category are autapomorphic sites.

Then I further filtered these positions as follows in order to remove variation that is solely attributable to geographic structure:
removed autapomorphic variation
removed variation fixed within species
removed gapped and constant positions
removed variation that was only present within only one species
removed variation that was present only in papio and one macaque individual, population, or species
removed variation that was present in only tonk + tog, even if not fixed in tonk (IBD)
kept variation that was present in some heck plus tonk 511
removed variation where papio was divergent and remaining polymorphism was in only one macaque population or species
removed variation fixed on Sulawesi, fixed in borneo, or fixed in nem, including if this variation was shared with papio.

In other words, only variation that was present in more than one macaque species was considered


Plink
Associations between SNPs and a phenotype (such as phenotypic sex) are easily calculated using plink.
First make plink files out of the vcf files:
```
module load nixpkgs/16.09
module load plink/1.9b_5.2-x86_64
plink --vcf FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz --recode --const-fid 0 --set-missing-var-ids @:# --allow-extra-chr --out chr19_plink
```
where this flag `--set-missing-var-ids @:# ` tells plink to set the SNP ids as "chr:bp".  This was important because the permutation test reports significant SNPs based on SNPID, which was blank without this flag.  The '--allow-extra-chr' flag was needed for chr02a and chr02b.

or, for chrX:
```
plink --vcf all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode.vcf.gz --recode --const-fid 0 --allow-extra-chr --out chrX_plink
```
Now test for associations for each SNP for each chromosome
```
plink --file chr01_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr01_ND1_
plink --file chr02a_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr02a_ND1_
plink --file chr02b_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr02b_ND1_
plink --file chr03_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr03_ND1_
plink --file chr04_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr04_ND1_
plink --file chr05_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr05_ND1_
plink --file chr06_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr06_ND1_
plink --file chr07_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr07_ND1_
plink --file chr08_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr08_ND1_
plink --file chr09_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr09_ND1_
plink --file chr10_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr10_ND1_
plink --file chr11_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr11_ND1_
plink --file chr12_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr12_ND1_
plink --file chr13_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr13_ND1_
plink --file chr14_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr14_ND1_
plink --file chr15_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr15_ND1_
plink --file chr16_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr16_ND1_
plink --file chr17_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr17_ND1_
plink --file chr18_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr18_ND1_
plink --file chr19_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chr19_ND1_
plink --file chrX_plink --pheno ./ND1/ND1_base.txt --assoc --allow-no-sex --allow-extra-chr --all-pheno --pfilter 1e-4 --out ./ND1/chrX_ND1_

```
where the "ND1_base.txt" file is a tab-delimited file that looks like this:
```
0	bru_PF707	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	download	2	1	2	1	1	0	0	1	1	1	1	2	2	2	1	1	1	2
0	hecki_PF505	1	1	2	1	1	1	1	1	1	1	1	1	0	1	0	1	2	1
0	hecki_PF643	1	1	1	1	1	1	2	2	1	1	2	1	2	1	2	1	1	1
0	hecki_PF644	1	1	1	1	1	1	2	2	1	1	2	1	2	1	2	1	1	1
0	hecki_PF647	1	1	2	1	1	1	1	1	1	1	1	1	0	1	0	1	2	1
0	hecki_PF648	1	1	2	1	1	1	1	1	1	1	1	1	0	1	0	1	2	1
0	maura_PF615	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	maura_PF713	2	1	1	1	1	1	1	1	1	2	1	2	1	2	1	1	1	1
0	maura_PM613	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	maura_PM614	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	maura_PM616	2	1	1	1	1	1	1	1	1	1	1	2	1	2	1	1	1	1
0	nem_GumGum_female	1	2	1	2	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nem_Ngsang_sumatra_female	1	1	1	0	2	1	1	1	1	2	1	2	2	2	1	2	2	1
0	nem_PM1206	1	1	2	1	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nem_PM664	1	1	1	1	1	2	1	1	2	1	1	1	1	1	1	1	1	2
0	nem_PM665	1	1	2	1	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nem_Sukai_male	1	2	1	2	1	2	1	1	2	1	1	1	1	1	1	1	1	1
0	nigra_PF1001	1	1	1	1	1	1	1	1	1	1	1	1	2	0	1	1	1	1
0	nigra_PF660	1	1	1	1	1	1	1	1	1	1	1	1	2	0	1	1	2	1
0	nigra_PM1003	1	1	1	1	1	1	1	1	1	1	1	1	2	0	1	1	1	1
0	nigrescens_PM1011	1	1	1	1	1	1	1	1	1	1	1	1	1	0	1	1	1	2
0	nigrescens_PM654	1	1	1	1	1	1	1	1	1	1	1	1	1	0	1	1	1	2
0	tog_PF549	1	1	1	2	2	2	1	1	1	1	1	1	1	1	1	2	1	1
0	tonk_PF511	1	1	1	1	1	1	2	2	1	1	2	1	2	1	2	1	1	1
0	tonk_PF559	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	tonk_PF563	1	1	1	1	2	1	1	1	1	1	1	1	1	1	2	1	1	1
0	tonk_PF597	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	tonk_PF626	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	tonk_PM592	2	2	1	1	1	1	1	1	2	1	1	1	1	1	1	1	1	1
```
The first column is the family ID (just zeros here). The second column is the sample name - this is the same as in the vcf file. The third column is the phenotype - should use 1 and 2, NOT 1 and 0, because 0 might be interpreted as a missing phenotye.  In positions where I have three amino acids, I encoded the 3rd one with a zero if it was population specific.  In some cases I had to split the variation up when three amino acids were present but none seemed to be based on population structure (e.g. ATP8 amino acid 9).

Also this flag "--const-fid 0" sets the family id to zero and tells plink to use the vcf sample name as the sample ID irrespective of whether there is an underscore in the name.

This flag "--chr-set 36" allows extra chrs. They will be numbers in the order they are encountered in the vcf file (I think - this will need to be confirmed...)

This --all-pheno tells plink to iterate over all phenotypes

This flag --pfilter 1e-4 tells plink to only output SNPs with Pvals <= 0.0001 (otherwise the files are huge!)
