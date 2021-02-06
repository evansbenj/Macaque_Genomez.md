# Testing for LD between nonsynonymous SNPs in the mtDNA and variation in the nuclear genome

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

module load plink/1.9b_5.2-x86_64
module load StdEnv/2020
module load r/4.0.2
plink --vcf temp.vcf.gz --recode --const-fid 0 --chr-set 36 no-y no-xy no-mt --out myplink
plink --file myplink --pheno sex_phenotype --assoc --allow-no-sex
where the "sex_phenotype" file is a tab-delimited file that looks like this:

0	sample1	1
0	sample2	2
0	sample3	1
0	sample4	2
0	sample5	2
0	sample6	1
The first column is the family ID (just zeros here). The second column is the sample name - this is the same as in the vcf file. The third column is the phenotype - should use 1 and 2, NOT 1 and 0, because 0 might be interpreted as a missing phenotye.

Also this flag "--const-fid 0" sets the family id to zero and tells plink to use the vcf sample name as the sample ID irrespective of whether there is an underscore in the name.

This flag "--chr-set 36" allows extra chrs. They will be numbers in the order they are encountered in the vcf file (I think - this will need to be confirmed...)
