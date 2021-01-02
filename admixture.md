Great directions are here: https://speciationgenomics.github.io/ADMIXTURE/


working directory: 
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data
```

```
mkdir ADMIXTURE
cd ADMIXTURE
```
use bcftools to combine phased SNPs into one file to feed into plink
```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.10.2
bcftools concat -o autosomes.vcf ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz 
```
compress and index
```
bgzip -c autosomes.vcf > autosomes.vcf.gz
tabix -p vcf autosomes.vcf.gz
``
convert to geno format using plink
first make a bed file and remove any SNP with no data for autosomes:

```
module load nixpkgs/16.09  intel/2016.4 plink/1.9b_5.2-x86_64

plink --vcf ./autosomes.vcf.gz --make-bed --geno 0.999 --out ./autosomes --allow-extra-chr --const-fid
```
and then separately for chrX
```
plink --vcf ../all_diploid_haploid_chrX_phased.vcf.gz.vcf.gz --make-bed --geno 0.999 --out ./chrX --allow-extra-chr --const-fid
```
now run admixture for k=2
```
module load StdEnv/2020 nixpkgs/16.09 admixture/1.3.0
admixture --cv chrX.bed 2 > chrXlog2.out
```
or try K=2 to 9
```
for i in {2..9}
do
 admixture --cv autsomes.bed $i > autosomeslog${i}.out
done
```
For plotting, first make a file with individual and species names
```
awk '{split($2,name,"_"); print name[2],name[1]}' chr01.nosex > chr01.list
```
edit to fix papio



download R script to plot admixture 
```
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r
chmod +x plotADMIXTURE.r
```
plot it
```
Rscript plotADMIXTURE.r -p chr01 -i chr01.list -k 2 -l papio,nem,nigra,nigrescens,hecki,tonk,tog,maura,bru 
```
