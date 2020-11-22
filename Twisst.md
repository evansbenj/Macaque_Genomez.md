# Twisst

A nice complement to a sliding window analysis is Twisst, which infers phylogenies in genomic windows.

A first step is to phase the data; I did this with Beagle 5.0 like this:
```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=def-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar beagle.18May20.d20.jar gt=../FandM_${1}_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.g
z.recode.vcf.gz out=../FandM_${1}_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz impute=t
rue 
```

Then I maked a phased geno file:
```
python ./genomics_general/VCF_processing/parseVCF.py -i FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | gzip > phased_genos/chr01.geno.gz
```
In this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst
```
and only for chrX, I had to convert the male genotypes to diploid like this before phasing:
```
echo "chrX 1 148935249 M 2" > ploidy.txt
bcftools +fixploidy ../all_diploid_haploid_chrX_BSQR_filtered3_noPAR_SNPsonly.vcf.gz.recode.vcf.gz.recode.vcf.gz -Ov -- -p ploidy.txt > chrX_diploid.vcf
```


I had to add a line to the 'phyml_sliding_windows.py' script to let it know where to look for the genomics.py file:

```
sys.path.append("/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dept
h_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general")
import genomics
```

and then I made nj trees for each window like this:
```
python genomics_general/phylo/phyml_sliding_windows.py -T 10 -g phased_genos/chr01.geno.gz --prefix phased_genos/output.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n
```
