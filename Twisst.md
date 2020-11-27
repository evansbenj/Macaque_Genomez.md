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

Then I made a phased geno file:
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
and then make a geno file from this...
```
python ../genomics_general/VCF_processing/parseVCF.py -i ../all_diploid_haploid_chrX_phased.vcf.gz.vcf.gz | gzip > ../phased_genos/chrX.geno.gz 
```


I had to add a line to the 'phyml_sliding_windows.py' script to let it know where to look for the genomics.py file:

```
sys.path.append("/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dept
h_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general")
import genomics
```
I also had to install phyml and add it to my path like this:
```
PATH=$PATH:/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/phyml/src
```
and also load several modules in some weird order...:
```
module load python/3.8
module load StdEnv/2020
module load scipy-stack/2020b
```

and then I made nj trees for each window like this:
```
python3 genomics_general/phylo/phyml_sliding_windows.py -T 10 -g phased_genos/chr02a.geno.gz --prefix phased_genos/chr02a_treez_w50 -w 50 --windType sites --model GTR
```

To get twisst to work I had to update the ete3 package first:
```
pip install --upgrade ete3
```

And now this seems to work (after changing the ../genomics_general/pops_twisst.txt file to not inlcude separate sum and tog populations):
```
python twisst.py -t ../phased_genos/chr01_treez_w50.trees.gz -w ../phased_genos/chr01_output.weights.csv.gz --outputTopos ../phased_genos/topologies_chr01.trees --outgroup papio -g nig -g nge -g hec -g ton -g mau -g bru -g nem -g papio --method complete --groupsFile ../genomics_general/pops_twisst.txt
```

And I'm running twisst with this sbatch script for each chr:
```
#!/bin/sh
#SBATCH --job-name=twisst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=128gb
#SBATCH --output=twisst.%J.out
#SBATCH --error=twisst.%J.err
#SBATCH --account=def-ben

# sbatch Twisst.sh chr
module load StdEnv/2020
module load scipy-stack/2020b
# may have to type this before running: pip install --upgrade ete3

python twisst.py -t ../phased_genos/${1}_treez_w50.trees.gz -w ../phased_genos/${1}_twisstoutput.weights.csv.gz --outputTopos 
../phased_genos/${1}_twisst_topologies.trees --outgroup papio -g nig -g nge -g hec -g ton -g mau -g bru -g nem -g papio --meth
od complete --groupsFile ../genomics_general/pops_twisst.txt --verbose
```
