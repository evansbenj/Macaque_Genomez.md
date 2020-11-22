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
python ./genomics_general/VCF_processing/parseVCF.py -i FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 | gzip > phased_genos/chr01.geno.gz
```
