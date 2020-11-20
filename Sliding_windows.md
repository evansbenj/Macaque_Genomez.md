Sliding Windows Analyses

Simon Martin has nice software that calculates D and fdm stats in sliding windows.  I also wrote a program to do this but it makes more sense to use his since it is probably better vetted than mine (!).

I am working in this directory on graham:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing
```
Firsts step is to convert my filtered vcf files to geno format like this:
```
python parseVCF.py -i ../../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.recode.vcf.gz -o chr1.geno.gz
```

Then it is necessary to swap any astrisks with Ns:
```
gunzip chr18.geno.gz
sed -i 's/\*/N/g' chr18.geno 
gzip -c chr18.geno > chr18.geno.gz
```

for autosomes:
```
#!/bin/sh
#SBATCH --job-name=abba
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=8gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

# sbatch ABBABABA.sh chr H1 H2 H3 O
# sbatch ABBABABA.sh chr01 nig nge hec papio

# populations
# bru papio hec mau nem sum nig nge tog ton


module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2

echo python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5
}.csv -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFai
ledWindows --windType coordinate

python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5}.csv
 -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFailedWi
ndows --windType coordinate
```

for chrX:
```
#!/bin/sh
#SBATCH --job-name=abba
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

# sbatch ABBABABA.sh chr H1 H2 H3 O
# sbatch ABBABABA_chrX.sh chrX nga nge hec papio

# populations
# bru papio hec mau nem sum nig nge tog ton


module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2

echo python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5
}.csv -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFai
ledWindows --windType coordinate --haploid maura_PM613,maura_PM614,maura_PM616,nem_PM1206,nem_PM664,nem_PM665,nem_Suka
i_male,nigra_PM1003,nigrescens_PM1011,nigrescens_PM654,tonk_PM592

python3 ABBABABAwindows.py -g ./VCF_processing/${1}.geno.gz -f phased -o ./VCF_processing/${1}_${2}_${3}_${4}_${5}.csv
 -w 100000 -m 100 -s 100000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops.txt --writeFailedWi
ndows --windType coordinate --haploid maura_PM613,maura_PM614,maura_PM616,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_mal
e,nigra_PM1003,nigrescens_PM1011,nigrescens_PM654,tonk_PM592
```
