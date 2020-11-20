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
