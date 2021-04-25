FOr ROH, I am working in this directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/ROH
```

# bcftools
```
module load bcftools
```

Batch process all chrs for each population.  Here is bor:
```
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
```
for mau
```
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
```
for ton
```
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
```
fir hec
```
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_hec_roh.txt
```
for nigra
```
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_nig_roh.txt
```

Now isolate the ROHs for each population
```
cat bor_chr03_roh.txt | grep 'RG' > RG_bor_chr03_roh.txt
cat mau_chr03_roh.txt | grep 'RG' > RG_mau_chr03_roh.txt
cat ton_chr03_roh.txt | grep 'RG' > RG_ton_chr03_roh.txt
cat hec_chr03_roh.txt | grep 'RG' > RG_hec_chr03_roh.txt
cat nig_chr03_roh.txt | grep 'RG' > RG_nig_chr03_roh.txt
```

Here is a script to do permutations and calculate statistics from these:
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# N_interact genes and their acronyms, and all other genes

# The other file is a file with ROH blocks from a population with lengths

# First identify how which ROH blocks have N_interact genes and which do not.
# test stat is the difference in the summed length of ROH blocks that do and do not have N_interact genes
# then randomly assign N_interact genes and 

# run like this:
# ./POH_permutation.pl FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez.txt RG_bor_roh.txt

my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];

my @windowsites;
my @Fst_values;
my $sumsites=0;
my @temp;
my $y;
my $x;
my %N_interact_hash;
my @interact_perm;

# first open up the N_interact_hash gene info 
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'gene'){ 
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[0];
		$N_interact_hash{$temp[2]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[5];
		push(@interact_perm,$temp[5]); # this will be used for the permutations later
	}	
}		
close DATAINPUT;

# now open up the ROH block data
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

my @temp1;
my $n_ROH_blocks_with_interacting_genez=0;
my $n_ROH_blocks_with_other_genez=0;
my $n_ROH_blocks_without_genez=0;

my $length_sum_ROH_blocks_with_interacting_genez=0;
my $length_sum_ROH_blocks_with_other_genez=0;
my $length_sum_ROH_blocks_without_genez=0;

my $n_Ninteract_genes_on_ROH_blockz=0;
my $n_genes_on_ROH_blockz=0;
my $n_genes_on_non_ROH_blockz=0;
my %ROH_blocks;
my %perm_blockz; 
my @length_array_for_perms;
my @N_array_for_perms;

while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split('\t',$line);
	# ignore first line
	if($temp[0] ne '# RG'){
			# intially set the block to be no genes (0)
			$ROH_blocks{$temp[2]."_".$temp[3]."_".$temp[4]}{"genes"} = 0; 
			# also assume that it does not have any interacting genes
			$ROH_blocks{$temp[2]."_".$temp[3]."_".$temp[4]}{"interacting"} = 0;
			# also keep track of lengths
			$ROH_blocks{$temp[2]."_".$temp[3]."_".$temp[4]}{"length"} = $temp[5];
			# Now cycle through all the genes to see if any are in this block
			foreach my $key (keys %N_interact_hash){
				@temp1=split('_',$key);
				# now check if this block contains any genes
				if(
					($temp1[0] eq $temp[2])&&($temp1[1] >= $temp[3])&&($temp1[1] <= $temp[4])
					# beginning of the gene is in this block
					|
					($temp1[0] eq $temp[2])&&($temp1[2] >= $temp[3])&&($temp1[2] <= $temp[4])
					# end is in this block	
					){
						# this ROH block has a gene
						$ROH_blocks{$temp[2]."_".$temp[3]."_".$temp[4]}{"genes"} += 1;
						$n_genes_on_ROH_blockz+=1;
						push(@length_array_for_perms,$temp[5]);
						# check if it is an interacting gene
						if($N_interact_hash{$key}{"mt_interact"} == 1){
							$ROH_blocks{$temp[2]."_".$temp[3]."_".$temp[4]}{"interacting"} = 1;
							print $N_interact_hash{$key}{"gene"}," is in ROH block ",$temp[2]."_".$temp[3]."_".$temp[4],"\n";
							$n_Ninteract_genes_on_ROH_blockz+=1;
						}
				}
			}
	}	
}

# ok now I have a hash that has information on whether or not each ROH block has any genes
# and whether or not any of these genes have any N_interact genes.
# print these numbers

foreach my $key (keys %ROH_blocks){
	if($ROH_blocks{$key}{"genes"} > 0){
		if($ROH_blocks{$key}{"interacting"} == 1){
			$n_ROH_blocks_with_interacting_genez+=1;
			$length_sum_ROH_blocks_with_interacting_genez+=$ROH_blocks{$key}{"length"};
		}
		else{
			$n_ROH_blocks_with_other_genez+=1;
			$length_sum_ROH_blocks_with_other_genez+=$ROH_blocks{$key}{"length"};
		}	
	}
	else{
		$n_ROH_blocks_without_genez+=1;
		$length_sum_ROH_blocks_without_genez+=$ROH_blocks{$key}{"length"};
	}	
}	

print "Number of ROH blocks with N_interact genes: ",$n_ROH_blocks_with_interacting_genez,"\n";
print "Number of ROH blocks with other genes: ",$n_ROH_blocks_with_other_genez,"\n";
print "Number of ROH blocks without genes: ",$n_ROH_blocks_without_genez,"\n";

print "Average length of ROH blocks with N_interact genes: ",$length_sum_ROH_blocks_with_interacting_genez/$n_ROH_blocks_with_interacting_genez,"\n";
print "Average length of ROH blocks with other genes: ",$length_sum_ROH_blocks_with_other_genez/$n_ROH_blocks_with_other_genez,"\n";
print "Average length of ROH blocks without genes: ",$length_sum_ROH_blocks_without_genez/$n_ROH_blocks_without_genez,"\n";

print "Number of genes on ROH blocks: ",$n_genes_on_ROH_blockz,"\n";
print "Proportion of ROH blocks with genes that have N_interact genes ",$n_ROH_blocks_with_interacting_genez/
($n_ROH_blocks_with_interacting_genez+$n_ROH_blocks_with_other_genez),"\n";
print "Proportion of genes on ROH blocks that are N_interact genes ",$n_Ninteract_genes_on_ROH_blockz/$n_genes_on_ROH_blockz,"\n";

# This is the test statistic
my $test_stat2 = printf("%.1f",$length_sum_ROH_blocks_with_interacting_genez/$n_ROH_blocks_with_interacting_genez-
$length_sum_ROH_blocks_with_other_genez/$n_ROH_blocks_with_other_genez);

######################
# Permutations
######################

my $perms=1000;
my $switch=0;
my $counter=0;
my $introg_interacter=0;
my $introg_withgenez=0;
my $introg_without_genez=0;
my $introg_withgenez_switch=0;
my $introg_interacter_switch=0;
my $genes_on_ROH_blocks=0;
my $Ninteract_genes_on_ROH_blocks=0;
my $n_ROH_blocks_with_interacting_genez_perms=0;
my $n_ROH_blocks_with_other_genez_perms=0;

my @quick_perm_genez;

for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@interact_perm );    # permutes the N_interact assignment for each of 16045 genes
		$counter=0;
		$switch=0;
		$length_sum_ROH_blocks_with_interacting_genez=0;
		$length_sum_ROH_blocks_with_other_genez=0;
		$length_sum_ROH_blocks_without_genez=0;
		foreach my $key (keys %ROH_blocks){ # go through all the blocks
			$switch=0; # assume no genes in block
			if($ROH_blocks{$key}{"genes"} > 0){ # does this block have any genes?
				for ($x = 0 ; $x < $ROH_blocks{$key}{"genes"}; $x++ ) {
					if($interact_perm[$counter] == 1){
						$switch=1; # if we find any 1s, this means there is an N_interact gene
					}
					$counter+=1;
				}
				if($switch == 1){
					$length_sum_ROH_blocks_with_interacting_genez+=$ROH_blocks{$key}{"length"};
					$n_ROH_blocks_with_interacting_genez_perms+=1;
				}	
				else{
					$length_sum_ROH_blocks_with_other_genez+=$ROH_blocks{$key}{"length"};
					$n_ROH_blocks_with_other_genez_perms+=1;
				}
			}
			else{ # block has no genes
				$length_sum_ROH_blocks_without_genez+=$ROH_blocks{$key}{"length"};
			}	
		}
		push(@quick_perm_genez,printf("%.1f",(($length_sum_ROH_blocks_with_interacting_genez/$n_ROH_blocks_with_interacting_genez_perms)-
			($length_sum_ROH_blocks_with_other_genez/$n_ROH_blocks_with_other_genez_perms))));	
}	


if($#quick_perm_genez != $perms-1){
	print "Hey, something wrong with perms\n";
}

my @quick_perm_genez_sorted = sort { $a <=> $b } @quick_perm_genez;
$switch=0;
my $pval=$perms; # this will make the pval be zero if the test stat is larger than all the perms
$counter=0;
print "@quick_perm_genez_sorted\n";
# now figure out where the test stat is
for ($y = 0 ; $y <= $#quick_perm_genez_sorted; $y++ ) {
	if(($test_stat2 <= $quick_perm_genez_sorted[$y])&&($switch==0)){ 
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	


print "Test stat:",$test_stat2,"\n";
print "QuickP = ",1-$pval/$perms,"\n";



# fisher_yates_shuffle( \@array ) : 
    # generate a random permutation of @array in place
    sub fisher_yates_shuffle {
        my $array = shift;
        my $i;
        for ($i = @$array; --$i; ) {
            my $j = int rand ($i+1);
            next if $i == $j;
            @$array[$i,$j] = @$array[$j,$i];
        }
    }

```

# Rehh  

I am going to try to use the R package rehh to detect long runs of homozygosity.

First step is to find some SNPs in each N_interact gene and name them in a vcf file.


First I made a bed file (`Ninteract_chr03.bed`)that has the ranges of all N_interact genes on chr03. I got this information from the gff file (more conveniently summarized in my table). Paste it in this file:
```
chr03	84594254	84613076	NDUFB5
chr03	104846075	104847927	NDUFAF3
chr03	148084830	148114936	ACAD9
chr03	156309713	156315708	NDUFB4
chr03	157411293	157436295	TIMMDC1
```
Then extract the SNPs like this:
```
module load tabix
tabix FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz -B Ninteract_chr03.bed > Ninteract_chr03_SNPs.vcf
```
Here are the first SNPs in each gene
```
chr03	84594254	84613076	NDUFB5 84594255
chr03	104846075	104847927	NDUFAF3  104846083
chr03	148084830	148114936	ACAD9  148085019
chr03	156309713	156315708	NDUFB4 156309746
chr03	157411293	157436295	TIMMDC1 157411323
```
Some (all?) of these are not variable in each species and the first one isn't even variable in the taxa at all (it just differs from the ref).  But this should be ok.

Now I am going to edit a version of the phased VCF file to add names to these SNPs (I'll call them each by the gene acronym).

Then I compressed it:
```
bgzip -c chr03_phased_named_SNPs.vcf > chr03_phased_named_SNPs.vcf.gz
```

Running on compute canada - issues with sh scrip were resolved with flags:
```
sbatch --account rrg-ben --time=72:00:00 ./ROH_permutation.pl FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez_orientation.txt RG_ton_roh.txt ROH_ton_perms.out ROH_ton_density.out
```

