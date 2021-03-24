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
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_bor_roh.txt
```
for mau
```
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_mau_roh.txt
```
for ton
```
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr01_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' > RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr02a_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr02b_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr03_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr04_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr05_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr06_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr07_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr08_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr09_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr10_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr11_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr12_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr13_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr14_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr15_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr16_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr17_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr18_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 ../FandM_chr19_BSQR_jointgeno_allsites_withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz | grep 'RG' >> RG_ton_roh.txt
```
```
bcftools roh -G30 --AF-dflt 0.4 -s nem_GumGum_female,nem_PM1206,nem_PM664,nem_PM665,nem_Sukai_male -M 100 chr03_phased_named_SNPs.vcf.gz > bor_chr03_roh.txt

bcftools roh -G30 --AF-dflt 0.4 -s maura_PF615,maura_PF713,maura_PM613,maura_PM614,maura_PM616 -M 100 chr03_phased_named_SNPs.vcf.gz > mau_chr03_roh.txt

bcftools roh -G30 --AF-dflt 0.4 -s tonk_PF511,tonk_PF559,tonk_PF563,tonk_PF597,tonk_PF626,tonk_PM592 -M 100 chr03_phased_named_SNPs.vcf.gz > ton_chr03_roh.txt

bcftools roh -G30 --AF-dflt 0.4 -s PF505,hecki_PF643,hecki_PF644,hecki_PF647,hecki_PF648 -M 100 chr03_phased_named_SNPs.vcf.gz > hec_chr03_roh.txt

bcftools roh -G30 --AF-dflt 0.4 -s nigra_PF1001,nigra_PF660,nigra_PM1003 -M 100 chr03_phased_named_SNPs.vcf.gz > nig_chr03_roh.txt
```

```
cat bor_chr03_roh.txt | grep 'RG' > RG_bor_chr03_roh.txt
cat mau_chr03_roh.txt | grep 'RG' > RG_mau_chr03_roh.txt
cat ton_chr03_roh.txt | grep 'RG' > RG_ton_chr03_roh.txt
cat hec_chr03_roh.txt | grep 'RG' > RG_hec_chr03_roh.txt
cat nig_chr03_roh.txt | grep 'RG' > RG_nig_chr03_roh.txt
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


