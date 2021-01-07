I am going to use an approach called XPCLR to perform a genomewide scan for selection.
working directory:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/XPCLR/xpclr
```

1. Download the code:
```
git clone https://github.com/hardingnj/xpclr.git
```
2. Load python module:
```
module load python
```
3. Create a virtual environment:
```
virtualenv /home/$USER/xpclr_venv
```
4. Activate your virtual environment:
```
source /home/$USER/xpclr_venv/bin/activate
```
5. Install the dependencies:

The package come with a file "requirements.txt" that contains the list of dependencies. Usually, the installation should work by running "pip install -r requirements.txt". However, in this case, numpy is required for installing "scikit-allel". Therefore, "numpy" should appear first on the list. You can just edit the file and change the order of the packages or install numpy first.

```
(xpclr_venv) [~]$ cd xpclr
(xpclr_venv) [xpclr]$ cat requirements.txt
scikit-allel>=1.2
numpy
pandas
scipy
h5py
zarr
```

To install the dependencies, I did run the following:

```
pip install numpy --no-index
pip install -r requirements.txt
```
6. Install the package "xpclr":
```
(xpclr_venv) [~]$ cd xpclr
(xpclr_venv) [~]$ python setup.py install
(xpclr_venv) [~]$ deactivate
```
7. Test the installation
```
[~]$ module load python
[~]$ source /home/$USER/xpclr_venv/bin/activate
(xpclr_venv) [~]$ xpclr --help
```
This seems to work for doing the analysis on all chrs (still haven't gotten an output yet though):
```
#!/bin/sh
#SBATCH --job-name=xpclr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=16gb
#SBATCH --output=xpclr.%J.out
#SBATCH --error=xpclr.%J.err
#SBATCH --account=def-ben

# sbatch xpclr.sh pop1 pop2 
module load nixpkgs/16.09 intel/2016.4 tabix/0.2.6
module load python
source /home/$USER/xpclr_venv/bin/activate


# chr 01, 2a 2b
xpclr --format vcf -Sa ${1} -Sb ${2} -I ../../FandM_chr01_BSQR_jointgeno_allsites_with
papio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz --rrate 0.448e-8 -C chr0
1 --size 100000 --step 100000 -O ${1}_vs_${2}_chr01_xpclr.out

xpclr --format vcf -Sa ${1} -Sb ${2} -I ../../FandM_chr02a_BSQR_jointgeno_allsites_wit
hpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz --rrate 0.448e-8 -C chr
02a --size 100000 --step 100000 -O ${1}_vs_${2}_chr0$2a_xpclr.out

xpclr --format vcf -Sa ${1} -Sb ${2} -I ../../FandM_chr02b_BSQR_jointgeno_allsites_wit
hpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz --rrate 0.448e-8 -C chr
02b --size 100000 --step 100000 -O ${1}_vs_${2}_chr02b_xpclr.out

# # chr 03 to 09
for i in {3..9}
do
 xpclr --format vcf -Sa ${1} -Sb ${2} -I ../../FandM_chr0${i}_BSQR_jointgeno_allsites_
withpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz --rrate 0.448e-8 -C 
chr0${i} --size 100000 --step 100000 -O ${1}_vs_${2}_chr0${i}_xpclr.out
done

for i in {10..19}
do
 xpclr --format vcf -Sa ${1} -Sb ${2} -I ../../FandM_chr${i}_BSQR_jointgeno_allsites_w
ithpapio_filtered2_coverage_SNPsonly.vcf.gz.phased.vcf.gz.vcf.gz --rrate 0.448e-8 -C c
hr${i} --size 100000 --step 100000 -O ${1}_vs_${2}_chr${i}_xpclr.out
done


xpclr --format vcf -Sa ${1} -Sb ${2} -I ../../all_diploid_haploid_chrX_phased.vcf.gz.v
cf.gz --rrate 0.448e-8 -C chrX --size 100000 --step 100000 -O ${1}_vs_${2}_chrX_xpclr.
out
```
