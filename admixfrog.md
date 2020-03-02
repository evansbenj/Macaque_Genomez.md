# Admixfrog

Vcf files are here on graham:
```
/home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/females_and_males/FandM*vcf.gz
```

On graham I installed admixfrog here `/home/ben/.local/bin/admixfrog` as follows:
```
module load scipy-stack/2019b
pip install cython scipy --upgrade --user
pip install git+https://github.com/benjaminpeter/admixfrog --user
```

Now I think it works based on this:
```
/home/ben/.local/bin/admixfrog --help
```
