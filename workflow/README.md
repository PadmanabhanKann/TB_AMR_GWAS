# workflow/

The eight pipeline stages shared by both organisms — nothing else lives
here. Alternative/exploratory methods are in `../alt_methods/`; pre-pipeline
data prep is in `../preprocessing/`.

```
01_annotation/    Prokka
02_pangenome/     Panaroo
03_phylogeny/     FastTree (N. gonorrhoeae) / IQ-TREE (M. tuberculosis)
04_kmers/         unitig-caller
05_gwas/
  ├── prep/            kinship/similarity matrix, references.txt, sample-name sanitizing
  ├── conditional/     conditional-GWAS covariate/k-mer prep + subsetting (M. tuberculosis)
  └── postprocessing/  hit annotation, PA-matrix → Phandango prep, k-mer filtering, variant calling
06_ld/            LD calculation (PLINK) + LD decay plotting (M. tuberculosis)
07_pangwes/       Cuttlefish / SpydrPick / ARACNE support scripts
08_plotting/      Phandango, iTOL, GWES / bubble plots, gene-hit distribution plots
```

Files here are written to be reusable across organisms — generic CLI flags
or a `--config <organism>.yaml` (see `../config/`) rather than hardcoded
paths. The actual per-organism runs (hardcoded HPC paths, SBATCH headers)
live in `../jobs/{gonorrhoeae,tuberculosis}/`.

There is no `05_gwas/standard/` logic script — for the standard (non-
conditional) PySeer run, the job scripts in `jobs/` are the only copy of
that logic (see `jobs/*/pyseer*.sh`).

Three files are vendored from the [PySeer](https://github.com/mgalardini/pyseer)
project (`count_patterns.py`, `phylogeny_distance.py`, `summarise_annotations.py`,
© Marco Galardini & John Lees) rather than original code — copyright
headers are preserved as-is.
