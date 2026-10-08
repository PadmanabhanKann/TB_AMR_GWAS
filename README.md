# Bacterial GWAS for Antimicrobial Resistance

Analysis code for an integrated genome-wide association study of
antimicrobial resistance across two bacterial pathogens:

1. **_Neisseria gonorrhoeae_** — azithromycin resistance, 1,029 genomes.
   Prokka → Panaroo → FastTree → unitig-caller → PySeer (standard +
   conditional on 23S rRNA) → PANGWES (Cuttlefish → SpydrPick → ARACNE) →
   Phandango / iTOL.
2. **_Mycobacterium tuberculosis_ Lineage 4** — 13 antibiotics, 3,080
   CRyPTIC isolates. Prokka → Panaroo → IQ-TREE → unitig-caller → PySeer
   (standard + conditional) → PLINK LD → PANGWES → plotting / heatmaps.

Both organisms run through the same eight logical pipeline stages
(`workflow/01`–`08`). Where the underlying logic is genuinely
organism-agnostic, it's written once and parameterized (CLI flags, a
`--config` file) rather than duplicated; organism-specific runs — hardcoded
paths, SLURM headers, one particular execution — live under `jobs/`.

## Highlight: automated multi-antibiotic PySeer GWAS

The *M. tuberculosis* arm of this project runs a full PySeer GWAS — run and
post-processed — across all 13 antibiotics without a single per-antibiotic
edit:

- **[`jobs/tuberculosis/standard/pyseer_array.sh`](jobs/tuberculosis/standard/pyseer_array.sh)**
  is a SLURM array job that discovers every antibiotic from the folder
  structure (`find $BASE -mindepth 2 -maxdepth 2 -name "*.tsv"`) and runs
  `pyseer --lmm` against each one as a separate array task.
- **[`jobs/tuberculosis/conditional/post_pyseer.sh`](jobs/tuberculosis/conditional/post_pyseer.sh)**
  sweeps the same antibiotic folders afterward: Bonferroni threshold,
  significant k-mer filtering, Phandango-ready plot files, and gene-level
  annotation + summarisation — for both the significant hits and the full
  result set.

Add a 14th antibiotic and both scripts pick it up automatically; nothing
in either script names an antibiotic directly. This pair is the core
deliverable of the TB pipeline.

**Prerequisites, exact input files, and the expected directory layout are
documented in [`workflow/05_gwas/README.md`](workflow/05_gwas/README.md)**
— read that before trying to run either script.

## Repository structure

```
config/          per-organism parameters (loci, antibiotics, tool choices)
preprocessing/   everything before 01_annotation: SRA download, assembly QC, phenotype prep
workflow/        the 8 shared pipeline stages, organism-agnostic logic
alt_methods/     alternative/exploratory GWAS methods (M. tuberculosis)
jobs/            HPC (SLURM) submission scripts, split by organism
```

Each folder has its own README describing what belongs there.

## Environment

All analyses were run on an HPC cluster via SLURM (`module load` +
`conda`/`source activate`), with paths hardcoded to that cluster's
filesystem in the job scripts under `jobs/`. This repository is the
organized record of that pipeline, not a locally runnable working
directory — reproducing a run means adapting the paths and environment
modules in the relevant `jobs/` script to your own cluster.

Core dependencies across the pipeline: Prokka, Panaroo, FastTree / IQ-TREE,
[unitig-caller](https://github.com/bacpop/unitig-caller),
[PySeer](https://pyseer.readthedocs.io/), PLINK,
[PANGWES](https://github.com/santeripuranen/PANGWES) (Cuttlefish +
SpydrPick), Python 3 (pandas, BioPython, pysam, pyfaidx, gffutils), and R
(ggplot2, ggrepel, treeWAS, ape, data.table).

## Future Implementation
an end to end pipeline for pyseer 

