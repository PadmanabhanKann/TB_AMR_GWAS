# jobs/

HPC (SLURM) submission scripts — hardcoded cluster paths, `#SBATCH`
headers, one specific run rather than reusable logic. Split by organism;
the reusable logic these jobs call into (where it exists independently)
lives in `../workflow/`, `../preprocessing/`, or `../alt_methods/`.

```
gonorrhoeae/    N. gonorrhoeae job submissions
tuberculosis/   M. tuberculosis job submissions
  ├── standard/     standard (unconditional) PySeer GWAS
  └── conditional/  conditional PySeer GWAS (conditioned on major resistance loci)
shared/         jobs run on generic SRA/assembly paths with no
                organism marker recoverable from the script itself
```

`tuberculosis/conditional/post_pyseer.sh` is the canonical post-processing
script for the GWAS stage — it annotates, summarises, and sorts hits for
both the significant *and* the full k-mer results. An earlier,
narrower version (significant hits only) has been removed in favor of this
one.

All paths inside these scripts point at cluster storage — this folder is
an organized record of what was submitted, not a live working directory.
