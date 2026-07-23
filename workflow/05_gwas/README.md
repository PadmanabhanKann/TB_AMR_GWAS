# 05_gwas/

PySeer GWAS, standard and conditional, automated across all 13 antibiotics
for *M. tuberculosis* in one array job plus one post-processing sweep. This
is the core of the pipeline — see below for what it takes to reproduce it.

```
prep/            kinship/similarity matrix, references.txt, sample-name sanitizing
conditional/     conditional-GWAS covariate/k-mer prep + subsetting (tuberculosis)
postprocessing/  hit annotation, significance filtering, PA-matrix → Phandango, variant calling
```

There's no `standard/` folder here — for the standard (non-conditional)
run, the job scripts in `jobs/tuberculosis/standard/` are the only copy of
that logic.

## The automated multi-antibiotic run

Rather than run PySeer once per antibiotic by hand, the TB pipeline drives
all 13 through the same two scripts:

- **`jobs/tuberculosis/standard/pyseer_array.sh`** — a SLURM array job. It
  discovers every antibiotic automatically (`find $BASE -mindepth 2
  -maxdepth 2 -name "*.tsv"`), so each array task picks up one antibiotic's
  phenotype file by index and runs `pyseer --lmm` against it. Adding a 14th
  antibiotic means dropping in one more phenotype file — no script changes.
- **`jobs/tuberculosis/conditional/post_pyseer.sh`** — loops over the same
  antibiotic folders afterward: computes the Bonferroni threshold, filters
  significant k-mers, builds Phandango-ready plot files, and annotates +
  summarises hits to the gene level — for both the significant k-mers *and*
  the full result set. One script handles all 13 antibiotics with no
  per-antibiotic branching.

## Prerequisites

**Software** (one conda env covers both scripts):
- [PySeer](https://pyseer.readthedocs.io/) — provides the `pyseer` command
  plus the `phandango_mapper` and `annotate_hits_pyseer` console scripts
  that `post_pyseer.sh` calls directly.
- Python 3 with `pandas`.
- The vendored PySeer helper scripts already in this repo:
  `05_gwas/postprocessing/count_patterns.py` and `summarise_annotations.py`,
  and `08_plotting/phandango_remapper.py`.
- Standard Unix tools (`awk`, `sort`) — no extra install.

**Inputs required before `pyseer_array.sh`:**

| File | Produced by |
|---|---|
| `unitig.pyseer.gz` — unitig presence/absence matrix | `workflow/04_kmers/unitig-caller.sh` |
| `phylogeny_K.tsv` — kinship/similarity matrix | `workflow/05_gwas/prep/phylogeny_distance.py <tree> --lmm`, run on the `workflow/03_phylogeny` IQ-TREE output |
| One phenotype `.tsv` per antibiotic, in its own subfolder | `preprocessing/phenotype_prep/` (`MIC_to_log.py`, `resistance_phenotype_file.sh`) |

**Additional input required before `post_pyseer.sh`:**

| File | Produced by |
|---|---|
| `references.txt` | `workflow/05_gwas/prep/make_ref.sh` — needs the H37Rv reference FASTA + GFF, plus the draft assemblies and their Panaroo-combined GFFs |

**Expected directory layout** (`$BASE` in both scripts):

```
BASE/
├── references.txt
├── unitig.pyseer.gz
├── phylogeny_K.tsv
├── rifampicin/
│   └── rifampicin.tsv
├── isoniazid/
│   └── isoniazid.tsv
├── ...                       (13 antibiotic folders total)
```

After `pyseer_array.sh` runs, each antibiotic folder gains
`<ab>_kmers_results.txt` and `<ab>_kmers_patterns.txt`; `post_pyseer.sh`
then fills in `significant_kmers.txt`, `annotated_kmers.txt` /
`annotated_allkmers.txt`, `gene_hits.txt` / `gene_all_hits.txt` (+
`.sorted.txt` variants), and a `plot/` folder per antibiotic.

To point this at a new project: change `BASE` in both scripts (or better,
factor it out into a `--config <organism>.yaml` flag — see the note in the
root README about parameterizing rather than forking). Everything else is
already antibiotic-agnostic.
