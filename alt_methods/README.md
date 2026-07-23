# alt_methods/

Alternative or exploratory GWAS / population-genetics methods run alongside
the main PySeer + PANGWES pipeline, on *M. tuberculosis* only. These are
**not** part of the eight-stage `workflow/` scheme — kept separate so that
scheme stays exactly the eight core stages.

```
gemma/          GEMMA linear mixed model GWAS (alternative to PySeer)
treewas/        treeWAS tree-based GWAS across all 13 antibiotics
genepi/         GenEpi epistasis detection (rifampicin)
hierbaps/       hierBAPS population structure / clustering
superdca/       SuperDCA coevolution analysis
clonalframeml/  recombination-corrected phylogeny
```
