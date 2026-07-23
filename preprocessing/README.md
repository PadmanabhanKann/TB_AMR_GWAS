# preprocessing/

Everything that happens before `workflow/01_annotation`: pulling raw reads,
QC-filtering assemblies, and building phenotype files. Organized by
function rather than by organism, since most of these are generic,
CLI-parameterized tools; organism-specific *invocations* (hardcoded paths,
one-time runs) live in `jobs/{gonorrhoeae,tuberculosis}/` instead — the
same job/logic split used throughout `workflow/`.

```
sra_download/     pull reads from ENA/SRA, fix paired FASTQs
assembly_qc/      QUAST-based assembly quality filtering
phenotype_prep/   MIC → log2 phenotype conversion, resistant/susceptible calls
genome_selection/ QC-threshold filtering, copying selected genomes
```

`download_sra.sh` takes an `-i <biosample_list.txt>` flag so the same
script covers every organism's read list rather than being copy-pasted per
organism.
