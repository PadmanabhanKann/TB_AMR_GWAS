# config/

One YAML per organism (`gonorrhoeae.yaml`, `tuberculosis.yaml`) holding the
parameters that differ between the two pipelines: genome counts, loci of
interest, antibiotic panel, and which tool variant is used at each stage
(e.g. FastTree vs. IQ-TREE for phylogeny).

Scripts under `workflow/` are designed to take a `--config` flag pointing
at one of these files rather than hardcoding organism-specific values —
that's the mechanism that keeps pipeline logic shared instead of forked.

`paths.raw_data` and `paths.results` are conventions for where a user's own
input data and pipeline outputs should live — this repository ships code
only, not data.
