# Alpha Bost-Connes Data and Scripts

This repository contains reproducibility materials for the manuscript:

**A Bost-Connes Arithmetic Hypothesis for the Fine-Structure Constant over
Q(sqrt(5))**

Author: Jian Zhou, Principal Investigator, JZ Institute of Science,
Hong Kong, China. Email: jack@jzis.org.

## Purpose

The repository is intended to support a Frontiers in Physics submission by
making the numerical scripts, processed outputs, and source inventories public.
It does not contain copyrighted paper PDFs, arXiv source bundles, or large raw
archives. Those files remain local and are represented here only by metadata
manifests.

## Contents

```text
scripts/
  frontiers_quadratic_field_screen.py
  frontiers_grammar_audit.py

data/processed/
  frontiers_quadratic_field_screen.csv
  frontiers_grammar_audit_summary.csv
  frontiers_grammar_audit_top20.csv

data/metadata/
  post2011_source_inventory.csv
  curated_source_inventory.csv
  reference_papers_manifest.csv
  post2011_direct_data_manifest.csv
  post2011_hepdata_exact_download_attempts.csv
  post2011_arxiv_source_table_candidates.csv
  local_hepdata_summary.csv

figures/
  alpha_inverse_residuals.png
  alpha_inverse_residuals_300dpi.jpg
  pdg_compas_r_ratio_landscape.png

docs/
  Frontiers_Contribution_To_The_Field_2026-04-29.md
```

## Reproduce the quadratic-field negative-control screen

The negative-control screen is dependency-free and uses Python's standard
library. By default it evaluates squarefree real quadratic fields through
`d <= 120`.

```bash
python3 scripts/frontiers_quadratic_field_screen.py > data/processed/frontiers_quadratic_field_screen.csv
```

The script compares the leading arithmetic scale

```text
A_d = (3! / zeta_{Q(sqrt(d))}(-3)) * epsilon_K^(-2)
```

and the same fixed three-term rule

```text
B_d = A_d - w_K * epsilon_K^(-3)
      + (Tr_{K/Q}(epsilon_K^2) * epsilon_K)^(-Delta_K)
```

over small squarefree real quadratic fields, reporting residuals against the
CODATA 2022 central value for alpha(0)^(-1). This is a diagnostic screen only;
it is not a full statistical proof of the manuscript hypothesis.

## Reproduce the finite grammar audit

The grammar audit enumerates a bounded expression space and records the closest
expressions to the CODATA 2022 central value.

```bash
python3 scripts/frontiers_grammar_audit.py \
  --summary-out data/processed/frontiers_grammar_audit_summary.csv \
  --top-out data/processed/frontiers_grammar_audit_top20.csv
```

The default audit uses squarefree `d <= 120`, exponent range `1..8`, torsion
coefficients from `-3` to `3`, and local signs `-1, 0, +1`. It also reports a
fixed-seed Monte Carlo sample from the same finite grammar. These fractions are
descriptive diagnostics, not formal p-values.

## Data policy

The CSV manifests record source coverage and download status for the associated
fine-structure constant, alpha(MZ), and g-2 literature survey. Some HEPData
downloads were blocked by HTTP 403 during the local curation run; this status is
recorded in `post2011_hepdata_exact_download_attempts.csv`.

## Licenses

Code is released under the MIT License. Metadata and generated CSV outputs are
released under CC BY 4.0 where the author has rights to do so. Third-party
paper metadata, experiment names, and bibliographic facts remain subject to
their original sources.
