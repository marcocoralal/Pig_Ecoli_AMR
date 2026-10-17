# Porcine *E. coli* antimicrobial non-susceptibility — reproducible analysis

This repository reproduces the antimicrobial non-susceptibility statistical analysis in 344 porcine *E. coli* isolates. It contains **no original producer names, addresses, cantons, sampling dates, or sample IDs**.

## Run

1. Install R (>= 4.2) and packages `dplyr`, `car`, `logistf`, `geepack`, `igraph`.
2. Open R in this repository root (the directory containing `R/`, `data/`, and `README.md`).
3. Run `source("R/00_master_analysis.R")`.
4. Verify `Analysis completed; validation checks passed.`. Generated tables and Figure 1 appear in `results/`.

The script runs **directly from the included public CSV**. The confidential Excel workbook is neither needed nor distributed.

## Methods

Twelve outcome-specific multivariable logistic models (11 predictors each), 132 Wald coefficient tests, Benjamini–Hochberg FDR correction, VIF, McFadden pseudo-R²; Firth penalized logistic regression (profile likelihood) and week-clustered GEE (exchangeable correlation) sensitivity analyses; network of GLM FDR-significant conditional associations. Directed arrows mean predictor-to-outcome model specification **only**, not biological causation.

The expected validation totals are 344 complete records, 132 tests, 40 nominal associations, 25 GLM FDR-significant associations (22 positive, 3 negative), 19 GEE FDR-significant associations, 17 GLM/GEE overlaps, 25 GLM/Firth overlaps, and 17 significant across all three methods. The study author also obtained these totals by running the public-CSV script in a fresh R session, with the final message `Analysis completed; validation checks passed.`.

For definitions of all columns, phenotype coding, and data provenance, see [DATA_DICTIONARY.md](DATA_DICTIONARY.md).\n\n## Variables and confidentiality

`data/AMR_deidentified_analysis.csv` contains 344 rows and 13 columns: `week_cluster` plus 12 `ns_` phenotypes. Week labels are recoded consecutive cluster identifiers, not dates. `ns_` replaces misleading source `sens_` headers. Interpretation is **non-susceptibility** (intermediate or resistant = 1; susceptible = 0), **pending independent confirmation from laboratory coding documentation**.

The file has been minimized to remove direct identifiers, but rare profiles and week clustering mean it is **de-identified, not guaranteed anonymous**. Public release remains subject to data-sharing authorization and disclosure-risk assessment. Do not upload the original workbook or any private data. No license is granted until the study owners choose one.

## Outputs

`results/` is populated by the script with GLM, Firth and GEE tables, model diagnostics, Table 4 CSVs, Figure 1 network PDF, and validation results. No Word or private Excel files are needed.

## Caveat

The script is a reproducibility companion for the analyses above, not a replacement for checking sampling denominators, clinical breakpoint interpretation, or raw laboratory provenance. A directed network is a visualization of outcome-specific model coefficients, not a causal network.
