# Data dictionary

## Dataset
`data/AMR_deidentified_analysis.csv` is the analysis-ready dataset for the antimicrobial non-susceptibility association analyses of porcine *Escherichia coli* isolates collected at a slaughter facility in Quito, Ecuador.

- **Unit of observation:** one isolate per row.
- **Records:** 344 complete observations used in the multivariable analyses.
- **Columns:** 13 (one sampling-week cluster identifier and 12 binary antimicrobial phenotypes).
- **Missing values:** none in the published analysis matrix.
- **Phenotype coding:** `0 = susceptible (S)`; `1 = non-susceptible (intermediate or resistant, I + R)`. **This interpretation must be checked against the original laboratory coding documentation before public release.**
- **Important:** `ns_` means *non-susceptibility*, not resistance alone.

| Column | Meaning |
|---|---|
| `week_cluster` | Recoded sampling-week cluster (1–25); used for GEE clustering. Values are not calendar dates. |
| `ns_sxt` | Trimethoprim–sulfamethoxazole non-susceptibility |
| `ns_gm` | Gentamicin non-susceptibility |
| `ns_cip` | Ciprofloxacin non-susceptibility |
| `ns_caz` | Ceftazidime non-susceptibility |
| `ns_tet` | Tetracycline non-susceptibility |
| `ns_s` | Streptomycin non-susceptibility |
| `ns_chl` | Chloramphenicol non-susceptibility |
| `ns_cf` | Cefoxitin non-susceptibility (**CF = cefoxitin**, not cefalotin) |
| `ns_ak` | Amikacin non-susceptibility |
| `ns_azm` | Azithromycin non-susceptibility |
| `ns_ff` | Florfenicol non-susceptibility |
| `ns_amc` | Amoxicillin–clavulanic acid non-susceptibility |

## Provenance and preprocessing

The source was a laboratory/field workbook containing antimicrobial susceptibility test interpretations and sampling metadata. The public analysis matrix retains only the twelve binary phenotypes and a recoded sampling-week identifier. Names, addresses, original sample identifiers, original dates, and geographic origin fields were not retained.

The 344 observations correspond to records with complete phenotype data used by the models. The dataset **does not independently document the complete original sampling, isolation, selection, or exclusion flow**; consult the manuscript and original study records for those details. The 344-row analysis matrix must not be interpreted as the number of all animals sampled or all isolates screened.

## Reproducibility and limitations

Run `source("R/00_master_analysis.R")` from the repository root using R. The included dataset supports GLM, Firth penalized regression, and sampling-week-clustered GEE, plus the related model diagnostics and network output. Model results describe statistical associations, not causality.

The dataset is **de-identified, not guaranteed anonymous**. Uncommon phenotype combinations and sampling-week membership remain. Publication requires authorization from the relevant data owners and an appropriate disclosure-risk review. The source workbook must not be uploaded.
