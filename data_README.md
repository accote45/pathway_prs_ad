# Data — pathway PRS for Alzheimer's disease

> Place this file at the top of the data directory as `README.md`.
> Layout matches the paths hardcoded in the scripts (see `scripts/README.md`).

Inputs and intermediate files for the AD pathway-PRS (PRSet) pipeline. Files are
grouped by role. The pipeline stages that produce/consume them are in
`scripts/README.md`.

```
data/
├── samples/         ancestry keep-lists + subset keep-file
├── pathways/        pathway gene sets, PRSet background, pruning intermediates
├── phenotypes/      residualized --pheno files for PRSet
├── raw_extracts/    raw UKB phenotype extracts (pre-residualization)
├── qc_diagnostics/  normality / QC plots
├── reference/       external reference inputs copied into the project
└── gwas_for_bnmf/   SEPARATE PROJECT — leave in place
```

## Folder contents

### `samples/`
| File | What it is | Produced by |
|---|---|---|
| `eur_sample_ids_80pc.txt` | European-ancestry IDs (ADMIXTURE ≥80%) | `awk` on the ancestry assignment file (see `create_ad_phenotype.R` header) |
| `afr_sample_ids_80pc.txt` | African-ancestry IDs (ADMIXTURE ≥80%) | same |
| `ad_subset1000_keep.txt` | FID/IID keep-list for the 1000/1000 subset | `subset_ad_phenotype_1000.R` |

### `pathways/`
| File | What it is | Produced by |
|---|---|---|
| `master.gmt` | pooled pathway gene sets (APOE retained) | `create_master_gmt.R` |
| `master_no_apoe.gmt` | same, APOE gene dropped (default runs) | `create_master_gmt.R` |
| `background_genes.txt` | PRSet competitive `--background` (APOE absent) | `create_master_gmt.R` |
| `background_genes_withapoe.txt` | background with APOE retained | `create_background_withapoe.R` |
| `pruned_pathway_output/` | redundancy-pruning outputs (retained gmt/tsv, maps) | `prune_pathways.R` |

### `phenotypes/` (PRSet `--pheno` files; residualized + RINT where noted)
| File | Phenotype | Produced by |
|---|---|---|
| `ad_phenotype_residuals.txt` | AD diagnosis (`AD_resid`) | `create_ad_phenotype.R` |
| `ad_phenotype_residuals_withapoe.txt` | AD diagnosis, APOE retained | `create_ad_phenotype_withapoe.R` |
| `ad_phenotype_residuals_subset1000.txt` | AD diagnosis, 1000/1000 subset | `subset_ad_phenotype_1000.R` |
| `ad_ageofonset_phenotype_residuals.txt` | AD age of onset (`age_of_onset_resid`) | `create_ageofonset_phenotype.R` |
| `residualized_NPX_withNA_RINT.txt` | plasma NEFL/GFAP (RINT) | `create_protein_phenotypes.R`¹ |
| `mri_phenotypes_resid_RINT.txt` | Hipp_Mean_Vol / WMH_Vol (RINT) | `create_mri_phenotypes.R` |

¹ `create_protein_phenotypes.R` writes to the UKB protein-QC directory; copy the
result here into `phenotypes/`.

### `raw_extracts/`
| File | What it is | Produced by |
|---|---|---|
| `UKB_AD_MRI_FSL.csv` | raw FSL MRI IDPs + covariates | `extract_ukb_hipp_wmh.sql` |

### `qc_diagnostics/`
Phenotype normality plots: `hist_GFAP.pdf`, `hist_NEFL.pdf`, `qq_GFAP.pdf`,
`qq_NEFL.pdf` (`ukb_protein_explore.R`), `mri_residual_qqplots.pdf`
(`check_mri_normality.R`).

### `reference/`
External reference inputs copied into the project so it is self-contained:
| File | Source | Used by |
|---|---|---|
| `celltype_specificity.gmt` | `psychgen/…/judit_revisions/…/intermediate_files/` | `create_master_gmt.R` |
| `Homo_sapiens.GRCh37.75.gtf.gz` | `paul_oreilly/lab/cotea02/project/data/reference/` | all `run_prset*.sh` (`--gtf`) |

Other external inputs (GO/MGI dicts, the source pathway `.gmt`s, `Ensembl.regions`,
the GWAS base, and UKB genotype/phenotype data) are read in place from
`/sc/arion/projects/paul_oreilly/data/…` and are **not** copied here.

### `gwas_for_bnmf/`
GWAS summary statistics for a **separate project** (bNMF); not part of this pipeline.
Left in place.
