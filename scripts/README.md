# Scripts — pathway PRS for Alzheimer's disease

Analysis code for the AD pathway-partitioned PRS (**PRSet**) project. Scripts are
grouped into the four pipeline stages, in run order. All paths inside the scripts are
absolute Minerva cluster paths (`/sc/arion/...`); they are meant to run on the HPC, not
locally.

```
scripts/
├── 1_pathway_prep/   build the pathway gene sets (GMT files)
├── 2_phenotypes/     pull UKB phenotypes + residualize them
├── 3_prset_runs/     PRSet LSF job scripts (produce the results/)
├── 4_figures/        downstream analysis + figures (read the .summary files)
├── apoe_included/    parallel APOE-retained sensitivity analysis (self-contained, own README)
└── sinai_adrc/       separate Mount Sinai ADRC replication cohort (WIP)
```

## Pipeline (run order)

### 1 · `1_pathway_prep/` — build the gene sets
Each `prep_*.R` converts one raw pathway/coexpression database into a standardized GMT
(gene sets in Ensembl IDs): `prep_syngo.R`, `prep_mitocarta.R`, `prep_dsigdb.R`,
`prep_pharmgkb.r`, `prep_proteins_phg.R`, `prep_proteins_johnson2022.R`,
`prep_snRNAseq_dlpfc.R`. Then:
- `prune_pathways.R` — removes redundant gene sets by gene-overlap (containment/Jaccard).
  Its full spec is `instructions_filter_pathways.rmd`.
- `create_master_gmt.R` — pools the pruned sets + cell-type-specificity into `master.gmt`,
  and writes `master_no_apoe.gmt` (APOE gene dropped) + `background_genes.txt` (the PRSet
  competitive background).

### 2 · `2_phenotypes/` — build the phenotypes
SQL extracts of UKB imaging phenotypes, then residualize + rank-inverse-normal transform:
- `extract_ukb_hipp_wmh.sql` → feeds `create_mri_phenotypes.R` (hippocampal vol, WMH vol).
- `extract_ukb_ad_mri.sql` — fuller FreeSurfer extract (alternative).
- `create_ad_phenotype.R` — AD case/control residuals (+ EUR/AFR ancestry keep-lists).
- `create_ageofonset_phenotype.R` — AD age-of-onset residuals (cases only).
- `create_protein_phenotypes.R` — plasma NEFL/GFAP (Olink) residualized + RINT.
- `create_mri_phenotypes.R` — Hipp_Mean_Vol / WMH_Vol residualized + RINT.
- `subset_ad_phenotype_1000.R` — 1000/1000 down-sample for the replication run.
- `check_mri_normality.R` — diagnostic (skew/kurtosis/Q-Q); `ukb_protein_explore.R` — EDA.

### 3 · `3_prset_runs/` — run PRSet
One LSF (`bsub`) job per analysis; all call PRSice.R with `--set-perm 10000`, the xukbb
EUR AD GWAS base, and `master_no_apoe.gmt`. Each writes a `.best/.log/.prsice/.summary`
quartet to `results/`.

| Script | Phenotype |
|---|---|
| `run_prset.sh` | AD diagnosis (main) |
| `run_prset_ageofonset.sh` | AD age of onset |
| `run_prset_NEFL.sh`, `run_prset_GFAP.sh` | plasma biomarkers |
| `run_prset_Hipp_Mean_Vol.sh`, `run_prset_WMH_Vol.sh` | MRI endophenotypes |
| `run_prset_subset1000.sh` | AD dx, 1000/1000 replication |

### 4 · `4_figures/` — analysis & figures
All read the small `.summary` files (drop `Base`/`Background`, handle the permutation
floor at 1/(perm+1)):
- `plot_prset_results.R` — main AD figures (lollipop, cross-phenotype concordance, table).
- `plot_prset_biomarkers.R` — pathway × biomarker/MRI panel.
- `plot_prset_subset_vs_full.R` — full vs 1000/1000 replication concordance.
- `count_cases_controls.R` — the only figure script that streams the large `.best` files.

## Side analyses
- **`apoe_included/`** — repeats the AD case/control run with APOE retained throughout.
  Self-contained with its own `README.md`, phenotype, background, run, and comparison scripts.
- **`sinai_adrc/`** — genotype QC + prep for a Mount Sinai ADRC replication cohort (WIP);
  not part of the UKB PRSet pipeline. `paths` is a scratch note of ADRC data locations.
