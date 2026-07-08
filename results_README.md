# PRSet results — pathway PRS for Alzheimer's disease

> Place this file at the top of the results directory as `README.md`.
> Layout produced by `scripts/reorganize_results.sh` (in the code repo).

This directory holds the **PRSet** (competitive pathway polygenic risk score) results
for the AD pathway-PRS project. Every run uses the same AD GWAS base
(`xukbb_case_control_eur_neff0.6_nsumstats1`), the same UK Biobank European target
genotypes, the same pruned pathway set (`master_no_apoe.gmt`), and `--set-perm 10000`
competitive permutations. Phenotypes are pre-residualized in R, so PRSice regresses a
scored pathway PRS against a residual column.

## Folder layout

| Folder | Contents |
|---|---|
| `ad_diagnosis/main/` | AD diagnosis (case/control), primary analysis |
| `ad_diagnosis/withapoe/` | AD diagnosis with APOE **retained** (sensitivity analysis) |
| `ad_diagnosis/subset1000/` | AD diagnosis on a 1000-case/1000-control subset (down-sampling, power analysis) |
| `ad_ageofonset/` | AD age of onset (cases only) |
| `biomarkers/` | Plasma protein endophenotypes: **NEFL**, **GFAP** (UKB Olink) |
| `mri/` | MRI endophenotypes: **Hipp_Mean_Vol** (hippocampal volume), **WMH_Vol** (white-matter hyperintensities) |
| `figures/` | Publication figures (concordance, lollipop, top-pathway table) |
| `qc_diagnostics/` | Phenotype normality diagnostics (`Rplots.pdf`, `hist_*.pdf`, `qq_*.pdf`) |
| `_superseded/` | An earlier NEFL/GFAP run, replaced |
| `bnmf_ad/` | **A SEPARATE ANALYSIS** bNMF analysis — see its own README |

## File naming key

Each run writes a quartet sharing one prefix:

| Extension | What it is | Size |
|---|---|---|
| `.summary` | Per-pathway competitive P, PRS R², coefficients | ~1–3 MB |
| `.prsice` | Per-pathway/threshold PRSice detail | ~1–2 MB |
| `.log` | PRSice run log (exact parameters) | ~5 KB |
| `.best` | Per-individual PRS scores —| ~94 GB |

Prefix tokens: `prset` = PRSet run · `nothreshold` = single P-threshold of 1 (all SNPs,
no tuning) · `eur` = European-ancestry samples · `withapoe` = APOE kept in covariates +
gene set + background and chr19 not excluded (default runs regress APOE out and use
`master_no_apoe.gmt`) · `subset1000` = 1000/1000 down-sample.

## Provenance — which script produced each analysis

Run scripts live in the code repo under `scripts/3_prset_runs/` (LSF `bsub` job files
calling PRSice.R). See `scripts/README.md` for the full pipeline.

| Folder | Produced by |
|---|---|
| `ad_diagnosis/main/` | `scripts/3_prset_runs/run_prset.sh` |
| `ad_diagnosis/withapoe/` | `scripts/apoe_included/run_prset_withapoe.sh` |
| `ad_diagnosis/subset1000/` | `scripts/3_prset_runs/run_prset_subset1000.sh` |
| `ad_ageofonset/` | `scripts/3_prset_runs/run_prset_ageofonset.sh` |
| `biomarkers/` | `scripts/3_prset_runs/run_prset_NEFL.sh`, `scripts/3_prset_runs/run_prset_GFAP.sh` |
| `mri/` | `scripts/3_prset_runs/run_prset_Hipp_Mean_Vol.sh`, `scripts/3_prset_runs/run_prset_WMH_Vol.sh` |

Figures in `figures/` are made by `scripts/4_figures/plot_prset_results.R` (main AD
figures), `scripts/4_figures/plot_prset_biomarkers.R`,
`scripts/4_figures/plot_prset_subset_vs_full.R`, and
`scripts/apoe_included/plot_prset_withapoe_vs_noapoe.R`. All of them read the small
`.summary` files.

> Note: the plot scripts have the *old flat* results paths hardcoded as defaults. After
> this reorganization, point them at the new per-folder `.summary` paths (each script
> takes CLI flags for its inputs).

## The large `.best` files

Every analysis folder contains a ~94 GB `.best` file of per-individual PRS scores. They
are:

- **Not needed for any figure** — the figures use `.summary`. Only
  `scripts/4_figures/count_cases_controls.R` streams `.best` (to tabulate the
  case/control composition of each regression sample).
- **Regenerable** by re-running the corresponding `scripts/3_prset_runs/run_prset*.sh`.
- **The intended final cleanup.** Once you have run `count_cases_controls.R` (or don't
  need the counts), delete them to reclaim ~650 GB:

  ```sh
  find . -name '*.best' -size +1G -print     # review the list first
  find . -name '*.best' -size +1G -delete    # then delete
  ```

## `_superseded/`

Holds an earlier combined NEFL/GFAP PRSet run (dotted names,
`prset_nothreshold_eur.NEFL.best` / `.GFAP.best` and the phenotype-less
`prset_nothreshold_eur.{log,prsice,summary}`). It was replaced by the per-biomarker
underscore runs now in `biomarkers/`, which the plot scripts use. Safe to delete.

## `bnmf_ad/`

A **different project** (its own Snakemake pipeline and README) that is nested here for
convenience. It is not part of the pathway-PRS results and is left untouched — see the
README inside that folder.
