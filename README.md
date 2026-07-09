# Pathway-partitioned PRS for Alzheimer's disease

Pathway-level polygenic risk scoring (**PRSet**) for Alzheimer's disease (AD) in the
UK Biobank. Instead of one genome-wide PRS, we partition AD genetic risk across curated
biological **pathways** (gene sets) and test — with competitive permutation — which
pathways carry risk signal for AD diagnosis and a set of AD-related endophenotypes.

- **Base GWAS:** European-ancestry AD case/control summary statistics.
- **Target:** UK Biobank European-ancestry genotypes.
- **Method:** [PRSet](https://choishingwan.github.io/PRSice/) (PRSice-2), competitive
  scoring against a matched gene background, `--set-perm 10000`.
- **Phenotypes:** AD diagnosis (primary), AD age of onset, plasma biomarkers (NEFL,
  GFAP), and MRI endophenotypes (hippocampal volume, white-matter hyperintensities).
- **APOE handling:** the primary analysis regresses out and excludes *APOE*; a parallel
  sensitivity analysis retains it throughout (`scripts/apoe_included/`).

> **Compute environment.** This is analysis code, not a packaged tool. All paths inside
> the scripts are absolute Minerva HPC paths (`/sc/arion/...`) and the PRSet runs are LSF
> (`bsub`) jobs — they are meant to run on the cluster, not locally.

## Pipeline

Scripts in [`scripts/`](scripts/) are grouped into four stages, in run order:

| Stage | Directory | Does |
|---|---|---|
| 1 | [`1_pathway_prep/`](scripts/1_pathway_prep/) | Build & prune the pathway gene sets → `master.gmt` + PRSet background |
| 2 | [`2_phenotypes/`](scripts/2_phenotypes/) | Extract UKB phenotypes, residualize + rank-inverse-normal transform |
| 3 | [`3_prset_runs/`](scripts/3_prset_runs/) | Run PRSet (one LSF job per phenotype) → `.summary`/`.prsice`/`.log`/`.best` |
| 4 | [`4_figures/`](scripts/4_figures/) | Downstream analysis + publication figures (read the `.summary` files) |

Two side analyses live alongside the main pipeline:

- [`scripts/apoe_included/`](scripts/apoe_included/) — the AD case/control run with *APOE*
  retained (self-contained, has its own README).
- [`scripts/sinai_adrc/`](scripts/sinai_adrc/) — genotype QC/prep for a Mount Sinai ADRC
  replication cohort (work in progress; not part of the UKB pipeline).

## Repository layout

```
pathway_prs_ad/
├── scripts/          all analysis code (see scripts/README.md)
├── data_README.md    describes the (git-ignored) data/ directory
└── results_README.md describes the (git-ignored) results/ directory
```

The bulky `data/` and `results/` directories are not tracked in git; their structure and
contents are documented in [data_README.md](data_README.md) and
[results_README.md](results_README.md).

## Documentation map

- [scripts/README.md](scripts/README.md) — full stage-by-stage pipeline description.
- [data_README.md](data_README.md) — every input / intermediate file and the script that
  makes it.
- [results_README.md](results_README.md) — PRSet output layout, the `.summary`/`.best`
  file quartet, and provenance.
- [scripts/apoe_included/README.md](scripts/apoe_included/README.md) — the APOE-retained
  sensitivity analysis.
</content>
</invoke>
