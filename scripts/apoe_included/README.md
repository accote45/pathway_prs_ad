# APOE-included PRSet AD-dx analysis

Repeat of the full-sample PRSet AD case/control association
(`scripts/run_prset.sh`) with **APOE retained** throughout, instead of
regressed out of the phenotype and stripped from the pathway file.

## What changes vs the no-APOE analysis

| Component | No-APOE (original) | APOE-included (here) |
|---|---|---|
| Phenotype residualization | covariates include `APOE_genotype` | `APOE_genotype` dropped from covariates |
| Pathway file (`--msigdb`) | `master_no_apoe.gmt` | `master.gmt` (APOE gene retained) |
| Background (`--background`) | `background_genes.txt` | `background_genes_withapoe.txt` |
| chr19 APOE region | `--x-range chr19:44409039-46412650` excluded | flag removed; region clumped & scored |
| Phenotype file | `ad_phenotype_residuals.txt` | `ad_phenotype_residuals_withapoe.txt` |
| Results prefix | `..._eur` | `..._eur_withapoe` |

`master.gmt` (APOE-included pathways) is already produced by
`scripts/create_master_gmt.R`, so no GMT rebuild is needed here.

## Run order

```sh
# 1. Phenotype residuals with APOE NOT regressed out
Rscript create_ad_phenotype_withapoe.R

# 2. APOE-included PRSet background gene list (from master.gmt)
Rscript create_background_withapoe.R

# 3. Full-sample PRSet run
bsub < run_prset_withapoe.sh

# 4. Compare APOE-included vs APOE-excluded pathway results
Rscript plot_prset_withapoe_vs_noapoe.R
```

## Comparison figures (`plot_prset_withapoe_vs_noapoe.R`)

Compares the APOE-included `.summary` against the original APOE-excluded run
(`scripts/run_prset.sh` output) for AD case/control. Unlike the subset-vs-full
comparison, divergence is the *expected* result: including APOE routes the
dominant chr19 signal into every pathway containing an APOE-region gene, so
those pathways gain competitive signal. Writes to `figures/withapoe_vs_noapoe/`:

- `prset_withapoe_vs_noapoe_logP.pdf` — −log10(competitive P) concordance scatter,
  APOE-excluded vs APOE-included, permutation floor drawn, top gainers in colour.
- `prset_withapoe_vs_noapoe_R2.pdf` — same for PRS R².
- `prset_withapoe_vs_noapoe_gainers_dumbbell.pdf` — the pathways that gain the most
  competitive signal when APOE is scored (the APOE-driven sets).
- `prset_withapoe_vs_noapoe_table.csv` — per-pathway P/R2/rank in each run + deltas,
  sorted by `logP_delta`.
- `prset_withapoe_vs_noapoe_stats.txt` — Kendall tau-b / Spearman rho and top-K
  hypergeometric overlap tests.

```sh
# defaults point at the cluster results dir; override for local runs
Rscript plot_prset_withapoe_vs_noapoe.R \
  --noapoe   <ad_case.control_prset_nothreshold_eur.summary> \
  --withapoe <ad_case.control_prset_nothreshold_eur_withapoe.summary> \
  --outdir figures/withapoe_vs_noapoe
```

GO/MP term names come from the `--godict`/`--mgidict` files on Minerva; run there
for readable pathway labels (locally, terms fall back to raw GO/MP IDs).

## New files written to `data/`

- `ad_phenotype_residuals_withapoe.txt`
- `background_genes_withapoe.txt`

All paths in the scripts are the absolute cluster paths under
`/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/`, matching the
originals.
