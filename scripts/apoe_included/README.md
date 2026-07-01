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
```

## New files written to `data/`

- `ad_phenotype_residuals_withapoe.txt`
- `background_genes_withapoe.txt`

All paths in the scripts are the absolute cluster paths under
`/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/`, matching the
originals.
