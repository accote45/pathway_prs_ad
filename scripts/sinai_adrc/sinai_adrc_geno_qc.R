



We then filter the genotypes for other QC based on the maximum genotype missingness (typically < 0.5; ≥ 0.95 call rate), Hardy Weinberg Equilibrium (HWE; default 1e-6, but PLINK recommends 1e-50) and minimum minor allele frequency (MAF; often ≥ 0.01 or ≥ 1%) specified in the configuration file.

Violations of Hardy Weinberg Equilibrium can indicate either the presence of population substructure, or the occurence of genotyping error. It is common practice to assume that violoations are indicative of genotyping error and remove SNPs in which the HWE test statistic is highly significant.

Filtering SNPs on MAF, HWE and call-rate can be done in `PLINK 1.9` by typing the following at the shell prompt:

plink --bfile [raw-GWA-data] \
  --geno [chosen missingness threshold] \
  --maf [chosen MAF threshold] \
  --keep-allele-order \
  --autosome \
  --hardy \
  --hwe 0.000001 \
  --make-bed --out [filtered-GWA-data]


Filtering samples on call-rate (95% here) can be done in `PLINK 1.9` by typing the following at the shell prompt:

```bash
plink --bfile raw-GWA-data  \
  --mind 0.05 \
  --make-bed --out filtered-GWA-data


Identification of individuals with discordent sex can be done in PLINK 1.9 by typing the following at the shell prompt, which will produce a list of individuals with discordent sex data.

```bash
plink --bfile raw-GWA-data  \
  --check-sex --out --out output.sexcheck



To construct a sample of unrelated individules, participants with the highest number of pairwise kinship coefficents > 0.1875 are then iterativly removed.


Individuals with outlying heterozygosity rates can be identified in PLINK 1.9 by typing the following command at the shell prompt:

```
plink --bfile raw-GWA-data \
  --extract snplist.prune.in \
  --het --out output.het
```

This produces a file containing Method-of-moments F coefficient estimates, which can be used to calculate the observed heterozygosity rate in each individual. Analysis is performed using an LD pruned snplist.

We calculate a heterozygocity similarly using observed and expected counts from the PLINK output [(Observed - Expected)/N) and exclude samples that are ± 3 sd from the cohort mean.


The following `PLINK 1.9` commands for performing the PCA analysis can be entered at the shell prompt to generate two files containing the principal component eigenvalues and Principal component eigenvectors.

```bash
plink --bfile merged-reference-sample-data \
  --pca 10 --within sample_population.txt \
  --pca-clusters population_clusters.txt \
  --out pca.output
```

We remove any samples not within 6 standard deviations of the chosen superpopulation on all 10 PCs.



To obtain the principal components for the sample dataset after population outliers have been removed, type the following `PLINK 1.9` commands at the shell prompt to generate the principal component eigenvalue and eigenvector files.

```bash
plink --bfile raw-GWA-data \
  --remove fail-ancestry-QC.txt \
  --pca 10 \
  --out filter-GWA-data
```

However, ancestry can bias the dimensionality reduction, and PC's may represent relatedness as well as gross population structure. To account for this, we have an option in the configuration (on by default) to use a similar algorithm to PC-AiR, from the GENESIS R package.


