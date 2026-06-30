
Brian:
-	Has APOE genotype, global ancestry inference, a few PRS, raw and TOPmed imputed genotypes
-	If we need IBD relatedness and QC data
-	If there are phenotypes but no genetic data  brian may have the microarray
o	Brian has pipelines for GWAS and imputed genetic


Please find the data here: /sc/arion/projects/load/data-int/MSSM_ADRC/queries/adrc_2026-01-06
And the sample pipeline here: /sc/arion/projects/load/etc/for_alanna



#ADRC_samples_pass.tsv
    # includes sample ID, PRS, APOE status, ancestry proportions, global ancestry assignment, batch, center, well





library(data.table)
library(tidyverse)

sinai <- read.table('/sc/arion/projects/load/data-int/MSSM_ADRC/queries/adrc_2026-01-06/ADRC_samples_pass.tsv', header = TRUE, sep = '\t')
nacc <- read.csv('/sc/arion/projects/paul_oreilly/data/nacc_adrc/mri_nacc71.csv')
phc <- read.csv('/sc/arion/projects/paul_oreilly/data/nacc_adrc/ADSP-PHC-122024-investigator/Neuropath/NACC_ADSP_PHC_Neuropath_2024.csv')



