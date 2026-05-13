

ml plink

awk '{print $1" "$1}' /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/afr_sample_ids_80pc.txt > /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/afr_sample_ids_80pc_plink.txt

# subset AFR samples
    plink \
    --bfile /sc/arion/projects/paul_oreilly/data/Biobanks/UKB/genotyped/ukb18177 \
    --keep /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/afr_sample_ids_80pc_plink.txt \
    --make-bed \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb18177_afr

# calculate MAF in AFR subset
  plink \
    --bfile /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb18177_afr \
    --freq \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_chr


# find SNPs with MAF > 1% in AFR subset
  awk 'NR > 1 && $5 > 0.01 {print $2}' /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_chr.frq > /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_chr_maf01.snps

# filter existing QC SNP list for SNPs with MAF > 1% in AFR subset

# create AFR-only PLINK files with QC'd variant list
for chr in {1..22}; do
  plink \
    --bfile /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb18177_afr \
    --extract /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_chr${chr}_maf01.snps \
    --make-bed \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_chr${chr}_maf01
done