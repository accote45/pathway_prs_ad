##### AFR ANCESTRY

#BSUB -L /bin/sh
#BSUB -n 48
#BSUB -J prset_afr_bmi_rppr
#BSUB -R rusage[mem=2048]
#BSUB -R select[mem>=2048]
#BSUB -R span[hosts=1]
#BSUB -q premium             # target queue for job execution
#BSUB -W 5:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o o.prset_afr_bmi_rppr
#BSUB -e e.prset_afr_bmi_rppr

ml R
Rscript /sc/arion/projects/psychgen/cotea02_prset/PRSice.R --prsice /sc/arion/projects/psychgen/cotea02_prset/PRSice_linux \
    --snp hg19_RefSNP_id \
    --pvalue pval \
    --beta \
    --stat beta \
    --fastscore \
    --bar-levels 1 \
    --a1 ea \
    --a2 ref \
    --ignore-fid \
    --base /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gwas/MVP_R4.1000G_AGR.BMI_Mean_INT.EUR.GIA.dbGaP.txt_hg19_maf0.05 \
    --binary-target F \
    --clump-kb 1000kb \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --extract /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/bmi_prset_nothreshold_afr.valid \
    --gtf /sc/arion/projects/paul_oreilly/lab/cotea02/project/data/reference/Homo_sapiens.GRCh37.75.gtf.gz \
    --keep /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/afr_sample_ids_80pc.txt \
    --msigdb /sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt \
    --num-auto 22 \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/bmi/bmi_prset_nothreshold_afr \
    --pheno /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forprset.txt \
    --pheno-col bmi_resid \
    --set-perm 1000 \
    --target /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_maf01 \
    --thread 48 \
    --ultra  \
    --background /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/background_genes.txt:gene \
    --wind-3 35kb \
    --wind-5 35kb



    --extract /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_chr_maf01_qc.snps \






##### EUR ANCESTRY

# downsample EUR target sample to match AFR sample 
set.seed(42)  # for reproducibility

# Read sample ID files
afr_ids <- read.table("/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/afr_sample_ids_80pc.txt", header = FALSE)
eur_ids <- read.table("/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/eur_sample_ids_80pc.txt", header = FALSE)
pheno   <- read.table("/sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forprset.txt", header = TRUE)

# IDs with non-missing phenotype
pheno_ids <- pheno$IID[!is.na(pheno$ad_resid)]

# Intersect both ancestry groups with phenotype-available IDs
afr_with_pheno <- afr_ids[afr_ids$V1 %in% pheno_ids, , drop = FALSE]
eur_with_pheno <- eur_ids[eur_ids$V1 %in% pheno_ids, , drop = FALSE]

n_afr <- nrow(afr_with_pheno)
cat("AFR N with phenotype:", n_afr, "\n")
cat("EUR N with phenotype:", nrow(eur_with_pheno), "\n")

# Downsample EUR (phenotype-available) to match AFR N
eur_downsampled <- eur_with_pheno[sample(nrow(eur_with_pheno), n_afr), , drop = FALSE]

write.table(eur_downsampled,
            "/sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/eur_sample_ids_downsampled.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)



#BSUB -L /bin/sh
#BSUB -n 48
#BSUB -J prset_eur_bmi_rppr
#BSUB -R rusage[mem=2048]
#BSUB -R select[mem>=2048]
#BSUB -R span[hosts=1]
#BSUB -q premium             # target queue for job execution
#BSUB -W 5:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o o.prset_eur_bmi_rppr
#BSUB -e e.prset_eur_bmi_rppr

ml R
Rscript /sc/arion/projects/psychgen/cotea02_prset/PRSice.R --prsice /sc/arion/projects/psychgen/cotea02_prset/PRSice_linux \
    --snp hg19_RefSNP_id \
    --pvalue pval \
    --beta \
    --stat beta \
    --fastscore \
    --bar-levels 1 \
    --a1 ea \
    --a2 ref \
    --ignore-fid \
    --base /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/gwas/MVP_R4.1000G_AGR.BMI_Mean_INT.EUR.GIA.dbGaP.txt_hg19_maf0.05 \
    --binary-target F \
    --clump-kb 1000kb \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --extract /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/bmi/bmi_prset_nothreshold_eur.valid \
    --gtf /sc/arion/projects/paul_oreilly/lab/cotea02/project/data/reference/Homo_sapiens.GRCh37.75.gtf.gz \
    --ignore-fid \
    --keep /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/eur_sample_ids_downsampled.txt  \
    --msigdb /sc/arion/projects/psychgen/cotea02_prset/pathway_overlap/data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt \
    --num-auto 22 \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/bmi/bmi_prset_nothreshold_eur \
    --pheno /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forprset.txt \
    --pheno-col bmi_resid \
    --set-perm 1000 \
    --target /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb18177_chr1.22 \
    --thread 48 \
    --ultra  \
    --background /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/background_genes.txt:gene \
    --wind-3 35kb \
    --wind-5 35kb