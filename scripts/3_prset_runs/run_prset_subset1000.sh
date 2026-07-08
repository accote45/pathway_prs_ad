#BSUB -L /bin/sh
#BSUB -n 48
#BSUB -J prsetexpr_subset1000
#BSUB -R rusage[mem=3072]
#BSUB -R select[mem>=3072]
#BSUB -R span[hosts=1]
#BSUB -q premium             # target queue for job execution
#BSUB -W 24:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o o.prsetexpr_subset1000
#BSUB -e e.prsetexpr_subset1000

# Repeat of the PRSet AD dx association (run_prset.sh) restricted to the
# 1000-case / 1000-control subset. Differences vs run_prset.sh:
#   --keep       -> subset keep file (built by subset_ad_phenotype_1000.R)
#   --pheno      -> subset residual phenotype file
#   --out        -> separate results prefix
# Run subset_ad_phenotype_1000.R first to generate the two subset files.

ml R
Rscript /sc/arion/projects/psychgen/cotea02_prset/PRSice.R --prsice /sc/arion/projects/psychgen/cotea02_prset/PRSice_linux \
    --snp rsid \
    --pvalue p_value \
    --beta \
    --stat beta \
    --fastscore \
    --bar-levels 1 \
    --a1 effect_allele \
    --a2 other_allele \
    --base /sc/arion/projects/paul_oreilly/data/GWASs/NonBiobanks/raw_data/ad/PGC3_Unpublished/xukbb/case_control/xukbb_case_control_eur_neff0.6_nsumstats1.txt.gz \
    --binary-target F \
    --clump-kb 1000kb \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --x-range chr19:44409039-46412650 \
    --extract /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb18177-qc.snplist \
    --gtf /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/reference/Homo_sapiens.GRCh37.75.gtf.gz \
    --keep /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/samples/ad_subset1000_keep.txt \
    --msigdb /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/master_no_apoe.gmt \
    --num-auto 22 \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results/ad_case.control_prset_nothreshold_eur_subset1000 \
    --pheno /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/phenotypes/ad_phenotype_residuals_subset1000.txt \
    --pheno-col AD_resid \
    --set-perm 10000 \
    --target /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb18177_chr1.22 \
    --thread 48 \
    --ultra  \
    --background /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/background_genes.txt:gene \
    --wind-3 35kb \
    --wind-5 35kb