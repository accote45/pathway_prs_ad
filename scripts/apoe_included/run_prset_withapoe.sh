#BSUB -L /bin/sh
#BSUB -n 48
#BSUB -J prsetexpr_nothreshold_withapoe
#BSUB -R rusage[mem=3072]
#BSUB -R select[mem>=3072]
#BSUB -R span[hosts=1]
#BSUB -q premium             # target queue for job execution
#BSUB -W 48:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o o.prsetexpr_nothreshold_withapoe
#BSUB -e e.prsetexpr_nothreshold_withapoe

# APOE-INCLUDED counterpart of scripts/run_prset.sh (full EUR sample).
# Differences vs run_prset.sh:
#   --msigdb      -> master.gmt                    (APOE gene retained)
#   --background  -> background_genes_withapoe.txt (APOE in the null background)
#   --pheno       -> ad_phenotype_residuals_withapoe.txt (APOE NOT regressed out)
#   --out         -> *_withapoe results prefix
#   (removed)     -> --x-range chr19:44409039-46412650  so APOE-region SNPs are
#                    clumped and scored instead of excluded.
# Run create_ad_phenotype_withapoe.R and create_background_withapoe.R first.

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
    --extract /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb18177-qc.snplist \
    --gtf /sc/arion/projects/paul_oreilly/lab/cotea02/project/data/reference/Homo_sapiens.GRCh37.75.gtf.gz \
    --keep /sc/arion/projects/paul_oreilly/data/ukb/genotyped/ukb18177-qc.fam \
    --msigdb /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/master.gmt \
    --num-auto 22 \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results/ad_case.control_prset_nothreshold_eur_withapoe \
    --pheno /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ad_phenotype_residuals_withapoe.txt \
    --pheno-col AD_resid \
    --set-perm 10000 \
    --target /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb18177_chr1.22 \
    --thread 48 \
    --ultra  \
    --background /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/background_genes_withapoe.txt:gene \
    --wind-3 35kb \
    --wind-5 35kb
