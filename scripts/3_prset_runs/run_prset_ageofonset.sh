


#BSUB -L /bin/sh
#BSUB -n 48
#BSUB -J prsetexpr_ageofonset
#BSUB -R rusage[mem=3072]
#BSUB -R select[mem>=3072]
#BSUB -R span[hosts=1]
#BSUB -q premium             # target queue for job execution
#BSUB -W 10:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o o.prsetexpr_ageofonset
#BSUB -e e.prsetexpr_ageofonset

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
    --keep /sc/arion/projects/paul_oreilly/data/ukb/genotyped/ukb18177-qc.fam \
    --msigdb /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/master_no_apoe.gmt \
    --num-auto 22 \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/results/ad_ageofonset_prset_nothreshold_eur \
    --pheno /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/phenotypes/ad_ageofonset_phenotype_residuals.txt \
    --pheno-col age_of_onset_resid \
    --set-perm 10000 \
    --target /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb18177_chr1.22 \
    --thread 48 \
    --ultra  \
    --background /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/pathways/background_genes.txt:gene \
    --wind-3 35kb \
    --wind-5 35kb
