


#BSUB -L /bin/sh
#BSUB -n 48
#BSUB -J prset_afr_t2d_rppr
#BSUB -R rusage[mem=2048]
#BSUB -R select[mem>=2048]
#BSUB -R span[hosts=1]
#BSUB -q premium             # target queue for job execution
#BSUB -W 5:00                # wall clock limit for job
#BSUB -P acc_psychgen             # project to charge time
#BSUB -o o.prset_afr_t2d_rppr
#BSUB -e e.prset_afr_t2d_rppr

ml R
Rscript /sc/arion/projects/psychgen/cotea02_prset/PRSice.R --prsice /sc/arion/projects/psychgen/cotea02_prset/PRSice_linux \
    --snp SNP \
    --pvalue P \
    --beta \
    --stat BETA \
    --fastscore \
    --bar-levels 1 \
    --a1 A1 \
    --a2 A2 \
    #--base /sc/arion/projects/paul_oreilly/data/GWASs/NonBiobanks/qced_data/t2d/t2d_hg19.txt.gz \
    --binary-target F \
    --clump-kb 1000kb \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --extract /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_chr_maf01_qc.snps \
    --gtf /sc/arion/projects/paul_oreilly/lab/cotea02/project/data/reference/Homo_sapiens.GRCh37.75.gtf.gz \
    --keep /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/afr_sample_ids_80pc.txt \
    --msigdb /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/master_no_apoe.gmt \
    --num-auto 22 \
    --out /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/t2d_case.control_prset_nothreshold_afr \
    --pheno /sc/arion/projects/psychgen/cotea02_prset/geneoverlap_nf/data/ukb/ukb_phenofile_forprset.txt \
    --pheno-col t2d_resid \
    --set-perm 10000 \
    --target /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/data/ukb_afr_maf01 \
    --thread 48 \
    --ultra  \
    --background /sc/arion/projects/paul_oreilly/lab/cotea02/pathway_prs_ad/scripts/output/background_genes.txt:gene \
    --wind-3 35kb \
    --wind-5 35kb

#     --x-range chr19:45000000-46000000 \
