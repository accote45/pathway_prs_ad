.header on
.mode csv
.output UKB_AD_MRI.csv

-- ============================================================================
-- UKB extraction: AD-relevant MRI endophenotypes for pathway PRS scan
-- ----------------------------------------------------------------------------
-- Schema: Choi Shing Wan SQLite layout
--   - One table per phenotype field, named fXXXX, with columns
--     (sample_id, pheno, instance, array)
--   - Participant table holds sample_id + withdrawn flag
--   - PCs stored as array in f22009 (1-indexed: array 1 = PC1, ... array 40 = PC40)
--
-- Run with:
--   sqlite3 ukb<ID>.db < extract_ukb_ad_mri.sql
--
-- Imaging fields are all instance = 2 (first imaging visit).
-- Baseline fields (sex, ethnicity, genetic QC) are instance = 0.
--
-- Phenotypes:
--   (1) Bilateral hippocampal volume        — FreeSurfer aseg
--   (2) Bilateral entorhinal cortex volume  — Desikan
--   (3) WMH volume                           — BIANCA
--   (4) AD-signature cortical thickness      — Jack 4-region bilateral mean
--                                              (entorhinal + inferior temporal
--                                              + middle temporal + fusiform)
-- ============================================================================


-- ─── Pivot the 20 genetic PCs from long → wide ────────────────────────────
CREATE TEMP TABLE PCs AS
SELECT sample_id,
    CAST(MAX(CASE WHEN array = 1  THEN pheno END) AS REAL) AS PC1,
    CAST(MAX(CASE WHEN array = 2  THEN pheno END) AS REAL) AS PC2,
    CAST(MAX(CASE WHEN array = 3  THEN pheno END) AS REAL) AS PC3,
    CAST(MAX(CASE WHEN array = 4  THEN pheno END) AS REAL) AS PC4,
    CAST(MAX(CASE WHEN array = 5  THEN pheno END) AS REAL) AS PC5,
    CAST(MAX(CASE WHEN array = 6  THEN pheno END) AS REAL) AS PC6,
    CAST(MAX(CASE WHEN array = 7  THEN pheno END) AS REAL) AS PC7,
    CAST(MAX(CASE WHEN array = 8  THEN pheno END) AS REAL) AS PC8,
    CAST(MAX(CASE WHEN array = 9  THEN pheno END) AS REAL) AS PC9,
    CAST(MAX(CASE WHEN array = 10 THEN pheno END) AS REAL) AS PC10,
    CAST(MAX(CASE WHEN array = 11 THEN pheno END) AS REAL) AS PC11,
    CAST(MAX(CASE WHEN array = 12 THEN pheno END) AS REAL) AS PC12,
    CAST(MAX(CASE WHEN array = 13 THEN pheno END) AS REAL) AS PC13,
    CAST(MAX(CASE WHEN array = 14 THEN pheno END) AS REAL) AS PC14,
    CAST(MAX(CASE WHEN array = 15 THEN pheno END) AS REAL) AS PC15,
    CAST(MAX(CASE WHEN array = 16 THEN pheno END) AS REAL) AS PC16,
    CAST(MAX(CASE WHEN array = 17 THEN pheno END) AS REAL) AS PC17,
    CAST(MAX(CASE WHEN array = 18 THEN pheno END) AS REAL) AS PC18,
    CAST(MAX(CASE WHEN array = 19 THEN pheno END) AS REAL) AS PC19,
    CAST(MAX(CASE WHEN array = 20 THEN pheno END) AS REAL) AS PC20
FROM f22009
GROUP BY sample_id;


-- ─── Final extract ────────────────────────────────────────────────────────
SELECT
    s.sample_id AS FID,
    s.sample_id AS IID,

    -- Demographics / visit info ──────────────────────────────────────────
    sex.pheno                       AS Sex,
    CAST(age.pheno     AS INT)      AS Age_Imaging,
    centre.pheno                    AS Centre_Imaging,
    imgdate.pheno                   AS Imaging_Date,
    eth.pheno                       AS Ethnic_Background,

    -- Genetic QC ─────────────────────────────────────────────────────────
    gsex.pheno                      AS Genetic_Sex,
    batch.pheno                     AS Genotyping_Batch,
    wbritish.pheno                  AS White_British,
    aneuploidy.pheno                AS Sex_Chr_Aneuploidy,
    kinship.pheno                   AS Genetic_Kinship,
    hetmiss.pheno                   AS Het_Missing_Outlier,

    -- Genetic PCs 1–20 ───────────────────────────────────────────────────
    pcs.PC1, pcs.PC2, pcs.PC3, pcs.PC4, pcs.PC5,
    pcs.PC6, pcs.PC7, pcs.PC8, pcs.PC9, pcs.PC10,
    pcs.PC11, pcs.PC12, pcs.PC13, pcs.PC14, pcs.PC15,
    pcs.PC16, pcs.PC17, pcs.PC18, pcs.PC19, pcs.PC20,

    -- Head size / ICV ────────────────────────────────────────────────────
    CAST(scaling.pheno AS REAL)     AS Head_Size_Scaling,
    CAST(etiv.pheno    AS REAL)     AS eTIV,

    -- ── PHENOTYPE 1: Hippocampal volume (aseg) ──────────────────────────
    CAST(hipp_l.pheno  AS REAL)     AS Hipp_L_Vol,
    CAST(hipp_r.pheno  AS REAL)     AS Hipp_R_Vol,
    (CAST(hipp_l.pheno AS REAL) + CAST(hipp_r.pheno AS REAL)) / 2.0
                                    AS Hipp_Mean_Vol,

    -- ── PHENOTYPE 2: Entorhinal volume (Desikan) ────────────────────────
    CAST(ento_l.pheno  AS REAL)     AS Ento_L_Vol,
    CAST(ento_r.pheno  AS REAL)     AS Ento_R_Vol,
    (CAST(ento_l.pheno AS REAL) + CAST(ento_r.pheno AS REAL)) / 2.0
                                    AS Ento_Mean_Vol,

    -- ── PHENOTYPE 3: WMH (BIANCA) ───────────────────────────────────────
    CAST(wmh.pheno     AS REAL)     AS WMH_Vol,

    -- ── PHENOTYPE 4: AD-signature cortical thickness ────────────────────
    -- Jack 4-region bilateral mean: entorhinal + inferior temporal +
    -- middle temporal + fusiform, averaged across L/R (8 values).
    -- Returns NULL if any constituent is missing (intended behavior).
    CAST(thk_ento_l.pheno     AS REAL) AS Thk_Ento_L,
    CAST(thk_ento_r.pheno     AS REAL) AS Thk_Ento_R,
    CAST(thk_inftemp_l.pheno  AS REAL) AS Thk_InfTemp_L,
    CAST(thk_inftemp_r.pheno  AS REAL) AS Thk_InfTemp_R,
    CAST(thk_midtemp_l.pheno  AS REAL) AS Thk_MidTemp_L,
    CAST(thk_midtemp_r.pheno  AS REAL) AS Thk_MidTemp_R,
    CAST(thk_fusiform_l.pheno AS REAL) AS Thk_Fusiform_L,
    CAST(thk_fusiform_r.pheno AS REAL) AS Thk_Fusiform_R,
    (
        CAST(thk_ento_l.pheno     AS REAL) + CAST(thk_ento_r.pheno     AS REAL) +
        CAST(thk_inftemp_l.pheno  AS REAL) + CAST(thk_inftemp_r.pheno  AS REAL) +
        CAST(thk_midtemp_l.pheno  AS REAL) + CAST(thk_midtemp_r.pheno  AS REAL) +
        CAST(thk_fusiform_l.pheno AS REAL) + CAST(thk_fusiform_r.pheno AS REAL)
    ) / 8.0                                AS AD_Sig_Thickness

FROM Participant s

    -- Baseline demographics ──────────────────────────────────────────────
    JOIN f31    sex     ON s.sample_id = sex.sample_id     AND sex.instance     = 0
    JOIN f21000 eth     ON s.sample_id = eth.sample_id     AND eth.instance     = 0

    -- Imaging visit anchor (instance 2). Restricts to imaged participants. ──
    JOIN f54    centre  ON s.sample_id = centre.sample_id  AND centre.instance  = 2
    JOIN f21003 age     ON s.sample_id = age.sample_id     AND age.instance     = 2
    JOIN f53    imgdate ON s.sample_id = imgdate.sample_id AND imgdate.instance = 2

    -- Hippocampus — LEFT JOIN so missing hippocampal data does not drop the
    -- subject from the other phenotypes (per-phenotype completeness handled in R)
    LEFT JOIN f26562 hipp_l  ON s.sample_id = hipp_l.sample_id  AND hipp_l.instance  = 2
    LEFT JOIN f26593 hipp_r  ON s.sample_id = hipp_r.sample_id  AND hipp_r.instance  = 2

    -- Other FreeSurfer measures — LEFT JOIN to tolerate per-region QC drops
    LEFT JOIN f26793 ento_l  ON s.sample_id = ento_l.sample_id  AND ento_l.instance  = 2
    LEFT JOIN f26894 ento_r  ON s.sample_id = ento_r.sample_id  AND ento_r.instance  = 2
    LEFT JOIN f26521 etiv    ON s.sample_id = etiv.sample_id    AND etiv.instance    = 2
    LEFT JOIN f25000 scaling ON s.sample_id = scaling.sample_id AND scaling.instance = 2

    -- WMH (different pipeline) — LEFT JOIN
    LEFT JOIN f25781 wmh ON s.sample_id = wmh.sample_id AND wmh.instance = 2

    -- AD-signature thickness components — LEFT JOIN
    LEFT JOIN f26760 thk_ento_l     ON s.sample_id = thk_ento_l.sample_id     AND thk_ento_l.instance     = 2
    LEFT JOIN f26861 thk_ento_r     ON s.sample_id = thk_ento_r.sample_id     AND thk_ento_r.instance     = 2
    LEFT JOIN f26763 thk_inftemp_l  ON s.sample_id = thk_inftemp_l.sample_id  AND thk_inftemp_l.instance  = 2
    LEFT JOIN f26864 thk_inftemp_r  ON s.sample_id = thk_inftemp_r.sample_id  AND thk_inftemp_r.instance  = 2
    LEFT JOIN f26769 thk_midtemp_l  ON s.sample_id = thk_midtemp_l.sample_id  AND thk_midtemp_l.instance  = 2
    LEFT JOIN f26870 thk_midtemp_r  ON s.sample_id = thk_midtemp_r.sample_id  AND thk_midtemp_r.instance  = 2
    LEFT JOIN f26761 thk_fusiform_l ON s.sample_id = thk_fusiform_l.sample_id AND thk_fusiform_l.instance = 2
    LEFT JOIN f26862 thk_fusiform_r ON s.sample_id = thk_fusiform_r.sample_id AND thk_fusiform_r.instance = 2

    -- Genetic QC (all baseline, instance 0) ──────────────────────────────
    LEFT JOIN f22001 gsex       ON s.sample_id = gsex.sample_id       AND gsex.instance       = 0
    LEFT JOIN f22000 batch      ON s.sample_id = batch.sample_id      AND batch.instance      = 0
    LEFT JOIN f22006 wbritish   ON s.sample_id = wbritish.sample_id   AND wbritish.instance   = 0
    LEFT JOIN f22019 aneuploidy ON s.sample_id = aneuploidy.sample_id AND aneuploidy.instance = 0
    LEFT JOIN f22021 kinship    ON s.sample_id = kinship.sample_id    AND kinship.instance    = 0
    LEFT JOIN f22027 hetmiss    ON s.sample_id = hetmiss.sample_id    AND hetmiss.instance    = 0

    -- Genetic PCs (pivoted in temp table above) ──────────────────────────
    LEFT JOIN PCs pcs ON s.sample_id = pcs.sample_id

WHERE s.withdrawn = 0
GROUP BY s.sample_id;
.quit

-- ============================================================================
-- DOWNSTREAM NOTES (R-side processing)
-- ----------------------------------------------------------------------------
-- Standard PRS-association covariate set (your usual UKB convention):
--   Age_Imaging + I(Age_Imaging^2) + Sex + Age_Imaging:Sex
--   + PC1..PC20 + factor(Genotyping_Batch) + factor(Centre_Imaging)
--   + eTIV    (volume phenotypes only — drop for AD_Sig_Thickness)
--
-- Phenotype transforms:
--   Hipp_Mean_Vol, Ento_Mean_Vol : residualize on eTIV, then rank-INT
--   WMH_Vol                       : log1p(WMH_Vol), then rank-INT
--   AD_Sig_Thickness              : no eTIV adjustment; rank-INT
--
-- Standard QC exclusions (apply in R for transparent N-flow logging):
--   Sex == Genetic_Sex
--   is.na(Sex_Chr_Aneuploidy)
--   is.na(Het_Missing_Outlier)
--   White_British == 1   (if EUR-only analysis)
--
-- AD-signature variant: this uses an unweighted mean (Jack et al. style).
-- For the surface-area-weighted Dickerson 2009 original, pull the matching
-- Desikan area fields (f26726/f26827 entorhinal, f26729/f26830 inferior
-- temporal, f26735/f26836 middle temporal, f26727/f26828 fusiform — confirm
-- against your BIG40 sheet) and compute
--   SUM(thickness_i * area_i) / SUM(area_i)
--
-- Instance-3 pooling: ~5k repeat-imaging subset exists at instance = 3.
-- If you want to gain those, add a second INNER JOIN on f26562 instance = 3
-- and COALESCE in the SELECT — but BIG40 GWAS used instance 2 only, so for
-- consistency with summary stats / pathway PRS scoring, leave as-is.
-- ============================================================================
