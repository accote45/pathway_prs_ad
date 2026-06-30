.header on
.mode csv
.output UKB_AD_MRI_FSL.csv

-- ============================================================================
-- UKB extraction: AD-relevant MRI endophenotypes — FSL-pipeline subset
-- ----------------------------------------------------------------------------
-- Target DB: ukb18177.db (FSL FAST/FIRST IDPs loaded; FreeSurfer 26500+ NOT).
-- Schema: Choi Shing Wan SQLite layout
--   - One table per phenotype field, named fXXXX, columns
--     (sample_id, pheno, instance, array)
--   - participant table holds sample_id + withdrawn flag
--   - PCs stored as array in f22009 (array 0 = PC1, ... , array 39 = PC40)
--
-- Run with:
--   sqlite3 ukb18177.db < extract_ukb_hipp_wmh.sql
--
-- Imaging fields are all instance = 2 (first imaging visit).
-- Baseline fields (sex, ethnicity, genetic QC) are instance = 0.
--
-- Phenotypes available in this DB:
--   (1) Bilateral hippocampal volume — FSL FIRST  (f25019 / f25020)
--   (3) WMH volume                    — T1+T2_FLAIR (f25781)
-- NOT available here (need a FreeSurfer-loaded DB): entorhinal volume,
-- AD-signature cortical thickness. eTIV (f26521) absent too; head-size
-- correction uses the FSL head-size scaling factor f25000 instead.
-- ============================================================================


-- ─── Pivot the 20 genetic PCs from long → wide ────────────────────────────
CREATE TEMP TABLE PCs AS
SELECT sample_id,
    CAST(MAX(CASE WHEN array = 0  THEN pheno END) AS REAL) AS PC1,
    CAST(MAX(CASE WHEN array = 1  THEN pheno END) AS REAL) AS PC2,
    CAST(MAX(CASE WHEN array = 2  THEN pheno END) AS REAL) AS PC3,
    CAST(MAX(CASE WHEN array = 3  THEN pheno END) AS REAL) AS PC4,
    CAST(MAX(CASE WHEN array = 4  THEN pheno END) AS REAL) AS PC5,
    CAST(MAX(CASE WHEN array = 5  THEN pheno END) AS REAL) AS PC6,
    CAST(MAX(CASE WHEN array = 6  THEN pheno END) AS REAL) AS PC7,
    CAST(MAX(CASE WHEN array = 7  THEN pheno END) AS REAL) AS PC8,
    CAST(MAX(CASE WHEN array = 8  THEN pheno END) AS REAL) AS PC9,
    CAST(MAX(CASE WHEN array = 9  THEN pheno END) AS REAL) AS PC10,
    CAST(MAX(CASE WHEN array = 10 THEN pheno END) AS REAL) AS PC11,
    CAST(MAX(CASE WHEN array = 11 THEN pheno END) AS REAL) AS PC12,
    CAST(MAX(CASE WHEN array = 12 THEN pheno END) AS REAL) AS PC13,
    CAST(MAX(CASE WHEN array = 13 THEN pheno END) AS REAL) AS PC14,
    CAST(MAX(CASE WHEN array = 14 THEN pheno END) AS REAL) AS PC15,
    CAST(MAX(CASE WHEN array = 15 THEN pheno END) AS REAL) AS PC16,
    CAST(MAX(CASE WHEN array = 16 THEN pheno END) AS REAL) AS PC17,
    CAST(MAX(CASE WHEN array = 17 THEN pheno END) AS REAL) AS PC18,
    CAST(MAX(CASE WHEN array = 18 THEN pheno END) AS REAL) AS PC19,
    CAST(MAX(CASE WHEN array = 19 THEN pheno END) AS REAL) AS PC20
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

    -- Head size (FSL T1 → standard-space scaling factor; eTIV unavailable) ─
    CAST(scaling.pheno AS REAL)     AS Head_Size_Scaling,

    -- ── PHENOTYPE 1: Hippocampal volume (FSL FIRST) ─────────────────────
    CAST(hipp_l.pheno  AS REAL)     AS Hipp_L_Vol,
    CAST(hipp_r.pheno  AS REAL)     AS Hipp_R_Vol,
    (CAST(hipp_l.pheno AS REAL) + CAST(hipp_r.pheno AS REAL)) / 2.0
                                    AS Hipp_Mean_Vol,

    -- ── PHENOTYPE 3: WMH (T1 + T2_FLAIR) ────────────────────────────────
    CAST(wmh.pheno     AS REAL)     AS WMH_Vol

FROM participant s

    -- Baseline demographics ──────────────────────────────────────────────
    JOIN f31    sex     ON s.sample_id = sex.sample_id     AND sex.instance     = 0
    JOIN f21000 eth     ON s.sample_id = eth.sample_id     AND eth.instance     = 0

    -- Imaging visit anchor (instance 2). Restricts to imaged participants. ──
    JOIN f54    centre  ON s.sample_id = centre.sample_id  AND centre.instance  = 2
    JOIN f21003 age     ON s.sample_id = age.sample_id     AND age.instance     = 2
    JOIN f53    imgdate ON s.sample_id = imgdate.sample_id AND imgdate.instance = 2

    -- Hippocampus (FIRST) — LEFT JOIN so missing hippocampal data does not
    -- drop the subject from WMH (per-phenotype completeness handled in R)
    LEFT JOIN f25019 hipp_l  ON s.sample_id = hipp_l.sample_id  AND hipp_l.instance  = 2
    LEFT JOIN f25020 hipp_r  ON s.sample_id = hipp_r.sample_id  AND hipp_r.instance  = 2

    -- Head-size scaling factor — LEFT JOIN
    LEFT JOIN f25000 scaling ON s.sample_id = scaling.sample_id AND scaling.instance = 2

    -- WMH — LEFT JOIN
    LEFT JOIN f25781 wmh ON s.sample_id = wmh.sample_id AND wmh.instance = 2

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
--   + Head_Size_Scaling   (volume phenotype only — drop for WMH if log-ratio)
--
-- Phenotype transforms:
--   Hipp_Mean_Vol : residualize on Head_Size_Scaling (FSL norm factor), rank-INT
--                   -- alternatively multiply raw volume by Head_Size_Scaling
--                      to put it in standard space before residualizing.
--   WMH_Vol       : log1p(WMH_Vol), then rank-INT
--
-- Standard QC exclusions (apply in R for transparent N-flow logging):
--   Sex == Genetic_Sex
--   is.na(Sex_Chr_Aneuploidy)
--   is.na(Het_Missing_Outlier)
--   White_British == 1   (if EUR-only analysis)
--
-- NOTE: hippocampus here is FSL FIRST (f25019/f25020), not FreeSurfer aseg.
-- Entorhinal volume and AD-signature thickness are NOT in ukb18177.db; pull
-- those from a DB with the FreeSurfer fields (f265xx–f273xx, categories 190-196)
-- and merge on IID. BIG40 GWAS for FIRST hippocampus used instance 2.
-- ============================================================================
