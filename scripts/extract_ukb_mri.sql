

.header on
.mode csv
.output hippo.csv              -- Output to file named Height.csv

SELECT  s.sample_id AS FID,     -- select sample_id from table 's' and call it FID 
s.sample_id AS IID,             -- select sample_id from table 's' and call it IID
age.pheno AS Age,               -- select pheno from table 'age' and call it Age
sex.pheno AS Sex,               -- select pheno from table 'sex' and call it Sex
bmi.pheno AS BMI,               -- select pheno from table 'bmi' and call it BMI
centre.pheno AS Centre          -- select pheno from table 'centre' and call it Centre

FROM    Participant s                            -- using table 'participant', now named as 's'
            JOIN    f21001 bmi ON 
                    s.sample_id=bmi.sample_id    -- join the BMI table by sample ID
                    AND bmi.instance = 0         -- only getting the baseline phenotype
            JOIN    f31 sex ON
                    s.sample_id=sex.sample_id    -- join the Sex table by sample ID
                    AND sex.instance = 0         -- only getting the baseline phenotype
            JOIN    f21003 age ON 
                    s.sample_id=age.sample_id    -- join the Age table by sample ID
                    AND age.instance = 0         -- only getting the baseline phenotype
            JOIN    f54 centre ON 
                    s.sample_id=centre.sample_id -- join the UKB assessment centre table by sample ID
                    AND centre.instance = 0      -- only getting the baseline phenotype
            WHERE   s.withdrawn = 0;             -- Exclude any samples who withdrawn their consent
.qui