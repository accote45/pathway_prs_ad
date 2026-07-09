library(DBI)
library(dplyr)
library(data.table)
con <-
    DBI::dbConnect(RSQLite::SQLite(), dbname = "/sc/arion/projects/data-ark/ukb/application/ukb18177/phenotype/ukb18177-v2.db")

# Get maternal age --------------------------------------------------------
maternal_age <- as.data.table(tbl(con, "f1845"))
setnames(maternal_age, "pheno", "age")
maternal_age_death <- as.data.table(tbl(con, "f3526"))
setnames(maternal_age_death, "pheno", "death_age")

# Get maternal phenotype --------------------------------------------------
maternal_pheno <-
    as.data.table(tbl(con, "f20110"))[pheno %in% c(-17, -11, -13, 1, 2, 6, 8, 9, 10)]
# Basically, any uncertain answer (-11, -13) is sest to -1
maternal_pheno[, status := -1]
# -17 is none of the above (group 1)
maternal_pheno[pheno == -17, status := -2]
maternal_pheno[pheno > 0, status := 0]
maternal_pheno[pheno == 10, status := 1]
# Get the first instance where each maternal status is reported
maternal_pheno <-
    maternal_pheno[, .(instance = min(instance)), by = c("sample_id", "status")]
# We also want the largest phenotype from the data
maternal_pheno <-
    maternal_pheno[, .SD[status == max(status)], by = sample_id]
maternal_pheno <-
    merge(
        maternal_pheno,
        maternal_age,
        by = c("sample_id", "instance"),
        all.x = TRUE
    )
maternal_pheno <-
    merge(
        maternal_pheno,
        maternal_age_death,
        by = c("sample_id", "instance"),
        all.x = TRUE
    )

maternal_pheno[, maternal_age := pmin(age, death_age, na.rm = TRUE)]
setnames(maternal_pheno, "status", "maternal_status")
# Get paternal age --------------------------------------------------------
paternal_age <- as.data.table(tbl(con, "f2946"))
setnames(paternal_age, "pheno", "age")
paternal_age_death <- as.data.table(tbl(con, "f1807"))
setnames(paternal_age_death, "pheno", "death_age")

# Get paternal phenotype --------------------------------------------------
paternal_pheno <-
    as.data.table(tbl(con, "f20107"))[pheno %in% c(-17, -11, -13, 1, 2, 6, 8, 9, 10)]
# Basically, any uncertain answer (-11, -13) is sest to -1
paternal_pheno[, status := -1]
# -17 is none of the above (group 1)
paternal_pheno[pheno == -17, status := -2]
paternal_pheno[pheno > 0, status := 0]
paternal_pheno[pheno == 10, status := 1]
# Get the first instance where each maternal status is reported
paternal_pheno <-
    paternal_pheno[, .(instance = min(instance)), by = c("sample_id", "status")]
# We also want the largest phenotype from the data
paternal_pheno <-
    paternal_pheno[, .SD[status == max(status)], by = sample_id]
paternal_pheno <-
    merge(
        paternal_pheno,
        paternal_age,
        by = c("sample_id", "instance"),
        all.x = TRUE
    )
paternal_pheno <-
    merge(
        paternal_pheno,
        paternal_age_death,
        by = c("sample_id", "instance"),
        all.x = TRUE
    )
paternal_pheno[, paternal_age := pmin(age, death_age, na.rm = TRUE)]
setnames(paternal_pheno, "status", "paternal_status")

# Extract sample phenotype -----------------------------------------------

exclusion  <- c("\"G310\"", "\"F020\"")
required <- c("\"G30", "\"F00")
icd10 <-
    c(main = "f41202",
      secondary = "f41204",
      summary  = "f41270")
icd10_dat <- NULL
for (i in seq_along(icd10)) {
    cur_icd10 <-
        as.data.table(tbl(con, icd10[i]))
    cur_icd10 <- cur_icd10[pheno %in% exclusion |
                               pheno %like% paste(required, sep = "|", collapse = "|")]
    cur_icd10[, source := names(icd10)[i]]
    icd10_dat <- rbind(icd10_dat, cur_icd10)
}

icd10_death <- c(main = "f40001", secondary = "f40002")
icd10_death_dat <- NULL
for (i in seq_along(icd10_death)) {
    cur_icd10 <- as.data.table(tbl(con, icd10_death[i]))
    
    cur_icd10 <- cur_icd10[pheno %in% exclusion |
                               pheno %like% paste(required, sep = "|", collapse = "|")]
    cur_icd10[, source := names(icd10_death)[i]]
    cur_icd10[, array := NULL]
    icd10_death_dat <- rbind(icd10_death_dat, cur_icd10)
}

# Now do filtering --------------------------------------------------------
# This ignore the source of data, sometimes we might only want primary.
# If that is the case, filter here
icd10_dat <- unique(icd10_dat[, source := NULL])
icd10_death_dat <- unique(icd10_death_dat[, source := NULL])
# Get age at diagnoses ----------------------------------------------------
birth_year <- as.data.table(tbl(con, "f34"))
birth_year[, pheno := as.numeric(pheno)]
setnames(birth_year, "pheno", "birth_year")
icd10_date_diagnosed <- as.data.table(tbl(con, "f41280"))
setnames(icd10_date_diagnosed, "pheno", "date_diagnosed")
icd10_pheno <-
    merge(icd10_dat,
          icd10_date_diagnosed,
          by  = c("sample_id", "instance", "array"))
icd10_pheno[, date_diagnosed := as.Date(date_diagnosed)]
icd10_pheno[, c("array", "instance") := NULL]
icd10_pheno <- merge(icd10_pheno, birth_year, all.x = TRUE)
icd10_pheno[, age_of_diagnosis := year(date_diagnosed) - birth_year]
# Get age at death --------------------------------------------------------
age_of_death <- as.data.table(tbl(con, "f40007"))
setnames(age_of_death, "pheno", "age_of_death")
icd10_death_info <-
    merge(icd10_death_dat,
          age_of_death,
          by  = c("sample_id", "instance"))
icd10_death_info[, instance := NULL]
icd10_death_info[, age_of_death := as.numeric(age_of_death)]


# GP codes ----------------------------------------------------------------


read3_exclude <- c(
    '.E115',
    '.E116',
    '.G78.',
    'E004.',
    'E0040',
    'E0041',
    'E0042',
    'E0043',
    'E004z',
    'Eu01.',
    'Eu010',
    'Eu011',
    'Eu012',
    'Eu013',
    'Eu01y',
    'Eu01z',
    'Eu020',
    'Eu025',
    'F111.',
    'F116.',
    'F118.',
    'F11x2',
    'F11y2',
    'F21y2',
    'X0034',
    'X0035',
    'X0036',
    'X0037',
    'X0039',
    'X003A',
    'X003B',
    'X003C',
    'X003D',
    'X003m',
    'X003R',
    'X003T',
    'X003V',
    'X003W',
    'Xa0lH',
    'Xa0sC',
    'Xa0sE',
    'XaE74',
    'XaKyY',
    'XE1Xs'
)
read2_exclude <-
    c (
        'E004.',
        'E0040',
        'E0041',
        'E0042',
        'E0043',
        'E004z',
        'E012.',
        'Eu01.',
        'Eu010',
        'Eu011',
        'Eu012',
        'Eu013',
        'Eu01y',
        'Eu01z',
        'Eu020',
        'Eu025',
        'F111.',
        'F116.',
        'F118.',
        'F11x2',
        'F11y2',
        'F21y2'
    )
read3_include <- c (
    '.F21Z',
    'Eu00.',
    'Eu000',
    'Eu001',
    'Eu002',
    'Eu00z',
    'F110.',
    'F1100',
    'F1101',
    'Fyu30',
    'X002x',
    'X002y',
    'X002z',
    'X0030',
    'X0031',
    'X0032',
    'X0033',
    'X003G',
    'XaIKB',
    'XaIKC',
    'XE17j'
)
read2_include <- c('Eu00.',
                   'Eu000',
                   'Eu001',
                   'Eu002',
                   'Eu00z',
                   'F110.',
                   'F1100',
                   'F1101',
                   'Fyu30')

# Extract information from GP data ----------------------------------------

gp_exclude <-
    tbl(con, "gp_clinical") |>
    filter(read3 %in% read3_exclude |
               read2 %in% read2_exclude) |>
    select(sample_id) |>
    distinct() |>
    collect() |>
    as.data.table()
gp_include <-  tbl(con, "gp_clinical") |>
    filter(read3 %in% read3_include |
               read2 %in% read2_include) |>
    select(sample_id, date_event) |>
    distinct() |>
    collect() |>
    as.data.table()
gp_include[, date_event := as.Date(date_event)]
gp_include <- merge(gp_include, birth_year)
gp_include[, age_of_diagnosis := year(date_event) - birth_year]


# Get all excluded samples ------------------------------------------------
exclude_icd <-
    icd10_pheno[pheno %in% exclusion, sample_id] |> unique()
exclude_death <-
    icd10_death_info[pheno %in% exclusion, sample_id] |> unique()
exclude_samples <-
    unique(c(exclude_icd, exclude_death, gp_exclude[, sample_id]))
# Combine all data source -------------------------------------------------
ad_cases <-
    rbind(gp_include[!sample_id %in% exclude_samples, c("sample_id", "age_of_diagnosis")],
          icd10_pheno[!sample_id %in% exclude_samples, c("sample_id", "age_of_diagnosis")])
ad_cases <-
    ad_cases[, .(age_of_diagnosis = min(age_of_diagnosis)), by = sample_id]
ad_cases <-
    merge(ad_cases, icd10_death_info[, c("sample_id", "age_of_death")], all = TRUE)
ad_cases[, age_of_death := ceiling(age_of_death)]
ad_cases[, age_of_diagnosis := pmin(age_of_death, age_of_diagnosis, na.rm = TRUE)]
ad_cases <- ad_cases[, c("sample_id", "age_of_diagnosis")]

# Adoption information ----------------------------------------------------
adoption <- tbl(con, "f1767") |>
    filter(pheno == "1") |>
    select(sample_id) |>
    distinct() |>
    collect() |>
    as.data.table()
# Have checked, if someone state that they are adopted, they will always be adopted


# Get age at recruitment --------------------------------------------------

age_recruitment <- as.data.table(tbl(con, "f21022"))
age_recruitment <-
    age_recruitment[, .(last_recorded_age = max(as.numeric(pheno)), baseline_age = min(as.numeric(pheno))), by = sample_id]

# Get assessment centre --------------------------------------------------
centre <- tbl(con, "f54") |>
    filter(instance=="0") |>
    select("sample_id", "pheno") |>
    distinct() |> 
    collect() |>
    as.data.table()
setnames(centre, "pheno", "assessment_centre")
centre[,record_system := "England"]
centre[assessment_centre %in% c(11004, 11005), record_system := "Scotland"]
centre[assessment_centre %in% c(11022, 11003, 11023), record_system := "Walsh/Irish"]
last_diagnose_date <- merge(centre, icd10_date_diagnosed)
last_diagnose_date <- last_diagnose_date[, .(date = max(date_diagnosed)), by = record_system]

sample_censor_date <- merge(centre, last_diagnose_date)
sample_censor_date[,date := as.Date(date)]
sample_censor_date <- merge(sample_censor_date, birth_year, by = "sample_id")
sample_censor_date[,final_diagnose_age := year(date) - birth_year]
sample_censor_date <- sample_censor_date[,c("sample_id", "final_diagnose_age")]
# Get samples -------------------------------------------------------------

all_samples <- tbl(con, "Participant") |>
    filter(withdrawn == 0) |>
    select(sample_id) |>
    as.data.table()
# Remove samples
all_samples <- all_samples[!sample_id %in% exclude_samples]

# Now construct required data.table ---------------------------------------
all_samples <- merge(all_samples, age_recruitment, all.x = TRUE)
# Likely due to withdrawn samples
all_samples <- all_samples[!is.na(last_recorded_age) & !is.na(baseline_age)]
all_samples <- merge(all_samples, ad_cases, all.x = TRUE)
all_samples <- merge(all_samples, sample_censor_date, all.x = TRUE)
all_samples[, AD := 0]
all_samples[!is.na(age_of_diagnosis), AD := 1]
all_samples[, censor_age := ifelse(AD == 1, age_of_diagnosis, pmax(last_recorded_age, final_diagnose_age))]
all_samples[, adopted := sample_id %in% adoption[, sample_id]]
all_samples <-
    merge(all_samples, maternal_pheno[, c("sample_id", "maternal_age", "maternal_status")], all.x = TRUE)
all_samples <-
    merge(all_samples, paternal_pheno[, c("sample_id", "paternal_age", "paternal_status")], all.x = TRUE)
all_samples[,prevalenc := ]
fwrite(all_samples, "AD-R.csv")
