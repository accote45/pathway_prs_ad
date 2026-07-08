# Protein phenotype exploration: missingness + normality of GFAP/NEFL NPX.
# Detection note: cannot assume everyone missing is actually zero.

library(tidyverse)
library(data.table)

dat <- fread('/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/residualized_NPX_withNA_RINT.txt')

# Missingness
table(is.na(dat$GFAP))
table(is.na(dat$NEFL))

# Histograms
pdf('hist_GFAP.pdf'); hist(dat$GFAP, main = "GFAP", xlab = "NPX"); dev.off()
pdf('hist_NEFL.pdf'); hist(dat$NEFL, main = "NEFL", xlab = "NPX"); dev.off()

# Q-Q plots
pdf('qq_GFAP.pdf'); qqnorm(dat$GFAP, main = "GFAP Q-Q"); qqline(dat$GFAP); dev.off()
pdf('qq_NEFL.pdf'); qqnorm(dat$NEFL, main = "NEFL Q-Q"); qqline(dat$NEFL); dev.off()

# Normality tests (Shapiro-Wilk max n = 5000 in R; KS on the full data)
shapiro.test(sample(dat$GFAP[!is.na(dat$GFAP)], 5000))
shapiro.test(sample(dat$NEFL[!is.na(dat$NEFL)], 5000))
ks.test(dat$GFAP[!is.na(dat$GFAP)], "pnorm",
        mean = mean(dat$GFAP, na.rm = TRUE),
        sd   = sd(dat$GFAP, na.rm = TRUE))

# GFAP/NEFL-only phenotype table (for a PRSet sensitivity check)
forprset <- dat[, .(sample_id, GFAP, NEFL)]
write.table(forprset,
            file = '/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/residualized_NPX_withNA_RINT_GFAP.NEFL.txt',
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
