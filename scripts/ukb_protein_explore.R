# detection 
# can't assume everyone missing is actually zero 

library(tidyverse)
library(data.table)

dat <- fread('/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/residualized_NPX_withNA.txt')

table(is.na(dat$GFAP))
FALSE  TRUE
50227  2795

table(is.na(dat$NEFL))
FALSE  TRUE
50155  2867

# Histograms
pdf('hist_GFAP.pdf')
hist(dat$GFAP, main = "GFAP", xlab = "NPX")
dev.off()

pdf('hist_NEFL.pdf')
hist(dat$NEFL, main = "NEFL", xlab = "NPX")
dev.off()

# Q-Q plots
pdf('qq_GFAP.pdf')
qqnorm(dat$GFAP, main = "GFAP Q-Q"); qqline(dat$GFAP)
dev.off()

pdf('qq_NEFL.pdf')
qqnorm(dat$NEFL, main = "NEFL Q-Q"); qqline(dat$NEFL)
dev.off()

# Shapiro-Wilk (best for smaller samples; max n=5000 in R)
shapiro.test(sample(dat$GFAP[!is.na(dat$GFAP)], 5000))
shapiro.test(sample(dat$NEFL[!is.na(dat$NEFL)], 5000))

# Kolmogorov-Smirnov (works on full dataset)
ks.test(dat$GFAP[!is.na(dat$GFAP)], "pnorm",
        mean = mean(dat$GFAP, na.rm = TRUE),
        sd   = sd(dat$GFAP, na.rm = TRUE))




dat <- fread('/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/residualized_NPX_withNA_RINT.txt')

forprset <- dat[, .(sample_id, GFAP, NEFL)]

write.table(forprset, file = '/sc/arion/projects/paul_oreilly/data/ukb/protein/qc/data/residualized_NPX_withNA_RINT_GFAP.NEFL.txt',
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
