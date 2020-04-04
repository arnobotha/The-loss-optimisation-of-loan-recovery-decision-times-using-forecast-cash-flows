# =================== Sampling settings
# This script assumes an adequately prepared longitudinal dataset of loan performance over time, called [dat.perf]

Cut.Off.Cohort <- 201512 # for some samples
Cohort.start <- 200404

# --- Sampling

# Sample 1. Exclude readvances and earlier cohorts
dat.sample <- subset(dat.perf, Advance_Count == 1 & Cohort>=Cohort.start)
(1 - NROW(unique(dat.sample$Account)) / NROW(unique(dat.perf$Account))) * 100 #Lost 8% of base portfolio sample


# --- Subsampling
# - Sample 6. Exclude all accounts with a term significantly different to 240 months
dat.sample6 <- subset(dat.sample, Term_Mode >=180 & Term_Mode <= 240 )
(1 - NROW(unique(dat.sample6$Account)) / NROW(unique(dat.perf$Account))) * 100
# Lost 16.15% of base portfolio sample
(1 - NROW(unique(dat.sample6$Account)) / NROW(unique(dat.sample$Account))) * 100
# Lost 8.76% of sample1


# --- Data prep for samples: creating unique indexes
dat.temp <- unique(dat.sample6[,list(Account)])
dat.temp[, LoanID := 1:.N]; dat.sample6$LoanID <- NULL
dat.sample6 <- merge(dat.sample6,dat.temp,by="Account"); rm(dat.temp)
