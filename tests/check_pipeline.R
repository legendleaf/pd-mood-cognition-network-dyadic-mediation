# Synthetic-data checks; no study data are read.
repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
for (file in list.files("script", pattern = "[.]R$", full.names = TRUE)) parse(file)
root <- tempfile("pd_test_")
dir.create(root)
input <- file.path(root, "inputs")
dir.create(input)
set.seed(123)
n <- 40L
domains <- paste0("PDQ39_", c("mobility", "adl", "emotional", "stigma", "social",
                             "cognition", "communication", "bodily"))
p <- data.frame(ID = seq_len(n), NMS_mood = sample(0:5, n, TRUE), WOQ9_total = sample(0:9, n, TRUE),
  age = sample(35:80, n, TRUE), sex_male = sample(0:1, n, TRUE), disease_duration = runif(n, 1, 15))
p[domains] <- replicate(8, runif(n, 0, 75), simplify = FALSE)
p$PDQ39_SI <- rowMeans(p[domains])
path <- file.path(input, "analysis_patients.csv")
write.csv(p, path, row.names = FALSE)
before <- tools::md5sum(path)
Sys.setenv(PD_DATA_DIR = input, PD_OUTPUT_DIR = file.path(root, "results"))
Sys.unsetenv(c("PD_ITEM_WORKBOOK", "PD_CAREGIVER_WORKBOOK", "PD_FUNCTION_CSV", "PD_MODEL_ARCHIVE_DIR"))
source("script/Step0_DataPreparation.R")
stopifnot(max(abs(pd_patients$PDQ39_SI - rowMeans(p[domains]))) < 1e-10,
  max(abs(pd_patients$PDQ39_six_domain_mean - rowMeans(p[setdiff(domains,
    c("PDQ39_emotional", "PDQ39_cognition"))]))) < 1e-10)
source("script/Step4_PatientSensitivity.R")
results <- read.csv(file.path(PD_OUTPUT_DIR, "tables", "patient_regressions_HC3.csv"))
stopifnot(nrow(results) == 10L, all(results$n == n), all(results$df == n - 5L))
stopifnot(identical(before, tools::md5sum(path)))
legacy <- p
legacy$PDQ39_SI <- legacy$PDQ39_SI + 1
write.csv(legacy, path, row.names = FALSE)
failure <- try(source("script/Step0_DataPreparation.R"), silent = TRUE)
stopifnot(inherits(failure, "try-error"))
write.csv(p, path, row.names = FALSE)
unsafe <- file.path(input, "must_not_be_created")
Sys.setenv(PD_OUTPUT_DIR = unsafe)
failure <- try(source("script/Step0_DataPreparation.R"), silent = TRUE)
stopifnot(inherits(failure, "try-error"), !dir.exists(unsafe))
cat("Synthetic scoring, HC3 regressions, legacy rejection and output guard passed.\n")
