# Shared input checks and scoring. Source from the repository root.

extra_library <- Sys.getenv("PD_R_LIBRARY")
if (nzchar(extra_library)) .libPaths(c(strsplit(extra_library, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))

PD_DATA_DIR <- normalizePath(Sys.getenv("PD_DATA_DIR", "private"), winslash = "/", mustWork = TRUE)
PD_OUTPUT_DIR <- Sys.getenv("PD_OUTPUT_DIR", "results_network_analysis")
inside <- function(path, parent) identical(tolower(path), tolower(parent)) ||
  startsWith(tolower(path), paste0(tolower(parent), "/"))
if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", PD_OUTPUT_DIR)) PD_OUTPUT_DIR <- file.path(getwd(), PD_OUTPUT_DIR)
PD_OUTPUT_DIR <- normalizePath(PD_OUTPUT_DIR, winslash = "/", mustWork = FALSE)
if (inside(PD_OUTPUT_DIR, PD_DATA_DIR)) stop("Choose an output directory outside the input directory.")
for (variable in c("PD_ITEM_WORKBOOK", "PD_CAREGIVER_WORKBOOK", "PD_FUNCTION_CSV")) {
  private_input_path <- Sys.getenv(variable)
  if (nzchar(private_input_path) && inside(PD_OUTPUT_DIR, dirname(normalizePath(private_input_path, winslash = "/", mustWork = TRUE)))) {
    stop("Private input files and outputs must use separate directories.")
  }
}
dir.create(PD_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
PD_OUTPUT_DIR <- normalizePath(PD_OUTPUT_DIR, winslash = "/", mustWork = TRUE)
for (sub in c("figures/main", "figures/supplement", "tables", "data_processed", "logs")) {
  dir.create(file.path(PD_OUTPUT_DIR, sub), recursive = TRUE, showWarnings = FALSE)
}

pd_domains <- paste0("PDQ39_", c("mobility", "adl", "emotional", "stigma", "social",
                                "cognition", "communication", "bodily"))
pd_network_nodes <- c("PDQ39_mobility", "PDQ39_adl", "PDQ39_emotional", "PDQ39_cognition",
  "PDQ39_social", "PDQ39_communication", "PDQ39_stigma", "NMS_sleep", "NMS_mood",
  "NMS_gi", "NMS_urinary", "NMS_cardiovascular", "NMS_perception", "NMS_sexual",
  "NMS_misc", "WOQ9_total")
pd_input_hashes <- character()
pd_require <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}
pd_z <- function(x) {
  if (any(!is.finite(x)) || sd(x) == 0) stop("Standardization requires finite, nonconstant values.")
  as.numeric(scale(x))
}
pd_count <- function(variable, default) {
  value <- as.integer(Sys.getenv(variable, as.character(default)))
  if (is.na(value) || value < 2L) stop(variable, " must be an integer >= 2.")
  value
}
pd_read <- function(path) {
  if (!file.exists(path)) stop("An approved private input file is missing.")
  pd_input_hashes[path] <<- unname(tools::md5sum(path))
  data <- read.csv(path, check.names = FALSE)
  pd_require(data, "ID")
  if (anyNA(data$ID) || anyDuplicated(data$ID)) stop("Input IDs must be nonmissing and unique.")
  data
}
pd_workbook <- function(variable) {
  path <- Sys.getenv(variable)
  if (!nzchar(path)) return(NULL)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (inside(PD_OUTPUT_DIR, dirname(path))) stop("Workbook and output directories must be separate.")
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read item workbooks.")
  pd_input_hashes[path] <<- unname(tools::md5sum(path))
  data <- as.data.frame(readxl::read_excel(path, sheet = 1))
  pd_require(data, "ID")
  if (anyNA(data$ID) || anyDuplicated(data$ID)) stop("Workbook IDs must be nonmissing and unique.")
  data
}
pd_complete <- function(data, columns, analysis) {
  pd_require(data, columns)
  selected <- data[columns]
  keep <- complete.cases(selected)
  numeric <- vapply(selected, is.numeric, logical(1))
  if (any(numeric)) keep <- keep & apply(selected[numeric], 1, function(x) all(is.finite(x)))
  entry <- data.frame(analysis = analysis, available = nrow(data), analyzed = sum(keep),
                      excluded_missing_or_nonfinite = sum(!keep))
  path <- file.path(PD_OUTPUT_DIR, "tables", "analysis_samples.csv")
  if (file.exists(path)) {
    previous <- read.csv(path)
    entry <- rbind(previous[previous$analysis != analysis, ], entry)
  }
  write.csv(entry, path, row.names = FALSE)
  data[keep, , drop = FALSE]
}
pd_finish <- function(step) {
  after <- tools::md5sum(names(pd_input_hashes))
  if (!identical(unname(after), unname(pd_input_hashes))) stop("An input file changed during execution.")
  write.csv(data.frame(file = basename(names(after)), md5 = unname(after)),
    file.path(PD_OUTPUT_DIR, "logs", paste0(step, "_input_hashes.csv")), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(PD_OUTPUT_DIR, "logs", paste0(step, "_session.txt")))
}

pd_patients <- pd_read(file.path(PD_DATA_DIR, "analysis_patients.csv"))
pd_require(pd_patients, c(pd_domains, "NMS_mood", "WOQ9_total", "age", "sex_male", "disease_duration"))
pd_items <- pd_workbook("PD_ITEM_WORKBOOK")
if (!is.null(pd_items)) {
  index <- match(pd_patients$ID, pd_items$ID)
  if (anyNA(index)) stop("Some patient IDs are absent from the item workbook.")
  pd_items <- pd_items[index, , drop = FALSE]
  item_names <- paste0("PDQ39_", 1:39)
  pd_require(pd_items, c(item_names, "PDQ39"))
  item_matrix <- as.matrix(pd_items[item_names])
  if (any(!is.finite(item_matrix)) || any(!item_matrix %in% 1:5)) stop("PDQ item workbook must use codes 1-5.")
  scored <- item_matrix - 1
  if (max(abs(rowSums(scored) - pd_items$PDQ39)) > 1e-10) stop("PDQ items and the finalized workbook total disagree.")
  spans <- list(1:10, 11:16, 17:22, 23:26, 27:29, 30:33, 34:36, 37:39)
  domains <- as.data.frame(lapply(spans, function(j) rowMeans(scored[, j, drop = FALSE]) * 25))
  names(domains) <- pd_domains
  other <- setdiff(pd_domains, "PDQ39_social")
  if (max(abs(as.matrix(domains[other]) - as.matrix(pd_patients[other]))) > 1e-8) {
    stop("Finalized workbook and analysis CSV disagree outside the documented social-domain correction.")
  }
  social_delta <- pd_patients$PDQ39_social - domains$PDQ39_social
  if (any(abs(social_delta) > 1e-8 & abs(social_delta - 100/12) > 1e-8)) {
    stop("Unexpected social-domain discrepancy; check the approved inputs.")
  }
  pd_patients[pd_domains] <- domains
  pd_patients$PDQ39_total <- pd_items$PDQ39
  pd_patients$PDQ39_total_0_100 <- pd_items$PDQ39 / 156 * 100
} else {
  pd_require(pd_patients, "PDQ39_SI")
  expected <- rowMeans(pd_patients[pd_domains], na.rm = FALSE)
  if (any(abs(pd_patients$PDQ39_SI - expected) > 1e-8, na.rm = TRUE)) {
    stop("Legacy SI detected. Supply PD_ITEM_WORKBOOK or approved R1 domain-score CSVs.")
  }
}
if (any(as.matrix(pd_patients[pd_domains]) < 0 | as.matrix(pd_patients[pd_domains]) > 100, na.rm = TRUE)) {
  stop("PDQ domain scores must lie between 0 and 100.")
}
pd_patients$PDQ39_SI <- rowMeans(pd_patients[pd_domains], na.rm = FALSE)
pd_patients$PDQ39_six_domain_mean <- rowMeans(pd_patients[setdiff(pd_domains,
  c("PDQ39_emotional", "PDQ39_cognition"))], na.rm = FALSE)
pd_patients$FMCS <- 100 - (10 * pd_patients$PDQ39_mobility + 6 * pd_patients$PDQ39_adl)/16
if (any(!pd_patients$sex_male %in% c(0, 1), na.rm = TRUE)) stop("sex_male must use 0/1 coding.")

pd_pairs <- NULL
pair_file <- file.path(PD_DATA_DIR, "analysis_pairs.csv")
if (file.exists(pair_file)) {
  pd_pairs <- pd_read(pair_file)
  index <- match(pd_pairs$ID, pd_patients$ID)
  if (anyNA(index)) stop("A paired patient is absent from the patient CSV.")
  shared <- intersect(names(pd_pairs), names(pd_patients))
  shared <- setdiff(shared, c(pd_domains, "PDQ39_SI", "PDQ39_total", "PDQ39_total_0_100"))
  if (!isTRUE(all.equal(pd_pairs[shared], pd_patients[index, shared], check.attributes = FALSE))) {
    stop("Shared non-PDQ patient variables disagree between patient and pair CSVs.")
  }
  derived <- c(pd_domains, "PDQ39_SI", "PDQ39_six_domain_mean", "FMCS")
  if ("PDQ39_total" %in% names(pd_patients)) derived <- c(derived, "PDQ39_total", "PDQ39_total_0_100")
  pd_pairs[derived] <- pd_patients[index, derived]
}

pd_function <- if (!is.null(pd_items) && all(c("independent_walk", "balance_obstacle") %in% names(pd_items)))
  pd_items[c("ID", "independent_walk", "balance_obstacle")] else NULL
if (nzchar(Sys.getenv("PD_FUNCTION_CSV"))) pd_function <- pd_read(Sys.getenv("PD_FUNCTION_CSV"))
if (!is.null(pd_function)) {
  pd_require(pd_function, c("independent_walk", "balance_obstacle"))
  for (v in c("independent_walk", "balance_obstacle")) {
    if (any(!pd_function[[v]] %in% c(1, 2), na.rm = TRUE)) stop("Functional items use 1=yes, 2=no.")
  }
  index <- match(pd_patients$ID, pd_function$ID)
  pd_patients$balance_difficulty <- as.numeric(pd_function$balance_obstacle[index] == 1)
  pd_patients$unable_independent_walk_stand <- as.numeric(pd_function$independent_walk[index] == 2)
  if (!is.null(pd_pairs)) {
    index <- match(pd_pairs$ID, pd_patients$ID)
    pd_pairs[c("balance_difficulty", "unable_independent_walk_stand")] <-
      pd_patients[index, c("balance_difficulty", "unable_independent_walk_stand")]
  }
}

pd_fit_sem <- function(model, data, name, resamples) {
  point <- identical(Sys.getenv("PD_POINT_ONLY"), "1")
  archive <- Sys.getenv("PD_MODEL_ARCHIVE_DIR")
  path <- file.path(archive, paste0(name, ".rds"))
  cached <- nzchar(archive) && file.exists(path) && !point
  if (cached) {
    pd_input_hashes[path] <<- unname(tools::md5sum(path))
    fit <- readRDS(path)
    old_data <- lavaan::lavInspect(fit, "data")
    pd_require(data, colnames(old_data))
    if (nrow(data) != nrow(old_data) || max(abs(as.matrix(data[colnames(old_data)]) - old_data)) > 1e-10) {
      stop("Archived model inputs differ from current standardized inputs.")
    }
    fresh <- lavaan::sem(model, data = data, se = "none")
    a <- lavaan::parameterEstimates(fit)
    b <- lavaan::parameterEstimates(fresh)
    key <- function(x) paste(x$lhs, x$op, x$rhs)
    selected <- a$op %in% c("~", ":=", "~~")
    index <- match(key(a)[selected], key(b))
    if (anyNA(index) || max(abs(a$est[selected] - b$est[index])) > 1e-8) stop("Archived model specification differs.")
  } else {
    set.seed(20260912)
    fit <- lavaan::sem(model, data = data, se = if (point) "none" else "bootstrap",
                       bootstrap = if (point) 0L else resamples)
  }
  if (!lavaan::lavInspect(fit, "converged") || !lavaan::lavInspect(fit, "post.check")) {
    stop("SEM did not converge to an admissible solution.")
  }
  B <- if (point) NULL else lavaan::lavInspect(fit, "boot")
  if (!point && (nrow(B) != resamples || any(!is.finite(B)) ||
                length(attr(B, "error.idx")) || length(attr(B, "nonadmissible")))) {
    stop("Bootstrap draws require review; no final table was written.")
  }
  check <- data.frame(model = name, n = nrow(data), bootstrap_draws = if (point) 0L else nrow(B),
                       source = if (cached) "verified_archive" else if (point) "point_check" else "new_bootstrap")
  write.csv(check, file.path(PD_OUTPUT_DIR, "tables", paste0(name, "_bootstrap_check.csv")), row.names = FALSE)
  saveRDS(fit, file.path(PD_OUTPUT_DIR, "data_processed", paste0(name, ".rds")))
  fit
}
