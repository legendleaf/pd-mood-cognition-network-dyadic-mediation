# Score distributions, retained-record completeness and functional descriptions.
source("script/Step0_DataPreparation.R")
frames <- list(patients = pd_patients)
if (!is.null(pd_pairs)) frames$dyads <- pd_pairs
descriptives <- frequencies <- functions <- list()
for (cohort in names(frames)) {
  data <- frames[[cohort]]
  columns <- intersect(c(pd_network_nodes, "PDQ39_bodily", "PDQ39_SI", "NMS_total", "AD8_total",
    "EQ5D_mobility", "EQ5D_selfcare", "EQ5D_activity", "EQ5D_anxiety", "age", "sex_male",
    "disease_duration", "LEDD", "FMCS", "ZBI_score", "DASS_depression", "DASS_anxiety", "DASS_stress"), names(data))
  if (all(c("DASS_depression", "DASS_anxiety", "DASS_stress") %in% names(data))) {
    data$DASS_total <- data$DASS_depression + data$DASS_anxiety + data$DASS_stress
    columns <- c(columns, "DASS_total")
  }
  for (v in columns) {
    x <- data[[v]][is.finite(data[[v]])]
    if (!length(x)) next
    q <- quantile(x, c(.25, .75), names = FALSE)
    descriptives[[paste(cohort, v)]] <- data.frame(cohort = cohort, variable = v, n = length(x),
      missing = nrow(data) - length(x), mean = mean(x), sd = sd(x), median = median(x),
      q1 = q[1], q3 = q[2], min = min(x), max = max(x))
    counts <- table(x)
    frequencies[[paste(cohort, v)]] <- data.frame(cohort = cohort, variable = v,
      value = as.numeric(names(counts)), n = as.integer(counts), denominator = length(x),
      percent_nonmissing = 100*as.integer(counts)/length(x))
  }
  for (v in intersect(c("balance_difficulty", "unable_independent_walk_stand",
                         "EQ5D_mobility", "EQ5D_selfcare", "EQ5D_activity"), names(data))) {
    x <- data[[v]]
    valid <- is.finite(x)
    limitations <- if (startsWith(v, "EQ5D")) x %in% 2:5 else x == 1
    n <- sum(limitations & valid)
    functions[[paste(cohort, v)]] <- data.frame(cohort = cohort, variable = v,
      limitations_n = n, denominator = sum(valid), missing = sum(!valid), percent = 100*n/sum(valid))
  }
}
write.csv(do.call(rbind, descriptives), file.path(PD_OUTPUT_DIR, "tables", "score_descriptives.csv"), row.names = FALSE)
write.csv(do.call(rbind, frequencies), file.path(PD_OUTPUT_DIR, "tables", "score_distributions.csv"), row.names = FALSE)
write.csv(do.call(rbind, functions), file.path(PD_OUTPUT_DIR, "tables", "functional_descriptors.csv"), row.names = FALSE)
if (!is.null(pd_items)) {
  groups <- list(NMSQuest = paste0("\u975e\u8fd0\u52a8_", 1:30),
    NMS_mood = paste0("\u975e\u8fd0\u52a8_", c(12, 13, 15, 16, 17)),
    PDQ39 = paste0("PDQ39_", 1:39), WOQ9 = paste0("WOQ9_", 1:9))
  completeness <- do.call(rbind, lapply(names(groups), function(name) {
    pd_require(pd_items, groups[[name]])
    items <- as.matrix(pd_items[groups[[name]]])
    data.frame(instrument = name, participants = nrow(items), items = ncol(items),
      expected_entries = length(items), missing_entries = sum(is.na(items)))
  }))
  pd_require(pd_items, c("WOQ9", groups$WOQ9))
  woq <- as.matrix(pd_items[groups$WOQ9])
  if (anyNA(woq) || any(!woq %in% 0:3)) stop("Unexpected WOQ item coding or missing entries.")
  if (any(rowSums(woq == 1) != pd_patients$WOQ9_total) ||
      any(pd_items$WOQ9 != pd_patients$WOQ9_total)) stop("WOQ item count does not match retained totals.")
  write.csv(completeness, file.path(PD_OUTPUT_DIR, "tables", "item_completeness.csv"), row.names = FALSE)
} else {
  warning("Item workbook unavailable: item completeness and WOQ item-total checks were not performed.")
}
pd_finish("step5")
