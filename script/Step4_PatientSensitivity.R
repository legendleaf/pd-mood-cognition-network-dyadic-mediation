# Full-cohort patient HRQoL regressions (Supplementary Tables S11A/S11B).
source("script/Step0_DataPreparation.R")
if (!requireNamespace("sandwich", quietly = TRUE)) stop("Install sandwich before running Step 4.")
predictors <- c("NMS_mood", "age", "sex_male", "disease_duration")
outcomes <- c("PDQ39_SI", "PDQ39_six_domain_mean")
tables <- summaries <- vifs <- list()
for (outcome in outcomes) {
  data <- pd_complete(pd_patients, c(outcome, predictors), paste0("patient_", outcome))
  fit <- lm(reformulate(predictors, outcome), data = data, na.action = na.fail)
  b <- coef(fit)
  se <- sqrt(diag(sandwich::vcovHC(fit, type = "HC3")))
  df <- df.residual(fit)
  scale_factor <- c(NA_real_, vapply(data[predictors], sd, numeric(1))/sd(data[[outcome]]))
  tables[[outcome]] <- data.frame(outcome = outcome, term = names(b), n = nobs(fit), df = df,
    estimate = unname(b), se_HC3 = unname(se), lower = unname(b - qt(.975, df)*se),
    upper = unname(b + qt(.975, df)*se), p = unname(2*pt(-abs(b/se), df)),
    beta = unname(b)*scale_factor)
  summaries[[outcome]] <- data.frame(outcome = outcome, n = nobs(fit), residual_df = df,
    R2 = summary(fit)$r.squared, adjusted_R2 = summary(fit)$adj.r.squared)
  vifs[[outcome]] <- do.call(rbind, lapply(predictors, function(v) {
    auxiliary <- lm(reformulate(setdiff(predictors, v), v), data = data)
    data.frame(outcome = outcome, predictor = v, VIF = 1/(1-summary(auxiliary)$r.squared))
  }))
}
result <- do.call(rbind, tables)
mood <- result[result$term == "NMS_mood", ]
mood$p_Holm_two_outcomes <- p.adjust(mood$p, method = "holm")
write.csv(result, file.path(PD_OUTPUT_DIR, "tables", "patient_regressions_HC3.csv"), row.names = FALSE)
write.csv(mood, file.path(PD_OUTPUT_DIR, "tables", "patient_mood_associations.csv"), row.names = FALSE)
write.csv(do.call(rbind, summaries), file.path(PD_OUTPUT_DIR, "tables", "patient_model_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, vifs), file.path(PD_OUTPUT_DIR, "tables", "patient_regression_VIF.csv"), row.names = FALSE)
pd_finish("step4")
