# Step 2: dyadic and secondary analyses
# Part 1 - NCT for SES moderation on a 10-node core network
# Part 2 - patient-centered actor-partner model (NMS_mood -> ZBI / PDQ39_SI)
# Part 3 - caregiver-patient discrepancy (incremental, directional, residual)

suppressPackageStartupMessages({
  library(tidyverse)
  library(bootnet)
  library(NetworkComparisonTest)
  library(lavaan)
  library(boot)
  library(patchwork)
})

# ---- config ----

INPUT_PATIENTS <- "analysis_patients.csv"
INPUT_PAIRS    <- "analysis_pairs.csv"
STEP1_RESULTS  <- "results_network_analysis/data_processed/step1_results.RData"
OUTPUT_DIR     <- "results_network_analysis"

set.seed(2024)
nct_iter           <- 1000
nct_gamma          <- 0.5
nct_cor_method     <- "cor_auto"
sem_bootstrap      <- 2000
mismatch_bootstrap <- 2000
export_pdf         <- TRUE
png_dpi            <- 600

for (sub in c("figures/main", "figures/supplement", "tables", "data_processed")) {
  dir.create(file.path(OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

if (file.exists(STEP1_RESULTS)) load(STEP1_RESULTS)
patients <- read_csv(INPUT_PATIENTS, show_col_types = FALSE)
pairs    <- read_csv(INPUT_PAIRS,    show_col_types = FALSE)

# ---- 10-node core for NCT ----

core_10_nodes <- c(
  "NMS_mood", "NMS_sleep", "NMS_cardiovascular", "NMS_misc",
  "PDQ39_emotional", "PDQ39_cognition", "PDQ39_mobility",
  "PDQ39_adl", "PDQ39_communication",
  "WOQ9_total"
)

# ---- SES composite index ----

if (!"SES_Index" %in% names(patients)) {
  ses_indicators <- c("education_level", "pension_category", "income_category",
                      "has_insurance", "economic_pressure")
  ses_available  <- intersect(ses_indicators, names(patients))
  if (length(ses_available) < 3) stop("Need at least 3 SES indicators")

  ses_data <- patients %>%
    dplyr::select(ID, all_of(ses_available)) %>%
    na.omit()

  # Reverse code: higher economic_pressure means worse SES
  if ("economic_pressure" %in% names(ses_data)) {
    max_p <- max(ses_data$economic_pressure, na.rm = TRUE)
    ses_data <- ses_data %>%
      mutate(economic_pressure = (max_p + 1) - economic_pressure)
  }

  ses_data$SES_Index <- rowMeans(scale(ses_data %>% dplyr::select(-ID)), na.rm = TRUE)
  patients <- patients %>%
    left_join(ses_data %>% dplyr::select(ID, SES_Index), by = "ID")
}

patients <- patients %>%
  mutate(SES_Group = if_else(SES_Index >= median(SES_Index, na.rm = TRUE),
                             "High", "Low"))

# =========================================================================
# Part 1: NCT for SES moderation
# =========================================================================

nct_data <- patients %>%
  dplyr::select(ID, SES_Group, all_of(core_10_nodes)) %>%
  na.omit()

n_high <- sum(nct_data$SES_Group == "High")
n_low  <- sum(nct_data$SES_Group == "Low")

data_high <- nct_data %>% filter(SES_Group == "High") %>%
  dplyr::select(all_of(core_10_nodes)) %>% as.data.frame()
data_low  <- nct_data %>% filter(SES_Group == "Low")  %>%
  dplyr::select(all_of(core_10_nodes)) %>% as.data.frame()

nct_result <- NetworkComparisonTest::NCT(
  data_high, data_low,
  it = nct_iter,
  binary.data = FALSE, paired = FALSE, weighted = TRUE,
  test.edges = TRUE, edges = "all",
  progressbar = TRUE, make.positive.definite = TRUE,
  p.adjust.methods = "none",
  test.centrality = TRUE,
  centrality = c("strength", "betweenness", "closeness")
)

# Re-estimate group networks for global strength comparison and plotting
net_high <- estimateNetwork(data_high, default = "EBICglasso",
                            tuning = nct_gamma, corMethod = nct_cor_method)
net_low  <- estimateNetwork(data_low,  default = "EBICglasso",
                            tuning = nct_gamma, corMethod = nct_cor_method)
gs_high <- sum(abs(net_high$graph[upper.tri(net_high$graph)]))
gs_low  <- sum(abs(net_low$graph[upper.tri(net_low$graph)]))

edge_diffs_fdr <- if (!is.null(nct_result$einv.pvals))
  p.adjust(nct_result$einv.pvals$`p-value`, method = "fdr") else NULL

nct_summary <- tibble(
  Test       = c("Global Strength", "Edge Differences (FDR)"),
  High_SES   = c(gs_high, NA),
  Low_SES    = c(gs_low,  NA),
  Difference = c(gs_high - gs_low, NA),
  p_value    = c(nct_result$glstrinv.pval,
                 if (!is.null(edge_diffs_fdr)) sum(edge_diffs_fdr < 0.05) else NA)
)
write_csv(nct_summary, file.path(OUTPUT_DIR, "tables", "nct_ses_summary.csv"))
save(nct_result, file = file.path(OUTPUT_DIR, "data_processed",
                                  "nct_ses_moderation.RData"))

# =========================================================================
# Part 2: Patient-centered APIM (NMS_mood -> PDQ39_SI / ZBI)
# =========================================================================

sem_data <- pairs %>%
  dplyr::select(any_of(c("ID", "NMS_mood", "PDQ39_SI",
                         "ZBI_score", "disease_duration"))) %>%
  na.omit() %>%
  mutate(across(where(is.numeric) & !matches("^ID$"),
                ~ as.numeric(scale(.)),
                .names = "{.col}_z"))

apim_model <- '
  PDQ39_SI_z  ~ a1*NMS_mood_z + disease_duration_z   # Actor:   patient mood -> patient HRQoL
  ZBI_score_z ~ p1*NMS_mood_z + disease_duration_z   # Partner: patient mood -> caregiver burden
'

sem_apim <- sem(apim_model, data = sem_data,
                se = "bootstrap", bootstrap = sem_bootstrap, verbose = FALSE)

apim_params <- parameterEstimates(sem_apim, boot.ci.type = "perc") %>%
  filter(op == "~", lhs %in% c("PDQ39_SI_z", "ZBI_score_z")) %>%
  mutate(Effect_Type = if_else(lhs == "PDQ39_SI_z", "Actor", "Partner")) %>%
  dplyr::select(Effect_Type, Predictor = rhs, Beta = est, SE = se,
                CI_Lower = ci.lower, CI_Upper = ci.upper, p = pvalue)

apim_results <- apim_params %>%
  filter(Predictor == "NMS_mood_z") %>%
  mutate(Model = "APIM (NMS_mood -> ZBI, control: duration)", .before = 1)
write_csv(apim_results, file.path(OUTPUT_DIR, "tables", "apim_sem_results.csv"))

# =========================================================================
# Part 3: Caregiver-patient discrepancy
# =========================================================================

mismatch_data_z <- pairs %>%
  dplyr::select(any_of(c("ID", "ZBI_score",
                         "DASS_depression", "DASS_anxiety", "DASS_stress",
                         "PDQ39_communication", "PDQ39_emotional", "NMS_mood"))) %>%
  na.omit() %>%
  mutate(
    DASS_total       = DASS_depression + DASS_anxiety + DASS_stress,
    PDQ_psychosocial = (PDQ39_communication + PDQ39_emotional) / 2,
    ZBI_score_z        = as.numeric(scale(ZBI_score)),
    DASS_total_z       = as.numeric(scale(DASS_total)),
    PDQ_psychosocial_z = as.numeric(scale(PDQ_psychosocial)),
    NMS_mood_z         = as.numeric(scale(NMS_mood)),
    # Interaction-salient pathway: PDQ-39 communication+emotional vs DASS
    mismatch_psysoc = abs(DASS_total_z - PDQ_psychosocial_z),
    signed_psysoc   =     DASS_total_z - PDQ_psychosocial_z,
    # Internal-symptom pathway: NMS mood/cognition vs DASS
    mismatch_biosoc = abs(DASS_total_z - NMS_mood_z),
    signed_biosoc   =     DASS_total_z - NMS_mood_z
  )

# Hierarchical regression: incremental ΔR² of mismatch beyond main effects
incremental_test <- function(data, patient_var, caregiver_var, mismatch_var,
                             outcome_var = "ZBI_score_z",
                             n_boot = mismatch_bootstrap) {
  df <- data %>%
    dplyr::select(all_of(c(outcome_var, patient_var, caregiver_var, mismatch_var))) %>%
    na.omit()

  m1 <- lm(reformulate(c(patient_var, caregiver_var), outcome_var), data = df)
  m2 <- lm(reformulate(c(patient_var, caregiver_var, mismatch_var), outcome_var),
           data = df)

  anova_res   <- anova(m1, m2)
  delta_r2    <- summary(m2)$r.squared - summary(m1)$r.squared
  mismatch_co <- summary(m2)$coefficients[mismatch_var, ]

  boot_res <- boot(df, function(d, i) {
    coef(lm(reformulate(c(patient_var, caregiver_var, mismatch_var), outcome_var),
            data = d[i, ]))[mismatch_var]
  }, R = n_boot)
  boot_ci <- boot.ci(boot_res, type = "perc")

  tibble(
    N                  = nrow(df),
    R2_Model1          = summary(m1)$r.squared,
    R2_Model2          = summary(m2)$r.squared,
    Delta_R2           = delta_r2,
    F_statistic        = anova_res$F[2],
    p_Delta_R2         = anova_res$`Pr(>F)`[2],
    Beta_Mismatch      = mismatch_co["Estimate"],
    SE_Mismatch        = mismatch_co["Std. Error"],
    t_Mismatch         = mismatch_co["t value"],
    p_Mismatch         = mismatch_co["Pr(>|t|)"],
    CI_Lower_Bootstrap = boot_ci$percent[4],
    CI_Upper_Bootstrap = boot_ci$percent[5]
  )
}

incr_psysoc <- incremental_test(mismatch_data_z, "PDQ_psychosocial_z",
                                "DASS_total_z", "mismatch_psysoc") %>%
  mutate(Pathway = "Psychosocial", .before = 1)
incr_biosoc <- incremental_test(mismatch_data_z, "NMS_mood_z",
                                "DASS_total_z", "mismatch_biosoc") %>%
  mutate(Pathway = "Biological", .before = 1)

write_csv(bind_rows(incr_psysoc, incr_biosoc),
          file.path(OUTPUT_DIR, "tables", "mismatch_incremental_validity.csv"))

# Directional: signed difference vs caregiver burden (Spearman ρ)
directional_test <- function(data, signed_var, outcome_var = "ZBI_score_z") {
  df <- data %>%
    dplyr::select(all_of(c(outcome_var, signed_var))) %>%
    na.omit()
  ct <- cor.test(df[[signed_var]], df[[outcome_var]], method = "spearman")
  tibble(N = nrow(df), Spearman_r = ct$estimate, p_value = ct$p.value)
}

dir_psysoc <- directional_test(mismatch_data_z, "signed_psysoc") %>%
  mutate(Pathway = "Psychosocial", .before = 1)
dir_biosoc <- directional_test(mismatch_data_z, "signed_biosoc") %>%
  mutate(Pathway = "Biological",   .before = 1)

write_csv(bind_rows(dir_psysoc, dir_biosoc),
          file.path(OUTPUT_DIR, "tables", "mismatch_directional_effects.csv"))

# Sensitivity: residualized mismatch (caregiver distress regressed on patient symptom)
residual_mismatch <- function(data, patient_var, caregiver_var) {
  df <- data %>% dplyr::select(all_of(c(patient_var, caregiver_var))) %>% na.omit()
  abs(residuals(lm(reformulate(patient_var, caregiver_var), data = df)))
}

mismatch_data_z <- mismatch_data_z %>%
  mutate(
    mismatch_psysoc_resid = residual_mismatch(., "PDQ_psychosocial_z", "DASS_total_z"),
    mismatch_biosoc_resid = residual_mismatch(., "NMS_mood_z",         "DASS_total_z")
  )

resid_psysoc <- incremental_test(mismatch_data_z, "PDQ_psychosocial_z",
                                 "DASS_total_z", "mismatch_psysoc_resid",
                                 n_boot = 1000) %>%
  mutate(Pathway = "Psychosocial (Residual)", .before = 1)
resid_biosoc <- incremental_test(mismatch_data_z, "NMS_mood_z",
                                 "DASS_total_z", "mismatch_biosoc_resid",
                                 n_boot = 1000) %>%
  mutate(Pathway = "Biological (Residual)", .before = 1)

write_csv(bind_rows(resid_psysoc, resid_biosoc),
          file.path(OUTPUT_DIR, "tables", "mismatch_sensitivity.csv"))

# =========================================================================
# Figures
# =========================================================================

# ---- Figure 2: mismatch scatter plots ----

theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 13, hjust = 0),
      axis.title   = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )
}

mismatch_scatter <- function(xvar, yvar, title, xlab, ylab, dr2, p) {
  ggplot(mismatch_data_z, aes(.data[[xvar]], .data[[yvar]])) +
    geom_point(alpha = 0.6, size = 3, color = "#2E86AB") +
    geom_smooth(method = "lm", color = "#A23B72", fill = "#A23B72",
                alpha = 0.2, linewidth = 1.2) +
    labs(title = title, x = xlab, y = ylab) +
    theme_pub() +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Delta R^2 = %.3f, p = %.3f", dr2, p),
             hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic")
}

fig_psysoc <- mismatch_scatter("mismatch_psysoc", "ZBI_score_z",
                               "A. Psychosocial pathway",
                               "Dyadic mismatch (absolute difference)",
                               "Caregiver burden (z-score)",
                               incr_psysoc$Delta_R2, incr_psysoc$p_Delta_R2)
fig_biosoc <- mismatch_scatter("mismatch_biosoc", "ZBI_score_z",
                               "B. Biological pathway",
                               "Dyadic mismatch (absolute difference)",
                               NULL,
                               incr_biosoc$Delta_R2, incr_biosoc$p_Delta_R2)

fig_mismatch <- fig_psysoc + fig_biosoc
ggsave(file.path(OUTPUT_DIR, "figures/main", "Figure2_Dyadic_Mismatch.png"),
       fig_mismatch, width = 10, height = 4.5, dpi = png_dpi, bg = "white")
if (export_pdf) {
  ggsave(file.path(OUTPUT_DIR, "figures/main", "Figure2_Dyadic_Mismatch.pdf"),
         fig_mismatch, width = 10, height = 4.5, dpi = 300)
}
