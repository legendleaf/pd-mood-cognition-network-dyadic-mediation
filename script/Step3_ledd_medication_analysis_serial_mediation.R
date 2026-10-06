# Step 3: Exploratory cross-sectional serial model.
# The filename is retained for compatibility with the original workflow.
source("script/Step0_DataPreparation.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lavaan)
  library(ggplot2)
})

# ---- config ----

OUTPUT_DIR     <- PD_OUTPUT_DIR

set.seed(2024)
bootstrap_iter <- pd_count("PD_SERIAL_BOOTSTRAP", 5000)
export_pdf     <- TRUE
png_dpi        <- 600

for (sub in c("figures/main", "figures/supplement", "tables", "data_processed")) {
  dir.create(file.path(OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

# ---- data ----

patients <- pd_patients

required_vars <- c("LEDD", "WOQ9_total", "NMS_mood", "PDQ39_SI",
                   "age", "sex_male", "disease_duration")
missing_vars  <- setdiff(required_vars, names(patients))
if (length(missing_vars) > 0)
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))

analysis_data_raw <- pd_complete(patients, c("ID", required_vars), "serial_model") %>%
  dplyr::select(ID, all_of(required_vars))

n_complete <- nrow(analysis_data_raw)

# Z-standardize all continuous variables (sex_male stays binary)
analysis_data <- analysis_data_raw %>%
  mutate(across(c(LEDD, WOQ9_total, NMS_mood, PDQ39_SI, age, disease_duration),
                ~ as.numeric(scale(.)),
                .names = "{.col}_z"))

control_formula <- "age_z + sex_male + disease_duration_z"

# ---- multicollinearity in each regression equation ----

vif_formulas <- c(
  paste("WOQ9_total_z ~ LEDD_z +", control_formula),
  paste("NMS_mood_z ~ WOQ9_total_z + LEDD_z +", control_formula),
  paste("PDQ39_SI_z ~ WOQ9_total_z + NMS_mood_z + LEDD_z +", control_formula))
vif_table <- bind_rows(lapply(vif_formulas, function(spec) {
  predictors <- attr(terms(as.formula(spec)), "term.labels")
  values <- vapply(predictors, function(v) {
    auxiliary <- lm(reformulate(setdiff(predictors, v), v), data = analysis_data)
    1/(1-summary(auxiliary)$r.squared)
  }, numeric(1))
  data.frame(equation = sub(" ~.*", "", spec), predictor = names(values), VIF = as.numeric(values))
}))
write_csv(vif_table, file.path(OUTPUT_DIR, "tables", "serial_equation_VIF.csv"))

# ---- serial mediation SEM ----

serial_model_spec <- sprintf('
  WOQ9_total_z ~ a1*LEDD_z + %s
  NMS_mood_z   ~ a2*WOQ9_total_z + d21*LEDD_z + %s
  PDQ39_SI_z   ~ b1*WOQ9_total_z + b2*NMS_mood_z + c_prime*LEDD_z + %s

  ind1_simple    := a1 * b1            # LEDD -> WOQ9 -> PDQ39
  ind2_serial    := a1 * a2 * b2       # LEDD -> WOQ9 -> NMS_mood -> PDQ39 (target)
  ind3_direct_m2 := d21 * b2           # LEDD -> NMS_mood -> PDQ39

  total_indirect    := ind1_simple + ind2_serial + ind3_direct_m2
  total_effect      := c_prime + total_indirect
', control_formula, control_formula, control_formula)

fit_serial <- pd_fit_sem(serial_model_spec, analysis_data, "serial_corrected_si", bootstrap_iter)

params <- parameterEstimates(fit_serial, boot.ci.type = "perc", ci = TRUE)

get_path    <- function(lbl) params %>% filter(label == lbl)
path_a1     <- get_path("a1")
path_a2     <- get_path("a2")
path_d21    <- get_path("d21")
path_b1     <- get_path("b1")
path_b2     <- get_path("b2")
path_c      <- get_path("c_prime")
ind1        <- get_path("ind1_simple")
ind2        <- get_path("ind2_serial")
ind3        <- get_path("ind3_direct_m2")
total_ind   <- get_path("total_indirect")
total_eff   <- get_path("total_effect")

# ---- save results ----

results_table <- data.frame(
  Path = c(
    "a1: LEDD -> WOQ9",
    "a2: WOQ9 -> NMS_mood",
    "d21: LEDD -> NMS_mood",
    "b1: WOQ9 -> PDQ39",
    "b2: NMS_mood -> PDQ39",
    "c': LEDD -> PDQ39 (direct)",
    "Indirect 1: LEDD -> WOQ9 -> PDQ39",
    "Indirect 2: LEDD -> WOQ9 -> NMS -> PDQ39 (SERIAL)",
    "Indirect 3: LEDD -> NMS -> PDQ39",
    "Total indirect association",
    "Total association"
  ),
  Beta_Standardized = c(path_a1$est, path_a2$est, path_d21$est,
                        path_b1$est, path_b2$est, path_c$est,
                        ind1$est, ind2$est, ind3$est,
                        total_ind$est, total_eff$est),
  CI_Lower = c(path_a1$ci.lower, path_a2$ci.lower, path_d21$ci.lower,
               path_b1$ci.lower, path_b2$ci.lower, path_c$ci.lower,
               ind1$ci.lower, ind2$ci.lower, ind3$ci.lower,
               total_ind$ci.lower, total_eff$ci.lower),
  CI_Upper = c(path_a1$ci.upper, path_a2$ci.upper, path_d21$ci.upper,
               path_b1$ci.upper, path_b2$ci.upper, path_c$ci.upper,
               ind1$ci.upper, ind2$ci.upper, ind3$ci.upper,
               total_ind$ci.upper, total_eff$ci.upper),
  p_value = c(path_a1$pvalue, path_a2$pvalue, path_d21$pvalue,
              path_b1$pvalue, path_b2$pvalue, path_c$pvalue,
              ind1$pvalue, ind2$pvalue, ind3$pvalue,
              total_ind$pvalue, total_eff$pvalue),
  Significant = c(path_a1$pvalue < 0.05, path_a2$pvalue < 0.05,
                  path_d21$pvalue < 0.05, path_b1$pvalue < 0.05,
                  path_b2$pvalue < 0.05, path_c$pvalue < 0.05,
                  ind1$pvalue < 0.05, ind2$pvalue < 0.05, ind3$pvalue < 0.05,
                  total_ind$pvalue < 0.05, total_eff$pvalue < 0.05),
  N              = n_complete,
  Bootstrap_iter = bootstrap_iter
)
write_csv(results_table,
          file.path(OUTPUT_DIR, "tables", "mediation_serial_results.csv"))

# ---- Figure 3: path diagram ----

nodes_df <- data.frame(
  x     = c(1, 4, 7, 10),
  y     = c(3, 4.5, 4.5, 3),
  label = c("LEDD", "WOQ-9\nTotal", "NMSQuest\nMood/Cognition", "PDQ-39 SI"),
  color = c("#E8F4F8", "#FFF4E6", "#FFE6E6", "#E8F5E8")
)

star <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))

arrows_df <- data.frame(
  x    = c(1, 4, 7, 4, 1, 1),
  y    = c(3, 4.5, 4.5, 4.2, 3.2, 3.5),
  xend = c(4, 7, 10, 10, 7, 10),
  yend = c(4.5, 4.5, 3, 3.2, 4.3, 2.8),
  label = c(
    sprintf("a1 = %.2f%s",  path_a1$est,  star(path_a1$pvalue)),
    sprintf("a2 = %.2f%s",  path_a2$est,  star(path_a2$pvalue)),
    sprintf("b2 = %.2f%s",  path_b2$est,  star(path_b2$pvalue)),
    sprintf("b1 = %.2f%s",  path_b1$est,  star(path_b1$pvalue)),
    sprintf("d21 = %.2f%s", path_d21$est, star(path_d21$pvalue)),
    sprintf("c' = %.2f%s",  path_c$est,   star(path_c$pvalue))
  ),
  type = c("Serial", "Serial", "Serial", "Direct", "Direct M2", "Direct")
)

path_diagram <- ggplot() +
  geom_segment(data = arrows_df,
               aes(x = x, y = y, xend = xend, yend = yend,
                   color = type, linetype = type, linewidth = type),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  geom_rect(data = nodes_df,
            aes(xmin = x - 0.8, xmax = x + 0.8,
                ymin = y - 0.6, ymax = y + 0.6),
            fill = nodes_df$color, color = "black", linewidth = 1) +
  geom_text(data = nodes_df, aes(x, y, label = label),
            size = 4, fontface = "bold") +
  geom_text(data = arrows_df,
            aes(x = (x + xend) / 2, y = (y + yend) / 2 + 0.3, label = label),
            size = 3.5, fontface = "bold") +
  scale_color_manual(values = c(Serial = "#D32F2F", Direct = "#1976D2",
                                `Direct M2` = "#7B1FA2")) +
  scale_linetype_manual(values = c(Serial = 1, Direct = 2, `Direct M2` = 3)) +
  scale_linewidth_manual(values = c(Serial = 1.2, Direct = 0.8, `Direct M2` = 0.8)) +
  xlim(0, 11) + ylim(1.5, 6) +
  labs(title = "Exploratory cross-sectional serial model",
       subtitle = sprintf("Serial indirect association: beta = %.3f [%.3f, %.3f], p = %.3f",
                          ind2$est, ind2$ci.lower, ind2$ci.upper, ind2$pvalue),
       caption = "Arrows specify regression equations. * p < .05; ** p < .01; *** p < .001 (nominal Wald tests).") +
  theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
        legend.position = "bottom", legend.title = element_blank())

path_file <- file.path(OUTPUT_DIR, "figures/main",
                       if (export_pdf) "Fig3_Mediation_Path.pdf"
                       else            "Fig3_Mediation_Path.png")
ggsave(path_file, path_diagram, width = 14, height = 8,
       dpi = if (export_pdf) 300 else png_dpi,
       bg = if (export_pdf) "transparent" else "white")

# ---- Figure S2: forest plot of standardized paths ----

plot_data <- results_table %>%
  mutate(
    Path_Type = case_when(
      grepl("[Ii]ndirect|^Total association", Path) ~ "Indirect and total associations",
      TRUE ~ "Component paths"
    ),
    Path_Clean      = gsub(".*: ", "", Path),
    Significant_cat = ifelse(Significant, "p < 0.05", "p >= 0.05")
  )

forest_plot <- ggplot(plot_data, aes(reorder(Path_Clean, Beta_Standardized),
                                     Beta_Standardized)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                width = 0.2, linewidth = 0.8, color = "gray30") +
  geom_point(aes(fill = Significant_cat, size = abs(Beta_Standardized)),
             shape = 21, stroke = 1.2) +
  scale_fill_manual(values = c(`p < 0.05` = "#2166AC", `p >= 0.05` = "gray70"),
                    name = "Significance") +
  scale_size_continuous(range = c(4, 9), guide = "none") +
  geom_text(aes(label = sprintf("p = %.3f", p_value)),
            hjust = -0.15, size = 3, color = "gray20") +
  coord_flip() +
  facet_wrap(~ Path_Type, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = "Standardized beta (95% CI)",
       title    = "Cross-sectional serial model: LEDD -> WOQ-9 -> Mood/Cognition -> SI",
       subtitle = sprintf("N = %d, Bootstrap = %d", n_complete, bootstrap_iter)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle      = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.title         = element_text(face = "bold", size = 11),
    axis.text          = element_text(size = 10, color = "black"),
    legend.position    = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold", size = 11),
    panel.border       = element_rect(fill = NA, color = "gray20", linewidth = 0.5)
  )

forest_file <- file.path(OUTPUT_DIR, "figures/supplement",
                         if (export_pdf) "FigS2_Mediation_Forest.pdf"
                         else            "FigS2_Mediation_Forest.png")
ggsave(forest_file, forest_plot, width = 12, height = 10,
       dpi = if (export_pdf) 300 else png_dpi,
       bg = if (export_pdf) "transparent" else "white")
write_csv(params, file.path(OUTPUT_DIR, "tables", "serial_all_parameters.csv"))
pd_finish("step3")
