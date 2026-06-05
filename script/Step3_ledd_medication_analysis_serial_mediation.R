# Step 3: Serial mediation - LEDD -> WOQ-9 -> NMS_mood -> PDQ39_SI

suppressPackageStartupMessages({
  library(tidyverse)
  library(lavaan)
  library(car)
  library(ggplot2)
})

# ---- config ----

INPUT_PATIENTS <- "analysis_patients.csv"
OUTPUT_DIR     <- "results_network_analysis"

set.seed(2024)
bootstrap_iter <- 5000
export_pdf     <- TRUE
png_dpi        <- 600

for (sub in c("figures/main", "figures/supplement", "tables", "data_processed")) {
  dir.create(file.path(OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

# ---- data ----

patients <- read_csv(INPUT_PATIENTS, show_col_types = FALSE)

required_vars <- c("LEDD", "WOQ9_total", "NMS_mood", "PDQ39_SI",
                   "age", "sex_male", "disease_duration")
missing_vars  <- setdiff(required_vars, names(patients))
if (length(missing_vars) > 0)
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))

analysis_data_raw <- patients %>%
  dplyr::select(ID, all_of(required_vars)) %>%
  na.omit()

n_complete <- nrow(analysis_data_raw)

# Z-standardize all continuous variables (sex_male stays binary)
analysis_data <- analysis_data_raw %>%
  mutate(across(c(LEDD, WOQ9_total, NMS_mood, PDQ39_SI, age, disease_duration),
                ~ as.numeric(scale(.)),
                .names = "{.col}_z"))

control_formula <- "age_z + sex_male + disease_duration_z"

# ---- multicollinearity check (manuscript reports all VIF < 2.5) ----

vif_model  <- lm(as.formula(paste("PDQ39_SI_z ~ LEDD_z + WOQ9_total_z + NMS_mood_z +",
                                  control_formula)),
                 data = analysis_data)
vif_values <- vif(vif_model)

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
  prop_mediated     := total_indirect / total_effect
  serial_proportion := ind2_serial / total_indirect
', control_formula, control_formula, control_formula)

fit_serial <- sem(serial_model_spec, data = analysis_data,
                  se = "bootstrap", bootstrap = bootstrap_iter)

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
prop_med    <- get_path("prop_mediated")
serial_prop <- get_path("serial_proportion")

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
    "Total Indirect",
    "Total Effect",
    "Proportion Mediated",
    "Serial Proportion (of indirect)"
  ),
  Beta_Standardized = c(path_a1$est, path_a2$est, path_d21$est,
                        path_b1$est, path_b2$est, path_c$est,
                        ind1$est, ind2$est, ind3$est,
                        total_ind$est, total_eff$est,
                        prop_med$est, serial_prop$est),
  CI_Lower = c(path_a1$ci.lower, path_a2$ci.lower, path_d21$ci.lower,
               path_b1$ci.lower, path_b2$ci.lower, path_c$ci.lower,
               ind1$ci.lower, ind2$ci.lower, ind3$ci.lower,
               total_ind$ci.lower, total_eff$ci.lower,
               prop_med$ci.lower, serial_prop$ci.lower),
  CI_Upper = c(path_a1$ci.upper, path_a2$ci.upper, path_d21$ci.upper,
               path_b1$ci.upper, path_b2$ci.upper, path_c$ci.upper,
               ind1$ci.upper, ind2$ci.upper, ind3$ci.upper,
               total_ind$ci.upper, total_eff$ci.upper,
               prop_med$ci.upper, serial_prop$ci.upper),
  p_value = c(path_a1$pvalue, path_a2$pvalue, path_d21$pvalue,
              path_b1$pvalue, path_b2$pvalue, path_c$pvalue,
              ind1$pvalue, ind2$pvalue, ind3$pvalue,
              total_ind$pvalue, total_eff$pvalue, NA, NA),
  Significant = c(path_a1$pvalue < 0.05, path_a2$pvalue < 0.05,
                  path_d21$pvalue < 0.05, path_b1$pvalue < 0.05,
                  path_b2$pvalue < 0.05, path_c$pvalue < 0.05,
                  ind1$pvalue < 0.05, ind2$pvalue < 0.05, ind3$pvalue < 0.05,
                  total_ind$pvalue < 0.05, total_eff$pvalue < 0.05, NA, NA),
  N              = n_complete,
  Bootstrap_iter = bootstrap_iter
)
write_csv(results_table,
          file.path(OUTPUT_DIR, "tables", "mediation_serial_results.csv"))

# ---- Figure 3: path diagram ----

nodes_df <- data.frame(
  x     = c(1, 4, 7, 10),
  y     = c(3, 4.5, 4.5, 3),
  label = c("LEDD\n(Medication)", "WOQ9\n(Motor\nFluctuations)",
            "NMS\n(Mood-Cog\nHub)", "PDQ39\n(Quality\nof Life)"),
  color = c("#E8F4F8", "#FFF4E6", "#FFE6E6", "#E8F5E8")
)

star <- function(p) ifelse(p < 0.05, "***", "")

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
  labs(title    = "Serial Mediation Pathway",
       subtitle = sprintf("Indirect serial effect: beta = %.3f [%.3f, %.3f], p = %.3f",
                          ind2$est, ind2$ci.lower, ind2$ci.upper, ind2$pvalue)) +
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
  filter(!grepl("Proportion|Total Effect", Path)) %>%
  mutate(
    Path_Type = case_when(
      grepl("Indirect", Path) ~ "Indirect Effects",
      grepl("c'",       Path) ~ "Direct Effect",
      TRUE                    ~ "Component Paths"
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
       title    = "Serial mediation: LEDD -> WOQ9 -> NMS_mood -> PDQ39",
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
