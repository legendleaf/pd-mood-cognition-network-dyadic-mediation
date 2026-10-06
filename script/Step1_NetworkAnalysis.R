# Step 1: 16-node EBICglasso symptom network
source("script/Step0_DataPreparation.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(qgraph)
  library(bootnet)
  library(networktools)
})

# ---- config ----

OUTPUT_DIR     <- PD_OUTPUT_DIR

set.seed(2024)
boot_iter   <- pd_count("PD_NETWORK_BOOTSTRAP", 1000)
case_iter   <- pd_count("PD_CASE_BOOTSTRAP", 500)
ebic_gamma  <- 0.5
cor_method  <- "cor_auto"
export_pdf  <- TRUE
png_dpi     <- 600

for (sub in c("figures/main", "figures/supplement", "tables", "data_processed")) {
  dir.create(file.path(OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

# ---- data ----

patients <- pd_patients

network_variables <- c(
  "PDQ39_mobility", "PDQ39_adl", "PDQ39_emotional", "PDQ39_cognition",
  "PDQ39_social", "PDQ39_communication", "PDQ39_stigma",
  "NMS_sleep", "NMS_mood", "NMS_gi", "NMS_urinary",
  "NMS_cardiovascular", "NMS_perception", "NMS_sexual", "NMS_misc",
  "WOQ9_total"
)

network_data <- pd_complete(patients, c("ID", network_variables), "network_16") %>%
  dplyr::select(ID, all_of(network_variables))

n_complete <- nrow(network_data)

# ---- network estimation ----

network_matrix <- network_data %>% dplyr::select(-ID) %>% as.data.frame()

net_result <- estimateNetwork(
  network_matrix,
  default = "EBICglasso", tuning = ebic_gamma,
  corMethod = cor_method, missing = "pairwise"
)

adj_matrix <- net_result$graph
n_edges    <- sum(adj_matrix[upper.tri(adj_matrix)] != 0)
sparsity   <- 1 - n_edges / choose(length(network_variables), 2)

# ---- centrality and bridge ----

centrality_table <- centralityTable(net_result, standardized = FALSE)
write_csv(centrality_table, file.path(OUTPUT_DIR, "tables", "centrality_indices.csv"))

# Three communities: PDQ-39 / NMS / WOQ-9
communities <- rep(NA_integer_, length(network_variables))
communities[grepl("PDQ39", network_variables)] <- 1
communities[grepl("NMS",   network_variables)] <- 2
communities[grepl("WOQ",   network_variables)] <- 3

bridge_result <- bridge(adj_matrix, communities = communities)
bridge_df <- data.frame(
  Node            = names(bridge_result$`Bridge Strength`),
  Bridge_Strength = as.numeric(bridge_result$`Bridge Strength`),
  Community       = c("PDQ-39", "NMS", "WOQ-9")[communities]
)
write_csv(bridge_df, file.path(OUTPUT_DIR, "tables", "bridge_centrality.csv"))

# ---- bootstrap stability ----

archive_path <- Sys.getenv("PD_NETWORK_ARCHIVE")
if (nzchar(archive_path)) {
  archived <- new.env()
  pd_input_hashes[archive_path] <- unname(tools::md5sum(archive_path))
  load(archive_path, envir = archived)
  old_matrix <- as.matrix(archived$net_result$data)
  centered_current <- sweep(as.matrix(network_matrix), 2, colMeans(network_matrix))
  centered_old <- sweep(old_matrix, 2, colMeans(old_matrix))
  stopifnot(identical(network_variables, archived$network_variables),
    nrow(network_matrix) == nrow(archived$net_result$data),
    max(abs(centered_current - centered_old)) < 1e-10,
    max(abs(as.matrix(cor(network_matrix)) - cor(archived$net_result$data))) < 1e-10,
    max(abs(adj_matrix - archived$adj_matrix)) < 1e-6,
    length(archived$boot_result$boots) == boot_iter,
    length(archived$boot_case$boots) == case_iter)
  boot_result <- archived$boot_result
  boot_case <- archived$boot_case
} else {
  set.seed(2024)
  boot_result <- bootnet(net_result, nBoots = boot_iter, type = "nonparametric",
                         statistics = c("edge", "strength", "closeness", "betweenness"))
  boot_case <- bootnet(net_result, nBoots = case_iter, type = "case",
                       statistics = c("strength", "closeness", "betweenness"))
}

cs_coef <- corStability(boot_case)["strength"]

# ---- Figure 1A: network graph ----

official_labels_map <- c(
  PDQ39_mobility      = "Mobility",
  PDQ39_adl           = "ADL",
  PDQ39_emotional     = "Emotional",
  PDQ39_cognition     = "Cognition",
  PDQ39_social        = "Social",
  PDQ39_communication = "Communication",
  PDQ39_stigma        = "Stigma",
  NMS_sleep           = "Sleep-Fatigue",
  NMS_mood            = "Mood-Cognition",
  NMS_gi              = "Gastrointestinal",
  NMS_urinary         = "Urinary",
  NMS_cardiovascular  = "Cardiovascular",
  NMS_perception      = "Perceptual-Hallucinations",
  NMS_sexual          = "Sexual",
  NMS_misc            = "Pain-Misc",
  WOQ9_total          = "WOQ-9"
)

adj <- as.matrix(net_result$graph); diag(adj) <- 0
nodes <- colnames(adj)
final_labels <- official_labels_map[nodes]

ct <- centralityTable(net_result) %>% filter(measure == "Strength")
strength_vec <- setNames(ct$value, ct$node)[nodes]
strength_vec[is.na(strength_vec)] <- 0
GLOBAL_MAX_STRENGTH <- max(strength_vec, 1)
vsize_vec <- 8 + (strength_vec / GLOBAL_MAX_STRENGTH) * 6

group_vec <- case_when(
  grepl("PDQ39", nodes) ~ "PDQ-39",
  grepl("NMS",   nodes) ~ "NMS",
  grepl("WOQ",   nodes) ~ "WOQ-9"
)
groups_list   <- list(`PDQ-39` = which(group_vec == "PDQ-39"),
                      NMS      = which(group_vec == "NMS"),
                      `WOQ-9`  = which(group_vec == "WOQ-9"))
custom_colors <- c(`PDQ-39` = "#D0CADE", NMS = "#A7D2BA", `WOQ-9` = "#F7C6A8")
color_vec     <- unname(custom_colors[names(groups_list)])

set.seed(2024)
Layout_Fixed <- qgraph(adj, layout = "spring", DoNotPlot = TRUE)$layout

fig1a_file <- file.path(OUTPUT_DIR, "figures/main",
                        if (export_pdf) "Fig1A_Network_16nodes.pdf"
                        else            "Fig1A_Network_16nodes.png")
if (export_pdf) {
  cairo_pdf(fig1a_file, width = 10, height = 9, family = "Arial")
} else {
  png(fig1a_file, width = 10, height = 9, units = "in", res = png_dpi, bg = "white")
}

qgraph(
  adj, layout = Layout_Fixed,
  vsize = vsize_vec, groups = groups_list, color = color_vec, borders = FALSE,
  labels = final_labels, label.cex = 1.1, label.scale = FALSE,
  label.color = "black", label.font = 2,
  minimum = 0.05, cut = 0.1, details = FALSE,
  edge.width = 1.5, posCol = "#ABC6E4", negCol = "#C39398",
  fade = FALSE, curveAll = TRUE, curveScale = TRUE,
  legend = FALSE, mar = c(3, 3, 3, 3)
)
legend("topright",
       legend = c("PDQ-39", "NMS", "WOQ-9"),
       col = c("#D0CADE", "#A7D2BA", "#F7C6A8"),
       pch = 19, pt.cex = 2, bty = "n", cex = 0.9,
       text.font = 2, y.intersp = 1.02, inset = c(0.03, 0.03))
dev.off()

# ---- Figures 1B and 1C: centrality and bridge bar plots ----

get_group <- function(name) {
  if (name %in% c("Mobility", "ADL", "Emotional", "Cognition",
                  "Social", "Communication", "Stigma")) return("PDQ-39")
  if (name == "WOQ-9") return("WOQ-9")
  "NMS"
}

bar_theme <- theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    axis.text.y        = element_text(color = "black", size = 9, face = "bold"),
    axis.text.x        = element_text(color = "black", size = 10),
    axis.title.x       = element_text(face = "bold", margin = margin(t = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line.x        = element_line(color = "black", linewidth = 0.5),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 10)
  )

cent_df <- centrality_table %>%
  filter(measure == "Strength") %>%
  mutate(Official_Name = official_labels_map[node],
         Group         = sapply(Official_Name, get_group)) %>%
  arrange(value) %>%
  mutate(Official_Name = factor(Official_Name, levels = Official_Name))

p1b <- ggplot(cent_df, aes(Official_Name, value, fill = Group)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  labs(x = NULL, y = "Strength Centrality") + bar_theme

bridge_plot_df <- bridge_df %>%
  mutate(Official_Name = official_labels_map[Node],
         Community     = factor(Community, levels = c("PDQ-39", "NMS", "WOQ-9"))) %>%
  arrange(Bridge_Strength) %>%
  mutate(Official_Name = factor(Official_Name, levels = Official_Name))

p1c <- ggplot(bridge_plot_df, aes(Official_Name, Bridge_Strength, fill = Community)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  labs(x = NULL, y = "Bridge Strength") + bar_theme

save_fig <- function(plot, name, w, h) {
  out <- file.path(OUTPUT_DIR, "figures/main",
                   paste0(name, if (export_pdf) ".pdf" else ".png"))
  if (export_pdf) {
    ggsave(out, plot, width = w, height = h, device = cairo_pdf)
  } else {
    ggsave(out, plot, width = w, height = h, dpi = png_dpi, bg = "white")
  }
}
save_fig(p1b, "Fig1B_Centrality_16nodes", 6, 6)
save_fig(p1c, "Fig1C_Bridge_16nodes",     6, 6)

# ---- Figure S1: case-dropping stability ----

figS1_file <- file.path(OUTPUT_DIR, "figures/supplement",
                        if (export_pdf) "FigS1_Stability.pdf" else "FigS1_Stability.png")
if (export_pdf) {
  cairo_pdf(figS1_file, width = 10, height = 6, family = "Arial")
} else {
  png(figS1_file, width = 10, height = 6, units = "in", res = png_dpi, bg = "white")
}

print(
  plot(boot_case, statistics = "strength") +
    theme_minimal(base_family = "Arial", base_size = 14) +
    labs(title    = "Network Stability Analysis",
         subtitle = sprintf("Case-dropping bootstrap | CS-coefficient = %.3f", cs_coef),
         x = "Proportion of Cases Retained",
         y = "Average Correlation with Original Centrality")
)
dev.off()

# ---- save objects for Step 2 ----

save(net_result, GLOBAL_MAX_STRENGTH, network_variables, network_data,
     adj_matrix, official_labels_map, cs_coef, boot_result, boot_case,
     file = file.path(OUTPUT_DIR, "data_processed", "step1_results.RData"))
write.csv(data.frame(n = n_complete, nodes = length(network_variables), edges = n_edges,
  sparsity = sparsity, strength_CS = unname(cs_coef)),
  file.path(OUTPUT_DIR, "tables", "network_summary.csv"), row.names = FALSE)
pd_finish("step1")
