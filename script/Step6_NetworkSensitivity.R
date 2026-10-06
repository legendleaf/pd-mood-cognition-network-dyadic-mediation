# Network intervals and removal of the PDQ emotional/cognition nodes (S5-S7).
source("script/Step0_DataPreparation.R")
if (!requireNamespace("qgraph", quietly = TRUE)) stop("Install qgraph before running Step 6.")
archive_path <- file.path(PD_OUTPUT_DIR, "data_processed", "step1_results.RData")
if (!file.exists(archive_path)) stop("Run Step 1 first.")
archived <- new.env()
load(archive_path, envir = archived)
estimate <- function(data) {
  R <- qgraph::cor_auto(data, missing = "pairwise", verbose = FALSE)
  qgraph::EBICglasso(R, n = nrow(data), gamma = .5, nlambda = 100,
    lambda.min.ratio = .01, refit = FALSE, threshold = FALSE, checkPD = TRUE, verbose = FALSE)
}
metrics <- function(graph) {
  if (any(!is.finite(graph)) || max(abs(graph - t(graph))) > 1e-10) stop("Invalid graph.")
  diag(graph) <- 0
  nodes <- colnames(graph)
  community <- ifelse(startsWith(nodes, "PDQ39"), "PDQ39",
    ifelse(startsWith(nodes, "NMS"), "NMS", "WOQ9"))
  cbind(strength = rowSums(abs(graph)),
    bridge_strength = rowSums(abs(graph)*outer(community, community, "!=")))
}
summarize_graphs <- function(graph, draws, prefix) {
  point <- metrics(graph)
  intervals <- differences <- list()
  for (metric in colnames(point)) {
    matrix <- t(vapply(draws, function(x) x[, metric], point[, metric]))
    ci <- apply(matrix, 2, quantile, c(.025, .975), type = 6)
    intervals[[metric]] <- data.frame(metric = metric, node = rownames(point),
      estimate = point[, metric], lower = ci[1, ], upper = ci[2, ],
      rank = rank(-point[, metric], ties.method = "min"), n_bootstrap = nrow(matrix))
    differences[[metric]] <- do.call(rbind, lapply(setdiff(rownames(point), "NMS_mood"), function(other) {
      delta <- matrix[, "NMS_mood"] - matrix[, other]
      ci <- quantile(delta, c(.025, .975), type = 6, names = FALSE)
      data.frame(metric = metric, node1 = "NMS_mood", node2 = other,
        difference = point["NMS_mood", metric] - point[other, metric], lower = ci[1], upper = ci[2],
        n_bootstrap = length(delta), multiplicity = "unadjusted exploratory comparisons")
    }))
  }
  write.csv(do.call(rbind, intervals), file.path(PD_OUTPUT_DIR, "tables", paste0(prefix, "_intervals.csv")), row.names = FALSE)
  write.csv(do.call(rbind, differences), file.path(PD_OUTPUT_DIR, "tables", paste0(prefix, "_mood_differences.csv")), row.names = FALSE)
  point
}
data <- pd_complete(pd_patients, pd_network_nodes, "network_16_intervals")
if (max(abs(estimate(data[pd_network_nodes]) - archived$adj_matrix)) > 1e-6) {
  stop("Network estimator does not reproduce Step 1.")
}
original_draws <- lapply(archived$boot_result$boots, function(x) metrics(x$graph))
original <- summarize_graphs(archived$adj_matrix, original_draws, "network_16")
nodes <- setdiff(pd_network_nodes, c("PDQ39_emotional", "PDQ39_cognition"))
reduced_data <- pd_complete(pd_patients, nodes, "network_14")[nodes]
graph <- estimate(reduced_data)
B <- pd_count("PD_REDUCED_BOOTSTRAP", 1000)
set.seed(20260912)
indices <- replicate(B, sample.int(nrow(reduced_data), replace = TRUE), simplify = FALSE)
draws <- vector("list", B)
status <- data.frame(replicate = seq_len(B), success = FALSE)
for (i in seq_len(B)) {
  draws[[i]] <- tryCatch(metrics(estimate(reduced_data[indices[[i]], , drop = FALSE])),
                         error = function(e) NULL)
  status$success[i] <- !is.null(draws[[i]])
  if (i %% 100 == 0) message("Reduced-network bootstrap: ", i, "/", B)
}
write.csv(status, file.path(PD_OUTPUT_DIR, "tables", "network_14_bootstrap_status.csv"), row.names = FALSE)
if (!all(status$success)) stop("Reduced-network bootstrap failures require review.")
reduced <- summarize_graphs(graph, draws, "network_14")
comparison <- data.frame(node = nodes, strength_original = original[nodes, "strength"],
  strength_reduced = reduced[nodes, "strength"],
  strength_rank_original = rank(-original[, "strength"], ties.method = "min")[nodes],
  strength_rank_reduced = rank(-reduced[, "strength"], ties.method = "min")[nodes],
  bridge_original = original[nodes, "bridge_strength"], bridge_reduced = reduced[nodes, "bridge_strength"],
  bridge_rank_original = rank(-original[, "bridge_strength"], ties.method = "min")[nodes],
  bridge_rank_reduced = rank(-reduced[, "bridge_strength"], ties.method = "min")[nodes])
write.csv(comparison, file.path(PD_OUTPUT_DIR, "tables", "network_redundancy_sensitivity.csv"), row.names = FALSE)
saveRDS(list(graph = graph, metrics = draws), file.path(PD_OUTPUT_DIR, "data_processed", "network_14_bootstrap.rds"))
pd_finish("step6")
