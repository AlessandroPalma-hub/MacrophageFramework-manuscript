
# ============================================================
# MACROPHAGEFRAMEWORK
# Ontology robustness analysis
# Leave-one-gene-out and random gene-set subsampling
# ============================================================

library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(2026)

library(Seurat)
library(devtools)

devtools::install_github(
  "AlessandroPalma-hub/MacrophageFramework",
  ref = "manuscript-unweighted-v1.0.0",
  upgrade = "never"
)

library(MacrophageFramework)

# ------------------------------------------------------------
# 1. INPUTS AND SETTINGS
# ------------------------------------------------------------

# Path to the preprocessed ovarian Seurat object.
seurat_file <- "PATH_TO_PREPROCESSED_SEURAT_OBJECT.RDS"

output_dir <- "Ontology_Robustness"
plot_dir <- file.path(output_dir, "Plots")
table_dir <- file.path(output_dir, "Tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ontology_version <- "v0.4"
n_iterations <- 1000
sampling_fraction <- 0.80

macrophage_clusters <- c("2", "4", "11", "14", "15", "16")

# ------------------------------------------------------------
# 2. LOAD DATA AND ONTOLOGY
# ------------------------------------------------------------

ovarian <- readRDS(seurat_file)
ontology_v04 <- load_ontology(ontology_version)

macrophage_cells <- colnames(ovarian)[
  ovarian$seurat_clusters %in% macrophage_clusters
]

counts_macrophages <- GetAssayData(
  ovarian,
  assay = "RNA",
  layer = "counts"
)[, macrophage_cells, drop = FALSE]

cluster_id <- setNames(
  as.character(
    ovarian$seurat_clusters[
      match(macrophage_cells, colnames(ovarian))
    ]
  ),
  macrophage_cells
)

available_genes <- intersect(
  ontology_v04$Gene,
  rownames(counts_macrophages)
)

ontology_available <- ontology_v04[
  ontology_v04$Gene %in% available_genes,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 3. BASELINE SCORES
# ------------------------------------------------------------

baseline <- MacrophageFramework:::score_programs(
  expression = counts_macrophages,
  ontology = ontology_v04
)

programs <- rownames(baseline$scores)

aggregate_scores <- function(score_matrix) {
  result <- sapply(
    macrophage_clusters,
    function(cl) {
      Matrix::rowMeans(
        score_matrix[
          ,
          cluster_id == cl,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    }
  )

  result[
    programs,
    macrophage_clusters,
    drop = FALSE
  ]
}

baseline_cluster_scores <- aggregate_scores(baseline$scores)

baseline_dominant <- apply(
  baseline_cluster_scores,
  2,
  function(x) rownames(baseline_cluster_scores)[which.max(x)]
)

# ------------------------------------------------------------
# 4. LEAVE-ONE-GENE-OUT ANALYSIS
# ------------------------------------------------------------

loo_results <- vector(
  "list",
  length(available_genes)
)

names(loo_results) <- available_genes

for (gene in available_genes) {

  affected_program <- ontology_available$Program[
    match(gene, ontology_available$Gene)
  ]

  ontology_loo <- ontology_available[
    ontology_available$Gene != gene,
    ,
    drop = FALSE
  ]

  result_loo <- MacrophageFramework:::score_programs(
    expression = counts_macrophages,
    ontology = ontology_loo
  )

  scores_loo <- aggregate_scores(result_loo$scores)

  program_cor <- cor(
    baseline_cluster_scores[affected_program, ],
    scores_loo[affected_program, ],
    method = "spearman"
  )

  profile_cor <- sapply(
    macrophage_clusters,
    function(cl) {
      cor(
        baseline_cluster_scores[, cl],
        scores_loo[, cl],
        method = "spearman"
      )
    }
  )

  delta <- scores_loo - baseline_cluster_scores

  dominant_loo <- apply(
    scores_loo,
    2,
    function(x) rownames(scores_loo)[which.max(x)]
  )

  loo_results[[gene]] <- data.frame(
    Gene = gene,
    Program = affected_program,
    Program_Spearman = program_cor,
    Mean_Profile_Spearman = mean(profile_cor, na.rm = TRUE),
    Min_Profile_Spearman = min(profile_cor, na.rm = TRUE),
    MAE = mean(abs(delta)),
    RMSE = sqrt(mean(delta^2)),
    Dominant_program_agreement = mean(
      dominant_loo == baseline_dominant
    )
  )
}

loo_summary <- bind_rows(loo_results)

write.csv(
  loo_summary,
  file.path(table_dir, "Leave_one_gene_out_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. RANDOM GENE-SET SUBSAMPLING
# ------------------------------------------------------------

program_stability <- matrix(
  NA_real_,
  nrow = n_iterations,
  ncol = length(programs),
  dimnames = list(NULL, programs)
)

cluster_stability <- matrix(
  NA_real_,
  nrow = n_iterations,
  ncol = length(macrophage_clusters),
  dimnames = list(NULL, macrophage_clusters)
)

dominant_by_iteration <- matrix(
  NA_character_,
  nrow = n_iterations,
  ncol = length(macrophage_clusters),
  dimnames = list(NULL, macrophage_clusters)
)

score_perturbations <- array(
  NA_real_,
  dim = c(
    n_iterations,
    length(programs),
    length(macrophage_clusters)
  ),
  dimnames = list(
    NULL,
    programs,
    macrophage_clusters
  )
)

subsampling_results <- vector("list", n_iterations)

for (iteration in seq_len(n_iterations)) {

  # Randomly subsample genes independently within each program
  ontology_subsampled <- do.call(
    rbind,
    lapply(
      split(ontology_available, ontology_available$Program),
      function(program_df) {

        n_selected <- max(
          2,
          floor(sampling_fraction * nrow(program_df))
        )

        program_df[
          sample(seq_len(nrow(program_df)), n_selected),
          ,
          drop = FALSE
        ]
      }
    )
  )

  result_sub <- MacrophageFramework:::score_programs(
    expression = counts_macrophages,
    ontology = ontology_subsampled
  )

  scores_sub <- aggregate_scores(result_sub$scores)

  score_perturbations[iteration, , ] <- scores_sub

  # Program-level stability across clusters
  program_stability[iteration, ] <- sapply(
    programs,
    function(program) {
      cor(
        baseline_cluster_scores[program, ],
        scores_sub[program, ],
        method = "spearman"
      )
    }
  )

  # Cluster-level profile stability across programs
  cluster_stability[iteration, ] <- sapply(
    macrophage_clusters,
    function(cl) {
      cor(
        baseline_cluster_scores[, cl],
        scores_sub[, cl],
        method = "spearman"
      )
    }
  )

  # Dominant program identity
  dominant_sub <- apply(
    scores_sub,
    2,
    function(x) rownames(scores_sub)[which.max(x)]
  )

  dominant_by_iteration[iteration, ] <- dominant_sub

  # Quantitative differences
  delta <- scores_sub - baseline_cluster_scores

  subsampling_results[[iteration]] <- data.frame(
    Iteration = iteration,
    Mean_program_Spearman = mean(
      program_stability[iteration, ],
      na.rm = TRUE
    ),
    Min_program_Spearman = min(
      program_stability[iteration, ],
      na.rm = TRUE
    ),
    Mean_profile_Spearman = mean(
      cluster_stability[iteration, ],
      na.rm = TRUE
    ),
    Min_profile_Spearman = min(
      cluster_stability[iteration, ],
      na.rm = TRUE
    ),
    MAE = mean(abs(delta)),
    RMSE = sqrt(mean(delta^2)),
    Dominant_program_agreement = mean(
      dominant_sub == baseline_dominant
    )
  )

  if (iteration %% 100 == 0) {
    message("Completed ", iteration, "/", n_iterations)
  }
}

subsampling_summary <- bind_rows(subsampling_results)

write.csv(
  subsampling_summary,
  file.path(table_dir, "Subsampling_iteration_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. SUBSAMPLING SUMMARY STATISTICS
# ------------------------------------------------------------

metrics <- setdiff(
  colnames(subsampling_summary),
  "Iteration"
)

subsampling_statistics <- data.frame(
  Metric = metrics,
  Mean = sapply(
    subsampling_summary[metrics],
    mean,
    na.rm = TRUE
  ),
  SD = sapply(
    subsampling_summary[metrics],
    sd,
    na.rm = TRUE
  ),
  Median = sapply(
    subsampling_summary[metrics],
    median,
    na.rm = TRUE
  ),
  Q025 = sapply(
    subsampling_summary[metrics],
    quantile,
    probs = 0.025,
    na.rm = TRUE
  ),
  Q975 = sapply(
    subsampling_summary[metrics],
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
)

write.csv(
  subsampling_statistics,
  file.path(table_dir, "Subsampling_summary_statistics.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. DOMINANT PROGRAM FREQUENCY
# ------------------------------------------------------------

dominant_frequency <- do.call(
  rbind,
  lapply(
    macrophage_clusters,
    function(cl) {
      prop.table(
        table(
          factor(
            dominant_by_iteration[, cl],
            levels = programs
          )
        )
      )
    }
  )
)

rownames(dominant_frequency) <- macrophage_clusters

write.csv(
  dominant_frequency,
  file.path(table_dir, "Dominant_program_frequency.csv")
)

# ------------------------------------------------------------
# 8. PROGRAM- AND CLUSTER-LEVEL STABILITY SUMMARIES
# ------------------------------------------------------------

program_stability_summary <- data.frame(
  Program = programs,
  Mean = apply(program_stability, 2, mean, na.rm = TRUE),
  SD = apply(program_stability, 2, sd, na.rm = TRUE),
  Median = apply(program_stability, 2, median, na.rm = TRUE),
  Q025 = apply(
    program_stability, 2, quantile,
    probs = 0.025, na.rm = TRUE
  ),
  Q975 = apply(
    program_stability, 2, quantile,
    probs = 0.975, na.rm = TRUE
  )
)

cluster_stability_summary <- data.frame(
  Cluster = macrophage_clusters,
  Mean = apply(cluster_stability, 2, mean, na.rm = TRUE),
  SD = apply(cluster_stability, 2, sd, na.rm = TRUE),
  Median = apply(cluster_stability, 2, median, na.rm = TRUE),
  Q025 = apply(
    cluster_stability, 2, quantile,
    probs = 0.025, na.rm = TRUE
  ),
  Q975 = apply(
    cluster_stability, 2, quantile,
    probs = 0.975, na.rm = TRUE
  )
)

write.csv(
  program_stability_summary,
  file.path(table_dir, "Program_level_stability.csv"),
  row.names = FALSE
)

write.csv(
  cluster_stability_summary,
  file.path(table_dir, "Cluster_level_stability.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. PLOT: CORRELATION DISTRIBUTIONS
# ------------------------------------------------------------

correlation_data <- subsampling_summary %>%
  select(
    Iteration,
    Mean_program_Spearman,
    Mean_profile_Spearman
  ) %>%
  pivot_longer(
    cols = -Iteration,
    names_to = "Metric",
    values_to = "Correlation"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      Mean_program_Spearman = "Program-level",
      Mean_profile_Spearman = "Multidimensional profile"
    )
  )

p1 <- ggplot(
  correlation_data,
  aes(x = Correlation, fill = Metric)
) +
  geom_density(alpha = 0.45) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey40"
  ) +
  facet_wrap(~Metric, ncol = 1) +
  labs(
    x = "Spearman correlation",
    y = "Density",
    title = "Stability of macrophage program scores"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 10. PLOT: QUANTITATIVE SCORE DIFFERENCES
# ------------------------------------------------------------

difference_data <- subsampling_summary %>%
  select(Iteration, MAE, RMSE) %>%
  pivot_longer(
    cols = -Iteration,
    names_to = "Metric",
    values_to = "Difference"
  )

p2 <- ggplot(
  difference_data,
  aes(x = Difference)
) +
  geom_density(
    fill = "steelblue",
    alpha = 0.5
  ) +
  facet_wrap(~Metric, scales = "free", ncol = 1) +
  labs(
    x = "Score difference",
    y = "Density",
    title = "Quantitative deviation from baseline scores"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 11. PLOT: DOMINANT PROGRAM AGREEMENT
# ------------------------------------------------------------

agreement_data <- subsampling_summary %>%
  count(Dominant_program_agreement) %>%
  mutate(
    Percentage = 100 * n / sum(n),
    Agreement = factor(
      Dominant_program_agreement,
      levels = c(2/3, 5/6, 1),
      labels = c("4/6 clusters", "5/6 clusters", "6/6 clusters")
    )
  )

p3 <- ggplot(
  agreement_data,
  aes(x = Agreement, y = Percentage)
) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%% (n=%d)", Percentage, n)),
    vjust = -0.4,
    size = 3.5
  ) +
  labs(
    x = "Clusters retaining the baseline dominant program",
    y = "Iterations (%)",
    title = "Dominant program agreement"
  ) +
  scale_y_continuous(
    limits = c(0, max(agreement_data$Percentage) * 1.2)
  ) +
  theme_classic(base_size = 12)

# ------------------------------------------------------------
# 12. PLOT: DOMINANT PROGRAM FREQUENCY HEATMAP
# ------------------------------------------------------------

heatmap_data <- as.data.frame(
  as.table(dominant_frequency)
)

colnames(heatmap_data) <- c(
  "Cluster", "Program", "Frequency"
)

heatmap_data <- heatmap_data %>%
  mutate(
    Cluster = factor(Cluster, levels = macrophage_clusters),
    Program = factor(Program, levels = programs),
    Frequency_percent = Frequency * 100
  )

p4 <- ggplot(
  heatmap_data,
  aes(x = Program, y = Cluster, fill = Frequency_percent)
) +
  geom_tile(
    color = "grey75",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.0f", Frequency_percent)),
    size = 3.5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "steelblue",
    limits = c(0, 100),
    name = "Dominance\nfrequency (%)"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL,
    y = "Macrophage cluster",
    title = "Robustness of dominant macrophage programs",
    subtitle = paste(
      "Frequency of each program being dominant across",
      n_iterations,
      "gene-set perturbations"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 13. PLOT: PROGRAM-LEVEL STABILITY
# ------------------------------------------------------------

program_stability_df <- as.data.frame(
  program_stability
) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Program",
    values_to = "Spearman"
  )

program_order <- names(
  sort(
    apply(program_stability, 2, median, na.rm = TRUE),
    decreasing = TRUE
  )
)

program_stability_df$Program <- factor(
  program_stability_df$Program,
  levels = program_order
)

p5 <- ggplot(
  program_stability_df,
  aes(x = Program, y = Spearman)
) +
  geom_boxplot(
    fill = "steelblue",
    alpha = 0.65,
    width = 0.65,
    outlier.shape = NA
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2,
    fill = "black"
  ) +
  coord_flip() +
  scale_y_continuous(
    breaks = seq(0.4, 1, 0.1)
  ) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(
    title = "Program-level stability",
    x = NULL,
    y = "Spearman correlation"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 14. PLOT: CLUSTER-LEVEL PROFILE STABILITY
# ------------------------------------------------------------

cluster_stability_df <- as.data.frame(
  cluster_stability
) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Cluster",
    values_to = "Spearman"
  )

cluster_order <- names(
  sort(
    apply(cluster_stability, 2, median, na.rm = TRUE),
    decreasing = TRUE
  )
)

cluster_stability_df$Cluster <- factor(
  cluster_stability_df$Cluster,
  levels = cluster_order
)

p6 <- ggplot(
  cluster_stability_df,
  aes(x = Cluster, y = Spearman)
) +
  geom_boxplot(
    fill = "darkseagreen3",
    alpha = 0.7,
    width = 0.65,
    outlier.shape = NA
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2,
    fill = "black"
  ) +
  coord_flip() +
  scale_y_continuous(
    breaks = seq(0.4, 1, 0.1)
  ) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(
    title = "Cluster-level profile stability",
    x = "Macrophage cluster",
    y = "Spearman correlation"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 15. SAVE ALL FIGURES
# ------------------------------------------------------------

plots <- list(
  Subsampling_correlations = p1,
  Subsampling_score_differences = p2,
  Subsampling_dominant_agreement = p3,
  Subsampling_dominant_program_heatmap = p4,
  Program_level_stability = p5,
  Cluster_level_stability = p6
)

for (plot_name in names(plots)) {
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(plot_name, ".pdf")
    ),
    plot = plots[[plot_name]],
    width = if (plot_name == "Subsampling_dominant_program_heatmap") 9 else 7,
    height = if (plot_name == "Subsampling_dominant_program_heatmap") 4.5 else 5
  )
}

# ------------------------------------------------------------
# 16. SAVE REPRODUCIBILITY INFORMATION
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

saveRDS(
  list(
    baseline_cluster_scores = baseline_cluster_scores,
    loo_results = loo_summary,
    subsampling_results = subsampling_summary,
    program_stability = program_stability,
    cluster_stability = cluster_stability,
    dominant_by_iteration = dominant_by_iteration,
    score_perturbations = score_perturbations
  ),
  file.path(output_dir, "robustness_analysis_results.rds")
)

cat("\n========================================\n")
cat("ONTOLOGY ROBUSTNESS ANALYSIS COMPLETED\n")
cat("========================================\n")
cat("Available ontology genes:", length(available_genes), "/", nrow(ontology_v04), "\n")
cat("Subsampling iterations:", n_iterations, "\n")
cat("Results directory:", output_dir, "\n")
