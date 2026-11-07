
# ============================================================
# MACanalyzeR vs MacrophageFramework benchmark
# Dataset: GSE155960 human adipose tissue macrophages
# ============================================================

# 1. Setup ----------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)
library(MACanalyzeR)
library(devtools)

set.seed(1982)

# Load the preprocessed Seurat object generated from GSE155960.
# This object is required to run the benchmark and is not included
# in this repository due to its file size.

data_file <- "PATH_TO_SEURAT_OBJECT.RDS"
framework_repo <- "PATH_TO_LOCAL_MACROPHAGEFRAMEWORK_REPOSITORY"
output_dir <- "GSE155960_Benchmark"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

load_all(framework_repo)


# 2. Load Seurat object ---------------------------------------

obj <- readRDS(data_file)
DefaultAssay(obj) <- "RNA"

obj$seurat_clusters <- as.character(obj$seurat_clusters)
cluster_levels <- as.character(0:8)
obj$seurat_clusters <- factor(obj$seurat_clusters, levels = cluster_levels)

umap_coordinates <- Embeddings(obj, reduction = "umap.harmony")
colnames(umap_coordinates) <- c("UMAP_1", "UMAP_2")

stopifnot(identical(rownames(umap_coordinates), colnames(obj)))
stopifnot(ncol(obj) == 3367)


# 3. MACanalyzeR ----------------------------------------------

mac_obj <- CreateMacObj(seurat_obj = obj, id = "seurat_clusters", reduction = "umap.harmony", org = "hs")
mac_obj <- MacPolarizeR(mac_obj)

mac_mtx <- mac_obj@CountData@NormalizedCount
pol_genes <- MACanalyzeR::human_PolGenes

m1_genes <- intersect(pol_genes$M1, rownames(mac_mtx))
m2_genes <- intersect(pol_genes$M2, rownames(mac_mtx))
polarization_genes <- c(m1_genes, m2_genes)

polarization_matrix <- t(as.data.frame(mac_mtx[polarization_genes, , drop = FALSE]))

set.seed(1982)
mac_km <- kmeans(polarization_matrix, centers = 3, iter.max = 5000)

mac_cell_table <- data.frame(
  Cell = rownames(polarization_matrix),
  Kmeans = paste0("k", mac_km$cluster),
  stringsAsFactors = FALSE
)

mac_cell_table$State <- recode(
  mac_cell_table$Kmeans,
  k1 = "Transitional",
  k2 = "Inflammatory",
  k3 = "Healing"
)

mac_cell_table$Cluster <- as.character(obj$seurat_clusters[
  match(mac_cell_table$Cell, colnames(obj))
])

mac_km_table <- mac_cell_table |>
  count(Cluster, State, name = "Cells") |>
  group_by(Cluster) |>
  mutate(Percentage = 100 * Cells / sum(Cells)) |>
  ungroup()

mac_table <- mac_km_table |>
  select(Cluster, State, Percentage) |>
  pivot_wider(
    names_from = State,
    values_from = Percentage,
    values_fill = 0,
    names_glue = "{State}_MACanalyzeR"
  ) |>
  arrange(factor(Cluster, levels = cluster_levels))


# 4. MacrophageFramework --------------------------------------

ontology <- load_ontology("v0.4")

framework_result <- analyze_macrophage_profile(
  object = obj,
  ontology = ontology,
  assay = "RNA",
  layer = "counts",
  reduction = "umap.harmony",
  detect_macrophages = FALSE,
  macrophage_clusters = cluster_levels,
  species = "Human",
  output_dir = output_dir,
  save_outputs = FALSE
)

framework_scores <- framework_result$profile$scores
framework_coverage <- framework_result$profile$coverage
framework_coherence <- framework_result$profile$coherence
framework_coherence_IQR <- framework_result$profile$coherence_IQR
framework_annotations <- framework_result$annotations

framework_scores <- framework_scores[, cluster_levels, drop = FALSE]
framework_coverage <- framework_coverage[, cluster_levels, drop = FALSE]
framework_coherence <- framework_coherence[, cluster_levels, drop = FALSE]
framework_coherence_IQR <- framework_coherence_IQR[, cluster_levels, drop = FALSE]

stopifnot(identical(colnames(framework_result$scores), colnames(obj)))
stopifnot(identical(rownames(framework_scores), rownames(framework_coverage)))


# 5. Integrated benchmark table --------------------------------

scores_table <- as.data.frame(t(framework_scores)) |>
  rownames_to_column("Cluster")

coverage_table <- as.data.frame(t(framework_coverage)) |>
  rownames_to_column("Cluster") |>
  rename_with(~ paste0(.x, "_Coverage"), -Cluster)

coherence_table <- as.data.frame(t(framework_coherence)) |>
  rownames_to_column("Cluster") |>
  rename_with(~ paste0(.x, "_Coherence"), -Cluster)

coherence_IQR_table <- as.data.frame(t(framework_coherence_IQR)) |>
  rownames_to_column("Cluster") |>
  rename_with(~ paste0(.x, "_Coherence_IQR"), -Cluster)

benchmark_table <- scores_table |>
  left_join(coverage_table, by = "Cluster") |>
  left_join(coherence_table, by = "Cluster") |>
  left_join(coherence_IQR_table, by = "Cluster") |>
  left_join(mac_table, by = "Cluster") |>
  left_join(framework_annotations, by = "Cluster") |>
  arrange(factor(Cluster, levels = cluster_levels))

write.csv(
  benchmark_table,
  file.path(output_dir, "GSE155960_Benchmark_table.csv"),
  row.names = FALSE
)

write.csv(
  mac_table,
  file.path(output_dir, "MACanalyzeR_cluster_percentages.csv"),
  row.names = FALSE
)

write.csv(
  framework_annotations,
  file.path(output_dir, "MacrophageFramework_annotations.csv"),
  row.names = FALSE
)


# 6. Framework heatmap -----------------------------------------

functional_programs <- c(
  "Inflammatory",
  "Repair-associated",
  "Lipid-associated",
  "ECM-remodeling",
  "IFN-responsive"
)

framework_heatmap <- benchmark_table |>
  select(Cluster, all_of(rownames(framework_scores))) |>
  pivot_longer(
    cols = -Cluster,
    names_to = "Program",
    values_to = "Score"
  ) |>
  group_by(Program) |>
  mutate(Z_score = as.numeric(scale(Score))) |>
  ungroup() |>
  mutate(
    Cluster = factor(Cluster, levels = cluster_levels),
    Program = factor(Program, levels = rev(rownames(framework_scores)))
  )

p_framework <- ggplot(framework_heatmap, aes(Cluster, Program, fill = Z_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B33B3B") +
  labs(x = "Macrophage cluster", y = NULL, fill = "Z-score",
       title = "MacrophageFramework program scores") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold"))


# 7. Framework functional-state composition -------------------

framework_composition <- benchmark_table |>
  select(Cluster, all_of(functional_programs)) |>
  pivot_longer(
    cols = -Cluster,
    names_to = "Program",
    values_to = "Score"
  ) |>
  group_by(Cluster) |>
  mutate(Percentage = 100 * Score / sum(Score, na.rm = TRUE)) |>
  ungroup() |>
  mutate(
    Cluster = factor(Cluster, levels = cluster_levels),
    Program = factor(Program, levels = functional_programs)
  )

p_framework_comp <- ggplot(
  framework_composition,
  aes(Cluster, Percentage, fill = Program)
) +
  geom_col(width = 0.85, color = "white", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Macrophage cluster", y = "Relative score contribution (%)",
       fill = "Functional program",
       title = "MacrophageFramework functional profile") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold"))


# 8. MACanalyzeR heatmap ---------------------------------------

mac_heatmap_data <- mac_table |>
  pivot_longer(
    cols = -Cluster,
    names_to = "State",
    values_to = "Percentage"
  ) |>
  mutate(
    State = sub("_MACanalyzeR", "", State),
    Cluster = factor(Cluster, levels = cluster_levels),
    State = factor(State, levels = c("Inflammatory", "Transitional", "Healing"))
  )

p_mac_heatmap <- ggplot(
  mac_heatmap_data,
  aes(Cluster, State, fill = Percentage)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#B33B3B",
                      limits = c(0, 100)) +
  labs(x = "Macrophage cluster", y = NULL, fill = "Cells (%)",
       title = "MACanalyzeR polarization states") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold"))


# 9. MACanalyzeR composition ----------------------------------

p_mac <- ggplot(
  mac_heatmap_data,
  aes(Cluster, Percentage, fill = State)
) +
  geom_col(width = 0.85, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    Inflammatory = "#C7473A",
    Transitional = "#E5B85C",
    Healing = "#4A9278"
  )) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  labs(x = "Macrophage cluster", y = "Cells (%)",
       fill = "MACanalyzeR state",
       title = "MACanalyzeR polarization states") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold"))


# 10. UMAP data ------------------------------------------------

umap_df <- as.data.frame(umap_coordinates) |>
  rownames_to_column("Cell") |>
  mutate(
    Cluster = as.character(obj$seurat_clusters[match(Cell, colnames(obj))]),
    MACanalyzeR = mac_cell_table$State[match(Cell, mac_cell_table$Cell)],
    Framework = framework_annotations$Full_annotation[
      match(Cluster, framework_annotations$Cluster)
    ]
  )

umap_df$Cluster <- factor(umap_df$Cluster, levels = cluster_levels)

p_umap_clusters <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = Cluster)) +
  geom_point(size = 0.25, alpha = 0.8) +
  labs(title = "Original macrophage clusters", color = "Cluster") +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"))

p_umap_mac <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = MACanalyzeR)) +
  geom_point(size = 0.25, alpha = 0.8) +
  scale_color_manual(values = c(
    Inflammatory = "#C7473A",
    Transitional = "#E5B85C",
    Healing = "#4A9278"
  )) +
  labs(title = "MACanalyzeR states", color = "MACanalyzeR state") +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"))

p_umap_framework <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = Framework)) +
  geom_point(size = 0.25, alpha = 0.8) +
  scale_color_brewer(palette = "Paired") +
  labs(title = "MacrophageFramework annotations", color = "Framework annotation") +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"))


# 11. Save figures ---------------------------------------------

figure_dir <- file.path(output_dir, "figures")

ggsave(file.path(figure_dir, "Figure_1A_Framework_heatmap.pdf"),
       p_framework, width = 7, height = 5.5, device = cairo_pdf)

ggsave(file.path(figure_dir, "Figure_1B_Framework_composition.pdf"),
       p_framework_comp, width = 7, height = 5.5, device = cairo_pdf)

ggsave(file.path(figure_dir, "Figure_1C_MACanalyzeR_heatmap.pdf"),
       p_mac_heatmap, width = 7, height = 4.5, device = cairo_pdf)

ggsave(file.path(figure_dir, "Figure_1D_MACanalyzeR_composition.pdf"),
       p_mac, width = 7, height = 5.5, device = cairo_pdf)

ggsave(file.path(figure_dir, "UMAP_01_Original_clusters.pdf"),
       p_umap_clusters, width = 6, height = 5.5, device = cairo_pdf)

ggsave(file.path(figure_dir, "UMAP_02_MACanalyzeR.pdf"),
       p_umap_mac, width = 6, height = 5.5, device = cairo_pdf)

ggsave(file.path(figure_dir, "UMAP_03_MacrophageFramework.pdf"),
       p_umap_framework, width = 7, height = 5.5, device = cairo_pdf)


# 12. Save R objects -------------------------------------------

saveRDS(mac_obj, file.path(output_dir, "MACanalyzeR_object.rds"))
saveRDS(mac_km, file.path(output_dir, "MACanalyzeR_kmeans.rds"))
saveRDS(framework_result, file.path(output_dir, "MacrophageFramework_result.rds"))
saveRDS(benchmark_table, file.path(output_dir, "Benchmark_table.rds"))

message("Benchmark completed. Results saved in: ", output_dir)
