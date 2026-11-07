# CellTypist benchmark: GSE155960 adipose macrophage-enriched subset
# CellTypist 1.7.1; Immune_All_Low.pkl

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Load the macrophage-enriched Seurat object and Framework results.
base_dir <- "benchmark_CellTypist/GSE155960"
input_dir <- file.path(base_dir, "CellTypist_input")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# Load the preprocessed Seurat object generated from GSE155960.
# This object is required for the benchmark and is not included
# in this repository due to its file size.
obj_adipose <- readRDS("PATH_TO_SEURAT_OBJECT.RDS")

framework_result <- readRDS(file.path(base_dir, "framework_result.rds"))

# Export raw counts and cell/gene identifiers for CellTypist.
counts <- GetAssayData(obj_adipose, assay = "RNA", layer = "counts")
writeMM(counts, file.path(input_dir, "matrix.mtx"))
write.table(rownames(counts), file.path(input_dir, "genes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(counts), file.path(input_dir, "cells.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Import CellTypist labels and the complete class-score matrix.
celltypist <- read.csv(file.path(base_dir, "CellTypist_predictions.csv"),
                       stringsAsFactors = FALSE)
celltypist_scores <- read.csv(file.path(base_dir, "CellTypist_probabilities.csv"),
                               check.names = FALSE, stringsAsFactors = FALSE)

celltypist <- celltypist[match(colnames(obj_adipose), celltypist$Cell), ]
celltypist_scores <- celltypist_scores[
  match(colnames(obj_adipose), celltypist_scores$Cell), ]

obj_adipose$CellTypist_label <- celltypist$CellTypist_label
obj_adipose$CellTypist_probability <- celltypist$CellTypist_probability

# Group CellTypist classes into broad lineage categories.
macrophage_labels <- c(
  "Macrophages", "Intermediate macrophages", "Mono-mac",
  "Alveolar macrophages", "Kupffer cells", "Intestinal macrophages",
  "Hofbauer cells", "Erythrophagocytic macrophages"
)
monocyte_labels <- c(
  "Classical monocytes", "Non-classical monocytes", "Monocytes"
)

class_scores <- as.matrix(celltypist_scores[, -1])
storage.mode(class_scores) <- "numeric"

monocyte_score <- rowSums(class_scores[, monocyte_labels, drop = FALSE])
macrophage_score <- rowSums(class_scores[, macrophage_labels, drop = FALSE])
other_labels <- setdiff(colnames(class_scores),
                        c(monocyte_labels, macrophage_labels))
other_score <- rowSums(class_scores[, other_labels, drop = FALSE])
total_score <- rowSums(class_scores)

obj_adipose$CT_Monocyte_relative <- monocyte_score / total_score
obj_adipose$CT_Macrophage_relative <- macrophage_score / total_score
obj_adipose$CT_Other_relative <- other_score / total_score

obj_adipose$CellTypist_broad <- ifelse(
  obj_adipose$CellTypist_label %in% macrophage_labels, "Macrophage-like",
  ifelse(obj_adipose$CellTypist_label %in% monocyte_labels,
         "Monocyte-like", "Other")
)

# Transfer cluster-level Framework lineage annotations to individual cells.
lineage_map <- setNames(
  framework_result$annotations$Lineage,
  as.character(framework_result$annotations$Cluster)
)
cell_lineage <- unname(
  lineage_map[as.character(obj_adipose$harmony_clusters)]
)
names(cell_lineage) <- colnames(obj_adipose)
obj_adipose$Framework_lineage <- cell_lineage

# Add cell-level Framework program scores.
framework_scores <- as.data.frame(t(framework_result$scores))
framework_scores <- framework_scores[
  match(colnames(obj_adipose), rownames(framework_scores)), , drop = FALSE
]
for (program in colnames(framework_scores)) {
  score_vector <- framework_scores[[program]]
  names(score_vector) <- rownames(framework_scores)
  obj_adipose[[program]] <- score_vector
}

# Panel A: CellTypist broad annotation on the Harmony UMAP.
pA <- DimPlot(
  obj_adipose, reduction = "umap.harmony",
  group.by = "CellTypist_broad", pt.size = 0.4
) +
  ggtitle("A. CellTypist") +
  theme_classic()

# Panel B: Framework lineage on the same UMAP.
pB <- DimPlot(
  obj_adipose, reduction = "umap.harmony",
  group.by = "Framework_lineage", pt.size = 0.4
) +
  ggtitle("B. MacrophageFramework") +
  theme_classic()

# Panel C: Relative functional-program enrichment across macrophage-related labels.
functional_programs <- c(
  "Inflammatory", "Repair-associated", "Lipid-associated",
  "ECM-remodeling", "IFN-responsive"
)

functional_data <- obj_adipose[[]] %>%
  filter(CellTypist_label %in% macrophage_labels) %>%
  group_by(CellTypist_label) %>%
  summarise(
    across(all_of(functional_programs), median, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(functional_programs),
    names_to = "Program",
    values_to = "Median_score"
  ) %>%
  group_by(Program) %>%
  mutate(Z_score = as.numeric(scale(Median_score))) %>%
  ungroup()

pC <- ggplot(
  functional_data,
  aes(x = Program, y = CellTypist_label, fill = Z_score)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#2878A5", mid = "white", high = "#C43C39",
    midpoint = 0, name = "Relative\nZ-score"
  ) +
  labs(
    title = "C. Functional heterogeneity",
    x = NULL, y = NULL
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel D: CellTypist relative lineage support, coloured by Framework lineage.
pD <- ggplot(
  obj_adipose[[]],
  aes(
    x = CT_Monocyte_relative,
    y = CT_Macrophage_relative,
    color = Framework_lineage
  )
) +
  geom_point(size = 0.7, alpha = 0.25) +
  stat_density_2d(
    aes(group = Framework_lineage),
    linewidth = 0.8, alpha = 0.8, bins = 5
  ) +
  scale_color_manual(
    values = c(
      "Monocyte-associated" = "#C43C39",
      "Resident-associated" = "#2878A5"
    )
  ) +
  labs(
    title = "D. CellTypist scores and Framework lineage",
    x = "CellTypist Monocyte-like support",
    y = "CellTypist Macrophage-like support",
    color = "MacrophageFramework lineage"
  ) +
  theme_classic(base_size = 11)

# Panel E: CellTypist label composition within Framework clusters.
label_counts <- table(obj_adipose$CellTypist_label)
top_labels <- names(sort(label_counts, decreasing = TRUE))[1:10]

plot_data <- obj_adipose[[]] %>%
  mutate(
    CellTypist_plot = ifelse(
      CellTypist_label %in% top_labels, CellTypist_label, "Other"
    ),
    harmony_clusters = factor(harmony_clusters, levels = as.character(0:8))
  ) %>%
  count(harmony_clusters, CellTypist_plot) %>%
  group_by(harmony_clusters) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  ungroup()

pE <- ggplot(
  plot_data,
  aes(x = harmony_clusters, y = Percentage, fill = CellTypist_plot)
) +
  geom_col(width = 0.85, color = "white", linewidth = 0.15) +
  scale_y_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = "E. CellTypist annotation across Framework clusters",
    x = "MacrophageFramework cluster",
    y = "Cells (%)",
    fill = "CellTypist annotation"
  ) +
  theme_classic(base_size = 11)

# Assemble and save the benchmark figure.
figure_benchmark <- (pA + pB) / (pC + pD) / pE +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(base_dir, "Figure_CellTypist_benchmark.pdf"),
  plot = figure_benchmark,
  width = 14, height = 16, device = cairo_pdf
)

# Save benchmark tables and annotated object.
write.csv(
  as.data.frame(table(obj_adipose$Framework_lineage,
                      obj_adipose$CellTypist_label)),
  file.path(base_dir, "CellTypist_by_Framework_lineage.csv"),
  row.names = FALSE
)
write.csv(
  functional_data,
  file.path(base_dir, "CellTypist_functional_program_summary.csv"),
  row.names = FALSE
)
saveRDS(obj_adipose, file.path(base_dir, "GSE155960_CellTypist_benchmark.rds"))
