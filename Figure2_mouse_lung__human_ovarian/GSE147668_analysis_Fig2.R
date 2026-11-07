
# ============================================================
# 01_preprocess_GSE147668.R
#
# Dataset: GSE147668
# Species: Mus musculus
# Tissue: Lung
#
# Purpose:
#   Import raw count data, extract sample metadata, perform
#   Seurat preprocessing and clustering, and run the
#   MacrophageFramework analysis.
#
# Framework version: ontology v0.4
# Random seed: 1982
# ============================================================


# 1. Load packages --------------------------------------------

library(data.table)
library(Matrix)
library(Seurat)

# Load the current MacrophageFramework source code
devtools::load_all(
  path = "/Users/alessandro/Desktop/Dati/PROGETTI ONGOING/Progetto Computational Framework Macrofagi/MacrophageFramework/MacrophageFramework",
  quiet = TRUE
)


# 2. Configuration --------------------------------------------

PROJECT_DIR <- paste0(
  "/Users/alessandro/Desktop/Dati/PROGETTI ONGOING/",
  "Progetto Computational Framework Macrofagi/",
  "MacrophageFramework/MacrophageFramework"
)

COUNTS_FILE <- file.path(
  PROJECT_DIR,
  "..",
  "datasets",
  "GSE147668_RAW",
  "GSE147668_table_counts_lungimmune.tsv"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "GSE147668_final"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

set.seed(1982)


# 3. Import raw count matrix ----------------------------------

lung_df <- data.table::fread(
  COUNTS_FILE,
  data.table = FALSE,
  check.names = FALSE
)

gene_names <- lung_df[[1]]

counts <- as.matrix(
  lung_df[, -1, drop = FALSE]
)

rownames(counts) <- make.unique(
  as.character(gene_names)
)

counts <- Matrix::Matrix(
  counts,
  sparse = TRUE
)

rm(lung_df)


# 4. Validate input data --------------------------------------

if (anyNA(counts@x)) {
  stop("The count matrix contains missing values.")
}

if (any(counts@x < 0)) {
  stop("The count matrix contains negative values.")
}

if (anyDuplicated(colnames(counts))) {
  stop("Duplicated cell identifiers detected.")
}


# 5. Create Seurat object -------------------------------------

lung <- CreateSeuratObject(
  counts = counts,
  project = "GSE147668"
)

rm(counts)


# 6. Extract sample metadata ---------------------------------

cell_ids <- colnames(lung)

parts <- strsplit(
  cell_ids,
  "_",
  fixed = TRUE
)

if (any(lengths(parts) != 7)) {
  stop("Unexpected cell identifier structure.")
}

lung$sex <- vapply(parts, `[`, character(1), 1)
lung$stage <- vapply(parts, `[`, character(1), 2)
lung$condition <- vapply(parts, `[`, character(1), 3)
lung$sorting_population <- vapply(parts, `[`, character(1), 4)
lung$replicate <- vapply(parts, `[`, character(1), 5)
lung$plate <- vapply(parts, `[`, character(1), 6)
lung$well <- vapply(parts, `[`, character(1), 7)

rm(parts, cell_ids)


# 7. Seurat preprocessing -------------------------------------

DefaultAssay(lung) <- "RNA"

# Library-size normalization
lung <- NormalizeData(
  lung,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# Highly variable genes
lung <- FindVariableFeatures(
  lung,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# Scaling
lung <- ScaleData(
  lung,
  features = VariableFeatures(lung),
  verbose = FALSE
)

# Principal component analysis
lung <- RunPCA(
  lung,
  features = VariableFeatures(lung),
  npcs = 30,
  verbose = FALSE
)

# Nearest-neighbour graph
lung <- FindNeighbors(
  lung,
  dims = 1:20,
  verbose = FALSE
)

# Graph-based clustering
lung <- FindClusters(
  lung,
  resolution = 0.3,
  random.seed = 1982,
  verbose = FALSE
)

# UMAP dimensionality reduction
lung <- RunUMAP(
  lung,
  dims = 1:20,
  seed.use = 1982,
  verbose = FALSE
)


# 8. Save preprocessing outputs -------------------------------

saveRDS(
  lung,
  file = file.path(
    OUTPUT_DIR,
    "GSE147668_lung_mouse_Seurat.rds"
  )
)

write.csv(
  lung@meta.data,
  file = file.path(
    OUTPUT_DIR,
    "GSE147668_cell_metadata.csv"
  ),
  row.names = TRUE
)

write.csv(
  as.data.frame(table(Idents(lung))),
  file = file.path(
    OUTPUT_DIR,
    "GSE147668_cluster_sizes.csv"
  ),
  row.names = FALSE
)


# 9. MacrophageFramework analysis -----------------------------

lung_result <- analyze_macrophage_profile(
  object = lung,
  assay = "RNA",
  reduction = "umap",
  species="Mouse",
  min_core_score = 0.1,
  output_dir = OUTPUT_DIR,
  save_outputs = TRUE,
  detect_macrophages = TRUE
)


# 10. Save computational environment -------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    OUTPUT_DIR,
    "GSE147668_sessionInfo.txt"
  )
)

# --------------------------------------------------
# Multidimensional macrophage profile plot
# --------------------------------------------------

score_df <- as.data.frame(lung_result$profile$scores)
score_df$Program <- rownames(score_df)

score_long <- tidyr::pivot_longer(
  score_df,
  cols = -Program,
  names_to = "Cluster",
  values_to = "Score"
)

ann_df <- lung_result$annotations |>
  dplyr::mutate(Cluster = as.character(Cluster))

profile_df <- dplyr::left_join(
  score_long,
  ann_df,
  by = "Cluster"
)

profile_plot_df <- profile_df |>
  dplyr::group_by(Full_annotation, Program) |>
  dplyr::summarise(
    Score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(Full_annotation) |>
  dplyr::mutate(
    Norm = sqrt(sum(Score^2)),
    Normalized_score = Score / Norm
  ) |>
  dplyr::ungroup()

program_order <- c(
  "Monocyte-associated",
  "Resident-associated",
  "Inflammatory",
  "Repair-associated",
  "Lipid-associated",
  "ECM-remodeling",
  "IFN-responsive"
)

profile_plot_df$Program <- factor(
  profile_plot_df$Program,
  levels = program_order
)

profile_plot <- ggplot2::ggplot(
  profile_plot_df,
  ggplot2::aes(
    x = Program,
    y = Normalized_score,
    group = Full_annotation,
    color = Full_annotation
  )
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::labs(
    title = "Multidimensional macrophage profiles",
    x = NULL,
    y = "Normalized program intensity",
    color = "Macrophage profile"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "figures",
    "Profile_plot_GSE147668.pdf"
  ),
  plot = profile_plot,
  width = 14,
  height = 7
)
