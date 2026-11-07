
# ==================================================
# GSE224392 — Human ovarian cancer
# MacrophageFramework analysis
# ==================================================

library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------------------------------------
# Project paths
# --------------------------------------------------

# --------------------------------------------------
# Project paths
# --------------------------------------------------

# Set this path to the directory containing the downloaded
# GSE224392 10x files (barcodes, features and matrix).
DATA_DIR <- "PATH_TO_GSE224392_RAW_DIRECTORY"

OUTPUT_DIR <- "GSE224392_final"

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(
  FIGURES_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

# --------------------------------------------------
# Load MacrophageFramework release used in the manuscript
# --------------------------------------------------

devtools::install_github(
  "AlessandroPalma-hub/MacrophageFramework",
  ref = "manuscript-unweighted-v1.0.0",
  upgrade = "never"
)

library(MacrophageFramework)

# --------------------------------------------------
# Reproducibility
# --------------------------------------------------

set.seed(1982)

# --------------------------------------------------
# Load and combine 10x samples
# --------------------------------------------------

load_10x_samples <- function(
    data_dir,
    samples = paste0("sample", 1:4)
) {

  matrices <- list()

  for (sample in samples) {

    barcode_file <- list.files(
      data_dir,
      pattern = paste0("_", sample, "_barcodes\\.tsv\\.gz$"),
      full.names = TRUE
    )

    feature_file <- list.files(
      data_dir,
      pattern = paste0("_", sample, "_features\\.tsv\\.gz$"),
      full.names = TRUE
    )

    matrix_file <- list.files(
      data_dir,
      pattern = paste0("_", sample, "_matrix\\.mtx\\.gz$"),
      full.names = TRUE
    )

    if (
      length(barcode_file) != 1 ||
      length(feature_file) != 1 ||
      length(matrix_file) != 1
    ) {
      stop("Could not identify files for ", sample)
    }

    mat <- Matrix::readMM(
      gzfile(matrix_file)
    )

    features <- read.delim(
      gzfile(feature_file),
      header = FALSE,
      stringsAsFactors = FALSE
    )

    barcodes <- read.delim(
      gzfile(barcode_file),
      header = FALSE,
      stringsAsFactors = FALSE
    )

    gene_names <- make.unique(features[[2]])
    cell_names <- barcodes[[1]]

    rownames(mat) <- gene_names
    colnames(mat) <- paste0(sample, "_", cell_names)

    matrices[[sample]] <- mat
  }

  expression <- do.call(cbind, matrices)

  return(expression)
}

ovarian_counts <- load_10x_samples(
  data_dir = DATA_DIR,
  samples = paste0("sample", 1:4)
)

# --------------------------------------------------
# Create Seurat object
# --------------------------------------------------

ovarian <- CreateSeuratObject(
  counts = ovarian_counts,
  project = "GSE224392"
)

ovarian$Sample <- sub(
  "_.*$",
  "",
  colnames(ovarian)
)

# --------------------------------------------------
# Seurat preprocessing
# --------------------------------------------------

ovarian <- NormalizeData(
  ovarian,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

ovarian <- FindVariableFeatures(
  ovarian,
  selection.method = "vst",
  nfeatures = 2000
)

ovarian <- ScaleData(
  ovarian,
  features = VariableFeatures(ovarian)
)

ovarian <- RunPCA(
  ovarian,
  features = VariableFeatures(ovarian),
  npcs = 30
)

ovarian <- FindNeighbors(
  ovarian,
  dims = 1:20
)

ovarian <- FindClusters(
  ovarian,
  resolution = 0.3,
  random.seed = 1982
)

ovarian <- RunUMAP(
  ovarian,
  dims = 1:20,
  seed.use = 1982
)

# --------------------------------------------------
# Save processed Seurat object and metadata
# --------------------------------------------------

saveRDS(
  ovarian,
  file = file.path(
    OUTPUT_DIR,
    "GSE224392_processed.rds"
  )
)

write.csv(
  ovarian@meta.data,
  file = file.path(
    OUTPUT_DIR,
    "cell_metadata.csv"
  )
)

write.csv(
  as.data.frame(table(Idents(ovarian))),
  file = file.path(
    OUTPUT_DIR,
    "cluster_sizes.csv"
  ),
  row.names = FALSE
)

# --------------------------------------------------
# MacrophageFramework analysis
# --------------------------------------------------

human_result <- analyze_macrophage_profile(
  object = ovarian,
  assay = "RNA",
  reduction = "umap",
  species = "Human",
  min_core_score = 0.1,
  output_dir = file.path(
    OUTPUT_DIR,
    "MacrophageFramework_output"
  ),
  save_outputs = TRUE,
  detect_macrophages = TRUE
)

# --------------------------------------------------
# Feature plots
# --------------------------------------------------

feature_genes <- c(
  "CD68",
  "CD3D",
  "CD79A"
)

for (gene in feature_genes) {

  p <- FeaturePlot(
    ovarian,
    features = gene,
    reduction = "umap"
  ) +
    ggtitle(gene) +
    theme_classic(base_size = 14)

  ggsave(
    filename = file.path(
      FIGURES_DIR,
      paste0("UMAP_", gene, ".pdf")
    ),
    plot = p,
    width = 7,
    height = 6
  )
}

# --------------------------------------------------
# Multidimensional macrophage profile plot
# --------------------------------------------------

score_df <- as.data.frame(
  human_result$profile$scores
)

score_df$Program <- rownames(score_df)

score_long <- pivot_longer(
  score_df,
  cols = -Program,
  names_to = "Cluster",
  values_to = "Score"
)

ann_df <- human_result$annotations %>%
  mutate(
    Cluster = as.character(Cluster)
  )

profile_df <- left_join(
  score_long,
  ann_df,
  by = "Cluster"
)

profile_plot_df <- profile_df %>%
  group_by(Full_annotation, Program) %>%
  summarise(
    Score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Full_annotation) %>%
  mutate(
    Norm = sqrt(sum(Score^2)),
    Normalized_score = Score / Norm
  ) %>%
  ungroup()

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

profile_plot <- ggplot(
  profile_plot_df,
  aes(
    x = Program,
    y = Normalized_score,
    group = Full_annotation,
    color = Full_annotation
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  theme_classic(base_size = 13) +
  labs(
    title = "Multidimensional macrophage profiles",
    x = NULL,
    y = "Normalized program intensity",
    color = "Macrophage profile"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "Profile_plot_GSE224392.pdf"
  ),
  plot = profile_plot,
  width = 14,
  height = 7
)

# --------------------------------------------------
# Session information
# --------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)
