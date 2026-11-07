# ============================================================
# Tabula Muris Senis macrophage atlas and aging analysis
# MacrophageFramework v0.4
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(glmmTMB)
  library(MacrophageFramework)
})

set.seed(1982)

# ------------------------------------------------------------
# 1. Paths and analysis settings
# ------------------------------------------------------------

data_dir <- "/data02/users_space/palma/MacrophageAtlas/datasets/TMS/MatrixMarket"
results_dir <- "/data02/users_space/palma/MacrophageAtlas/datasets/TMS/Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)

tissues <- c(
  "Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Large_Intestine",
  "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Marrow",
  "Pancreas", "Skin", "Spleen", "Thymus", "Tongue", "Trachea"
)

age_map <- c("1m" = 1, "3m" = 3, "18m" = 18, "21m" = 21, "24m" = 24, "30m" = 30)

# ------------------------------------------------------------
# 2. Read each tissue, cluster independently, and detect macrophages
# ------------------------------------------------------------

process_tissue <- function(tissue) {
  tissue_dir <- file.path(data_dir, tissue)

  # Read the sparse count matrix and its feature/cell metadata.
  counts <- readMM(gzfile(file.path(tissue_dir, "matrix.mtx.gz")))
  genes <- read.delim(file.path(tissue_dir, "genes.tsv"),
                      header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(file.path(tissue_dir, "barcodes.tsv"),
                         header = FALSE, stringsAsFactors = FALSE)
  metadata <- read.delim(file.path(tissue_dir, "metadata.tsv"),
                         stringsAsFactors = FALSE, check.names = FALSE)

  gene_names <- make.unique(as.character(genes[[ncol(genes)]]))
  cell_names <- as.character(barcodes[[1]])
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  rownames(metadata) <- as.character(metadata$cell)
  metadata <- metadata[cell_names, , drop = FALSE]

  # Independently cluster each tissue before macrophage detection.
  seu <- CreateSeuratObject(counts = counts, meta.data = metadata,
                            project = tissue, min.cells = 0, min.features = 0)
  seu$tissue <- tissue
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu, features = VariableFeatures(seu))
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30)
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:20, seed.use = 1982)

  # Apply the independent macrophage-core detector at cluster level.
  detector <- analyze_macrophage_profile(
    object = seu, assay = "RNA", layer = "counts", reduction = "umap",
    species = "Mouse", min_core_score = 0.1,
    detect_macrophages = TRUE, save_outputs = FALSE
  )

  keep_clusters <- detector$macrophage_detection$included_clusters
  if (length(keep_clusters) == 0) return(NULL)

  macrophages <- subset(seu, subset = seurat_clusters %in% keep_clusters)
  saveRDS(macrophages, file.path(tissue_dir, paste0(tissue, "_macrophages.rds")))
  macrophages
}

tissue_objects <- lapply(tissues, process_tissue)
names(tissue_objects) <- tissues
tissue_objects <- tissue_objects[!vapply(tissue_objects, is.null, logical(1))]

# ------------------------------------------------------------
# 3. Merge tissue-specific macrophage objects
# ------------------------------------------------------------

tms_macrophages <- merge(
  x = tissue_objects[[1]],
  y = tissue_objects[-1],
  add.cell.ids = names(tissue_objects),
  project = "TMS_macrophages"
)
tms_macrophages$atlas_tissue <- tms_macrophages$tissue

saveRDS(
  tms_macrophages,
  file.path(results_dir, "TMS_macrophages_merged_pre_atlas.rds")
)

# ------------------------------------------------------------
# 4. Global macrophage atlas: normalization, dimensionality reduction,
#    clustering, and visualization
# ------------------------------------------------------------

tms_macrophages <- JoinLayers(tms_macrophages, assay = "RNA")
DefaultAssay(tms_macrophages) <- "RNA"
tms_macrophages <- NormalizeData(
  tms_macrophages, normalization.method = "LogNormalize", scale.factor = 10000
)
tms_macrophages <- FindVariableFeatures(
  tms_macrophages, selection.method = "vst", nfeatures = 3000
)
tms_macrophages <- ScaleData(
  tms_macrophages, features = VariableFeatures(tms_macrophages)
)
tms_macrophages <- RunPCA(
  tms_macrophages, features = VariableFeatures(tms_macrophages), npcs = 50
)
tms_macrophages <- RunUMAP(
  tms_macrophages, dims = 1:30, seed.use = 1982,
  n.neighbors = 30, min.dist = 0.3
)
tms_macrophages <- FindNeighbors(tms_macrophages, dims = 1:30)
tms_macrophages <- FindClusters(
  tms_macrophages, resolution = 0.5, random.seed = 1982
)

saveRDS(
  tms_macrophages,
  file.path(results_dir, "TMS_macrophages_global_atlas.rds")
)

p_umap_clusters <- DimPlot(
  tms_macrophages, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE
) + ggtitle("TMS macrophage atlas — clusters")

p_umap_tissues <- DimPlot(
  tms_macrophages, reduction = "umap", group.by = "atlas_tissue"
) + ggtitle("TMS macrophage atlas — tissues")

ggsave(file.path(results_dir, "Figures", "Figure_1_UMAP_clusters.pdf"),
       p_umap_clusters, width = 8, height = 6)
ggsave(file.path(results_dir, "Figures", "Figure_1_UMAP_clusters.png"),
       p_umap_clusters, width = 8, height = 6, dpi = 300)
ggsave(file.path(results_dir, "Figures", "Figure_2_UMAP_tissues.pdf"),
       p_umap_tissues, width = 8, height = 6)
ggsave(file.path(results_dir, "Figures", "Figure_2_UMAP_tissues.png"),
       p_umap_tissues, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------
# 5. Global MacrophageFramework profiling and cluster annotation
# ------------------------------------------------------------

framework_result <- analyze_macrophage_profile(
  object = tms_macrophages,
  assay = "RNA",
  layer = "counts",
  reduction = "umap",
  species = "Mouse",
  detect_macrophages = FALSE,
  macrophage_clusters = as.character(0:35),
  save_outputs = TRUE,
  output_dir = file.path(results_dir, "MacrophageFramework")
)

framework_annotations <- framework_result$annotations %>%
  transmute(
    Cluster = as.character(Cluster),
    Lineage,
    FunctionalState = Functional_State_1,
    FrameworkAnnotation = paste(Lineage, Functional_State_1, sep = " / ")
  )

framework_meta <- tms_macrophages[[]] %>%
  mutate(Cluster = as.character(seurat_clusters)) %>%
  left_join(framework_annotations, by = "Cluster")

saveRDS(framework_result, file.path(results_dir, "TMS_MacrophageFramework_result.rds"))
write.csv(framework_annotations,
          file.path(results_dir, "TMS_Framework_cluster_annotations.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# 6. Cluster marker analysis and top-five marker visualization
# ------------------------------------------------------------

tms_macrophages <- JoinLayers(tms_macrophages, assay = "RNA")
Idents(tms_macrophages) <- "seurat_clusters"

markers <- FindAllMarkers(
  tms_macrophages,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

write.csv(markers, file.path(results_dir, "TMS_all_cluster_markers.csv"),
          row.names = FALSE)

markers_filtered <- markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25,
    pct.1 > 0.10,
    (pct.1 - pct.2) > 0.10
  )

top_markers <- markers_filtered %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(top_markers,
          file.path(results_dir, "TMS_top5_markers_per_cluster.csv"),
          row.names = FALSE)

top_features <- unique(top_markers$gene)
p_top_markers <- DotPlot(
  tms_macrophages, features = top_features,
  group.by = "seurat_clusters"
) + RotatedAxis() +
  ggtitle("Top five markers per macrophage cluster")

ggsave(file.path(results_dir, "Figures", "Figure_3_Top5_markers.pdf"),
       p_top_markers, width = 14, height = 9)
ggsave(file.path(results_dir, "Figures", "Figure_3_Top5_markers.png"),
       p_top_markers, width = 14, height = 9, dpi = 300)

# ------------------------------------------------------------
# 7. Mouse × tissue × Framework annotation composition
# ------------------------------------------------------------

annotation_order <- c(
  "Resident-associated / Lipid-associated",
  "Monocyte-associated / Lipid-associated",
  "Monocyte-associated / Inflammatory",
  "Monocyte-associated / IFN-responsive"
)

aging_counts <- framework_meta %>%
  count(mouse.id, tissue, age, sex, FrameworkAnnotation, name = "n_cells")

aging_totals <- framework_meta %>%
  count(mouse.id, tissue, age, sex, name = "total_macrophages")

aging_data <- aging_totals %>%
  crossing(FrameworkAnnotation = annotation_order) %>%
  left_join(
    aging_counts,
    by = c("mouse.id", "tissue", "age", "sex", "FrameworkAnnotation")
  ) %>%
  mutate(
    n_cells = replace_na(n_cells, 0L),
    failures = total_macrophages - n_cells,
    age_months = unname(age_map[age])
  )

# Fit the four prespecified composition models from the original aging analysis.
# ECM-remodeling was not part of the original four-category model set.
models_original <- lapply(annotation_order, function(annotation) {
  dat <- aging_data %>%
    filter(FrameworkAnnotation == annotation) %>%
    mutate(
      tissue = factor(tissue),
      sex = factor(sex),
      mouse.id = factor(mouse.id)
    )

  glmmTMB(
    cbind(n_cells, failures) ~ age_months + tissue + sex + (1 | mouse.id),
    family = betabinomial(link = "logit"),
    data = dat
  )
})
names(models_original) <- annotation_order

age_results <- bind_rows(lapply(names(models_original), function(annotation) {
  coef_age <- summary(models_original[[annotation]])$coefficients$cond["age_months", ]
  data.frame(
    Annotation = annotation,
    Beta_age = coef_age["Estimate"],
    SE = coef_age["Std. Error"],
    p_value = coef_age["Pr(>|z|)"],
    convergence = models_original[[annotation]]$fit$convergence,
    positive_hessian = models_original[[annotation]]$sdr$pdHess
  )
})) %>%
  mutate(
    OR_per_10_months = exp(Beta_age * 10),
    CI_low = exp((Beta_age - 1.96 * SE) * 10),
    CI_high = exp((Beta_age + 1.96 * SE) * 10),
    FDR = p.adjust(p_value, method = "BH", n = length(annotation_order))
  )

write.csv(age_results,
          file.path(results_dir, "TMS_aging_composition_model_results.csv"),
          row.names = FALSE)

# Display only estimable, converged, FDR-significant associations.
forest_data <- age_results %>%
  filter(
    convergence == 0,
    positive_hessian,
    is.finite(SE),
    is.finite(CI_low),
    is.finite(CI_high),
    is.finite(FDR),
    FDR < 0.05
  ) %>%
  mutate(Annotation = factor(Annotation, levels = rev(Annotation)))

p_aging_forest <- ggplot(forest_data, aes(x = OR_per_10_months, y = Annotation)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_segment(aes(x = CI_low, xend = CI_high, y = Annotation, yend = Annotation),
               linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_log10() +
  labs(
    title = "Age-associated changes in macrophage composition",
    subtitle = "Odds ratio per 10-month increase in age",
    x = "Odds ratio (95% CI, logarithmic scale)",
    y = NULL
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(results_dir, "Figures", "Figure_G_Aging_composition_forestplot.pdf"),
       p_aging_forest, width = 8, height = 4.5)
ggsave(file.path(results_dir, "Figures", "Figure_G_Aging_composition_forestplot.png"),
       p_aging_forest, width = 8, height = 4.5, dpi = 300)

saveRDS(
  list(
    seurat = tms_macrophages,
    framework = framework_result,
    annotations = framework_annotations,
    aging_data = aging_data,
    aging_models = models_original,
    aging_results = age_results
  ),
  file.path(results_dir, "TMS_macrophage_atlas_complete_analysis.rds")
)
