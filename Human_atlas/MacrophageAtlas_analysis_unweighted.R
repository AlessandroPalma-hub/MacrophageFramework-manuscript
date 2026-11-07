# Macrophage Atlas analysis pipeline
# Reproduces integration, framework summaries, clinical comparisons, and figures.

library(Seurat)
library(harmony)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

devtools::install_github(
  "AlessandroPalma-hub/MacrophageFramework",
  ref = "manuscript-unweighted-v1.0.0",
  upgrade = "never"
)

library(MacrophageFramework)

# Load the merged count matrix and aligned cell metadata.
# Input directory containing the merged count matrix and cell metadata.
DATA_DIR <- "PATH_TO_HUMAN_ATLAS_INPUT_DIRECTORY"

# Output directory for analysis results.
output_dir <- "Human_atlas_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(DATA_DIR, "merged_counts_sparse.rds"))

metadata <- readRDS(file.path(DATA_DIR, "atlas_merged_metadata.rds"))

rownames(metadata) <- sub("^[^.]+[.]", "", rownames(metadata))

metadata <- metadata[colnames(counts), , drop = FALSE]

atlas <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  meta.data = metadata,
  project = "MacrophageFramework",
  min.cells = 0,
  min.features = 0
)

saveRDS(
  atlas,
  file = file.path(output_dir, "atlas_seurat.rds"),
  compress = FALSE
)

# Normalize the atlas, identify variable features, and compute PCA.
DefaultAssay(atlas) <- "RNA"

atlas <- NormalizeData(
  atlas,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

atlas <- FindVariableFeatures(
  atlas,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

atlas <- ScaleData(
  atlas,
  features = VariableFeatures(atlas),
  verbose = FALSE
)

atlas <- RunPCA(
  atlas,
  features = VariableFeatures(atlas),
  npcs = 50,
  verbose = FALSE
)

pdf(file.path(output_dir, "ElbowPlot.pdf"))
print(ElbowPlot(atlas, ndims = 50))
dev.off()

# Correct dataset-associated variation with Harmony and cluster the integrated embedding.
atlas <- harmony::RunHarmony(
  object = atlas,
  group.by.vars = "dataset_id",
  reduction = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  verbose = FALSE
)

atlas <- FindNeighbors(
  atlas,
  reduction = "harmony",
  dims = 1:30,
  k.param = 20,
  verbose = FALSE
)

atlas <- FindClusters(
  atlas,
  resolution = 0.5,
  algorithm = 1,
  random.seed = 1234,
  verbose = FALSE
)

gc(verbose = FALSE)

atlas <- RunUMAP(
  atlas,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 1234,
  verbose = FALSE
)

pdf(file.path(output_dir, "UMAP_clusters_res0.5.pdf"), width = 12, height = 9)
print(DimPlot(
  atlas,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
))
dev.off()

pdf(file.path(output_dir, "UMAP_dataset_id.pdf"), width = 12, height = 9)
print(DimPlot(
  atlas,
  reduction = "umap",
  group.by = "dataset_id",
  raster = TRUE
))
dev.off()

tab_ds_cluster <- table(
  dataset_id = atlas$dataset_id,
  cluster = Idents(atlas)
)

prop_ds_cluster <- prop.table(tab_ds_cluster, margin = 2) * 100

cluster_batch_check <- data.frame(
  cluster = colnames(tab_ds_cluster),
  n_cells = as.integer(colSums(tab_ds_cluster)),
  n_datasets = colSums(prop_ds_cluster >= 5),
  dominant_dataset = rownames(tab_ds_cluster)[
    max.col(t(tab_ds_cluster), ties.method = "first")
  ],
  dominant_dataset_pct = apply(prop_ds_cluster, 2, max),
  row.names = NULL
)

cluster_batch_check$dominant_dataset <- apply(
  prop_ds_cluster, 2,
  function(x) names(which.max(x))
)

cluster_batch_check <- cluster_batch_check[
  order(-cluster_batch_check$dominant_dataset_pct),
]

write.csv(
  round(prop_ds_cluster, 2),
  file.path(output_dir, "dataset_composition_within_clusters.csv"), row.names = TRUE
)

DefaultAssay(atlas) <- "RNA"

Idents(atlas) <- "seurat_clusters"

markers <- FindAllMarkers(
  atlas,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  max.cells.per.ident = 5000,
  random.seed = 1234,
  verbose = FALSE
)

top10 <- do.call(rbind, lapply(
  split(markers, markers$cluster),
  function(x) {
    x <- x[order(x$p_val_adj, -x$avg_log2FC), ]
    head(x, 10)
  }
))

write.csv(markers, file.path(output_dir, "cluster_markers_all.csv"), row.names = FALSE)
write.csv(top10, file.path(output_dir, "cluster_markers_top10.csv"), row.names = FALSE)

checkpoint <- file.path(
  output_dir,
  "atlas_harmony_30PC_res0.5.rds"
)

saveRDS(
  atlas,
  file = checkpoint,
  compress = FALSE
)

analysis_result <- MacrophageFramework::analyze_macrophage_profile(
  object = atlas,
  ontology = load_ontology(),
  assay = "RNA",
  layer = "counts",
  detect_macrophages = FALSE,
  macrophage_clusters = as.character(unique(Seurat::Idents(atlas))),
  min_core_score = 0.1,
  species = "Human",
  output_dir = output_dir,
  save_outputs = TRUE
)

# -----------------------------------------------------------------------------
# UNWEIGHTED HUMAN ATLAS: sample-level summaries and manuscript figures
# The preprocessing and MacrophageFramework call above are retained unchanged.
# -----------------------------------------------------------------------------

functional_programs <- c(
  "Inflammatory", "Repair-associated", "Lipid-associated",
  "ECM-remodeling", "IFN-responsive"
)

# Cell-level scores -> sample-by-cluster means -> sample-level weighted means.
score_mat <- analysis_result$scores
md_use <- atlas@meta.data[colnames(score_mat), , drop = FALSE]

cell_scores <- as.data.frame(t(score_mat), check.names = FALSE)
cell_scores$dataset <- as.character(md_use$dataset_id)
cell_scores$sample <- as.character(md_use$sample_id)
cell_scores$tissue <- as.character(md_use$tissue)
cell_scores$cluster <- as.character(Seurat::Idents(atlas)[colnames(score_mat)])

programs <- rownames(score_mat)
score_columns <- programs

sample_cluster_scores <- cell_scores %>%
  group_by(dataset, sample, tissue, cluster) %>%
  summarise(
    across(all_of(score_columns), ~ mean(.x, na.rm = TRUE)),
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  sample_cluster_scores,
  file.path(output_dir, "17_sample_cluster_program_scores_unweighted.csv"),
  row.names = FALSE
)

sample_scores <- sample_cluster_scores %>%
  mutate(across(all_of(programs), ~ .x * n_cells, .names = "weighted_{.col}")) %>%
  group_by(dataset, sample, tissue) %>%
  summarise(
    n_cells = sum(n_cells),
    across(starts_with("weighted_"), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

for (program in programs) {
  sample_scores[[program]] <-
    sample_scores[[paste0("weighted_", program)]] / sample_scores$n_cells
}
sample_scores <- sample_scores %>%
  select(dataset, sample, tissue, n_cells, all_of(programs))

write.csv(
  sample_scores,
  file.path(output_dir, "18_sample_level_program_scores_unweighted.csv"),
  row.names = FALSE
)

# Panel B: sample-level scores across tissues.
panel_B_data <- sample_scores %>%
  select(dataset, sample, tissue, n_cells, all_of(functional_programs)) %>%
  pivot_longer(
    cols = all_of(functional_programs),
    names_to = "Program", values_to = "Score"
  )

p_panel_B <- ggplot(panel_B_data, aes(x = tissue, y = Score)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey35", width = 0.65) +
  geom_jitter(width = 0.14, alpha = 0.65, size = 1.5) +
  facet_wrap(~Program, scales = "free_y", ncol = 2) +
  labs(
    title = "Sample-level macrophage program scores across tissues",
    subtitle = "Each point represents one sample",
    x = NULL, y = "Mean program score"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1),
        strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "Panel_B_sample_level_program_scores_unweighted.pdf"),
       p_panel_B, width = 12, height = 10)
ggsave(file.path(output_dir, "Panel_B_sample_level_program_scores_unweighted.png"),
       p_panel_B, width = 12, height = 10, dpi = 300)

# Panel C: tissue-to-annotation circos, based on sample-level annotation proportions.
ann <- analysis_result$annotations
annotation_map <- setNames(as.character(ann$Full_annotation), as.character(ann$Cluster))
md_all <- atlas@meta.data
cell_annotation <- data.frame(
  dataset = as.character(md_all$dataset_id),
  sample = as.character(md_all$sample_id),
  tissue = as.character(md_all$tissue),
  cluster = as.character(Seurat::Idents(atlas)),
  stringsAsFactors = FALSE
)
cell_annotation$annotation <- unname(annotation_map[cell_annotation$cluster])
cell_annotation <- cell_annotation %>%
  filter(!is.na(dataset), !is.na(sample), !is.na(tissue), !is.na(annotation)) %>%
  mutate(sample_key = paste(dataset, sample, sep = "__"))

annotation_counts <- cell_annotation %>%
  count(sample_key, tissue, annotation, name = "n_cells")
annotation_props <- annotation_counts %>%
  group_by(sample_key, tissue) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()
tissue_annotation <- annotation_props %>%
  group_by(tissue, annotation) %>%
  summarise(fraction = mean(fraction), .groups = "drop")

annotations <- sort(unique(tissue_annotation$annotation))
annotation_codes <- setNames(sprintf("A%02d", seq_along(annotations)), annotations)
annotation_colors <- setNames(grDevices::hcl.colors(length(annotations), "Dark 3"), annotations)
tissue_nodes <- sort(unique(tissue_annotation$tissue))
annotation_nodes <- unname(annotation_codes[annotations])
links <- tissue_annotation %>%
  mutate(annotation_code = unname(annotation_codes[annotation]),
         col = grDevices::adjustcolor(annotation_colors[annotation], alpha.f = 0.48))

pdf(file.path(output_dir, "Panel_C_circos_tissues_annotations_unweighted.pdf"),
    width = 16, height = 16, useDingbats = FALSE)
circos.clear()
circos.par(start.degree = 90, gap.degree = 3, points.overflow.warning = FALSE)
chordDiagram(
  data.frame(from = links$tissue, to = links$annotation_code, value = links$fraction),
  order = c(tissue_nodes, annotation_nodes),
  grid.col = c(setNames(grDevices::hcl.colors(length(tissue_nodes), "Set 2"), tissue_nodes),
               setNames(annotation_colors[annotations], annotation_nodes)),
  col = links$col, transparency = 0.15,
  annotationTrack = "grid", preAllocateTracks = list(track.height = 0.13)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  circos.text(mean(xlim), 0.5, sector, facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = if (sector %in% tissue_nodes) 0.72 else 0.85)
}, bg.border = NA)
title("Tissue distribution of MacrophageFramework annotations")
circos.clear()
plot.new(); par(mar = c(1, 1, 3, 1))
title(main = "MacrophageFramework annotation legend", cex.main = 1.2)
legend("topleft", legend = paste0(unname(annotation_codes[annotations]), " | ", annotations),
       col = unname(annotation_colors[annotations]), pch = 15, pt.cex = 1.3,
       cex = 0.78, bty = "n", x.intersp = 0.7, y.intersp = 1.35)
dev.off()

# Definitive clinical grouping used for Panels D and E.
# Keep all accepted samples; exclude only ambiguous hashed mixed samples.
global_audit <- unique(atlas@meta.data[, c(
  "dataset_id", "sample_id", "condition", "condition_detail"
)])
names(global_audit) <- c("dataset", "sample", "condition", "condition_detail")
global_audit$global_group <- NA_character_

set_group <- function(ds, pattern, group, field = "condition_detail") {
  idx <- global_audit$dataset == ds &
    grepl(pattern, ifelse(is.na(global_audit[[field]]), "", global_audit[[field]]),
          ignore.case = TRUE)
  global_audit$global_group[idx] <<- group
}

# Healthy / non-tumor pathology assignments.
set_group("GSE128033", "^Normal$", "Healthy", "condition")
set_group("GSE128033", "^IPF$", "Non-tumor pathology", "condition")
set_group("GSE136103", "healthy", "Healthy")
set_group("GSE136103", "cirrhosis", "Non-tumor pathology")
set_group("GSE145154", "^Control$", "Healthy", "condition")
set_group("GSE145154", "^(DCM|ICM)$", "Non-tumor pathology", "condition")
set_group("GSE145926", "^control$", "Healthy", "condition")
set_group("GSE145926", "^(moderate|severe)$", "Non-tumor pathology", "condition")
set_group("GSE155960", "^lean$", "Healthy", "condition")
set_group("GSE155960", "^obese$", "Non-tumor pathology", "condition")
set_group("GSE156110", "lean", "Healthy")
set_group("GSE156110", "obese", "Non-tumor pathology")
set_group("GSE162183", "^(Healthy|Ctrl)", "Healthy", "condition")
set_group("GSE162183", "Psoriasis", "Non-tumor pathology", "condition")
set_group("GSE165045", "^Control$", "Healthy", "condition")
set_group("GSE165045", "^(Hereditary|Idiopathic)$", "Non-tumor pathology", "condition")
set_group("GSE196964", "^Normal$", "Healthy", "condition")
set_group("GSE196964", "^Tumor$", "Tumor", "condition")
set_group("GSE202109", "^Healthy$", "Healthy", "condition")
set_group("GSE235143", "^Control$", "Healthy", "condition")
set_group("GSE235143", "^PAD$", "Non-tumor pathology", "condition")
set_group("GSE252126", "^control$", "Healthy", "condition")
set_group("GSE252126", "gastroparesis", "Non-tumor pathology", "condition")

# Tumor datasets.
set_group("GSE141383", "GBM", "Tumor", "condition")
set_group("GSE224392", "HGSOC", "Tumor", "condition")
set_group("GSE261774", "breast_cancer", "Tumor", "condition")

global_audit$global_group[
  global_audit$dataset == "GSE196964" &
    grepl("hashed", global_audit$condition_detail, ignore.case = TRUE)
] <- NA_character_

global_primary <- global_audit %>%
  filter(!is.na(global_group)) %>%
  mutate(global_group = factor(global_group,
    levels = c("Healthy", "Non-tumor pathology", "Tumor")))

# Join clinical groups to sample-level program scores using dataset + sample.
sample_scores <- sample_scores %>%
  left_join(global_primary %>% select(dataset, sample, global_group),
            by = c("dataset", "sample"))

# Panel D: mean within-sample relative composition of the five programs.
stacked_data <- sample_scores %>%
  filter(!is.na(global_group)) %>%
  select(dataset, sample, global_group, all_of(functional_programs)) %>%
  pivot_longer(cols = all_of(functional_programs),
               names_to = "Program", values_to = "Score") %>%
  group_by(dataset, sample, global_group) %>%
  mutate(Relative_score = Score / sum(Score, na.rm = TRUE) * 100) %>%
  ungroup()
stacked_summary <- stacked_data %>%
  group_by(global_group, Program) %>%
  summarise(Mean_relative_score = mean(Relative_score, na.rm = TRUE), .groups = "drop")
p_panel_D <- ggplot(stacked_summary,
  aes(x = global_group, y = Mean_relative_score, fill = Program)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.35) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Relative composition of functional macrophage programs",
       subtitle = "Mean within-sample normalized program scores",
       x = NULL, y = "Relative program score (%)", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(colour = "grey35"), legend.position = "right")
ggsave(file.path(output_dir, "Panel_D_clinical_groups_unweighted.pdf"), p_panel_D, width = 8, height = 5.5)
ggsave(file.path(output_dir, "Panel_D_clinical_groups_unweighted.png"), p_panel_D, width = 8, height = 5.5, dpi = 300)

# Panel E: descriptive distributions; no inferential statistics.
panel_E_data <- sample_scores %>%
  filter(!is.na(global_group)) %>%
  select(dataset, sample, global_group, all_of(functional_programs)) %>%
  pivot_longer(cols = all_of(functional_programs), names_to = "Program", values_to = "Score")
clinical_colors <- c("Healthy" = "#9BB5D1", "Non-tumor pathology" = "#E8C77E", "Tumor" = "#E9A18F")
p_panel_E <- ggplot(panel_E_data, aes(x = global_group, y = Score, fill = global_group)) +
  geom_violin(trim = TRUE, scale = "width", color = "#333333", linewidth = 0.6) +
  geom_jitter(width = 0.08, size = 1.8, alpha = 0.65, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", color = "black") +
  facet_wrap(~Program, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = clinical_colors) +
  labs(title = "Functional macrophage programs across clinical groups",
       subtitle = "Integrated atlas; descriptive comparison across heterogeneous tissues and datasets",
       x = NULL, y = "Program score", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(output_dir, "Panel_E_functional_programs_violin_unweighted.pdf"), p_panel_E, width = 12, height = 8)
ggsave(file.path(output_dir, "Panel_E_functional_programs_violin_unweighted.png"), p_panel_E, width = 12, height = 8, dpi = 300)

# IPF panel: normal lung tissue versus IPF, excluding the two normal BAL samples.
ipf_data <- sample_scores %>%
  filter(dataset == "GSE128033") %>%
  left_join(global_audit %>%
      filter(dataset == "GSE128033") %>%
      select(dataset, sample, condition_detail),
    by = c("dataset", "sample")) %>%
  filter(condition_detail != "normal_BAL") %>%
  mutate(Group = ifelse(grepl("^IPF", condition_detail), "IPF", "Normal lung"),
         Group = factor(Group, levels = c("Normal lung", "IPF")))

ipf_long <- ipf_data %>%
  select(sample, Group, all_of(functional_programs)) %>%
  pivot_longer(cols = all_of(functional_programs), names_to = "Program", values_to = "Score")

ipf_stats <- ipf_long %>%
  group_by(Program) %>%
  summarise(
    n_normal = sum(Group == "Normal lung"),
    n_ipf = sum(Group == "IPF"),
    median_normal = median(Score[Group == "Normal lung"], na.rm = TRUE),
    median_ipf = median(Score[Group == "IPF"], na.rm = TRUE),
    p_value = wilcox.test(Score ~ Group, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))
write.csv(ipf_stats, file.path(output_dir, "Panel_IPF_functional_programs_statistics_unweighted.csv"), row.names = FALSE)

ipf_annotations <- ipf_stats %>%
  mutate(label = case_when(
    FDR < 0.001 ~ paste0("FDR < 0.001 ***"),
    FDR < 0.01 ~ paste0("FDR = ", sprintf("%.3f", FDR), " **"),
    FDR < 0.05 ~ paste0("FDR = ", sprintf("%.3f", FDR), " *"),
    TRUE ~ paste0("FDR = ", sprintf("%.2f", FDR), " ns")
  )) %>%
  left_join(ipf_long %>% group_by(Program) %>%
      summarise(y = max(Score, na.rm = TRUE) + 0.13 * diff(range(Score, na.rm = TRUE)),
                .groups = "drop"), by = "Program")

p_ipf <- ggplot(ipf_long, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "#333333", linewidth = 0.6) +
  geom_jitter(width = 0.08, size = 2.5, alpha = 0.9, color = "black") +
  facet_wrap(~Program, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Normal lung" = "#B5C9E0", "IPF" = "#E8B0A1")) +
  labs(title = "Functional macrophage programs in IPF",
       subtitle = "Human lung tissue: IPF versus normal lung",
       x = NULL, y = "Program score", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 13),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "top", axis.text = element_text(color = "black"))

p_ipf_final <- p_ipf + geom_text(data = ipf_annotations,
  aes(x = 1.5, y = y, label = label), inherit.aes = FALSE,
  fontface = "bold", size = 4)
ggsave(file.path(output_dir, "Panel_IPF_functional_programs_unweighted.pdf"), p_ipf_final, width = 12, height = 8, device = cairo_pdf)
ggsave(file.path(output_dir, "Panel_IPF_functional_programs_unweighted.png"), p_ipf_final, width = 12, height = 8, dpi = 300)

# Save the final sample-level table used by the manuscript panels.
write.csv(sample_scores, file.path(output_dir, "sample_level_program_scores_clinical_unweighted.csv"), row.names = FALSE)
