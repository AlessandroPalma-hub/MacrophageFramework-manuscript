# ============================================================
# MacrophageFramework — sequencing-depth sensitivity audit
# Dataset: Human ovarian (GSE224392)
#
# Uses the existing baseline_ovarian object.
# Does NOT rerun the baseline and does NOT copy the Seurat object.
# ============================================================

# ----------------------------
# 1. Settings
# ----------------------------

# Load the MacrophageFramework release used in the manuscript.
devtools::install_github(
  "AlessandroPalma-hub/MacrophageFramework",
  ref = "manuscript-unweighted-v1.0.0",
  upgrade = "never"
)

library(MacrophageFramework)

OUTPUT_DIR <- "SENSITIVITY/OVARIAN_SEQUENCING_DEPTH_SENSITIVITY"
FRACTIONS <- c(0.75, 0.50, 0.25)
N_REPLICATES <- 10L
BASE_SEED <- 20261101L
ASSAY <- "RNA"
COUNTS_LAYER <- "counts"
CLUSTER_COLUMN <- "seurat_clusters"
MIN_CORE_SCORE <- 0.1
SPECIES <- "Human"

# ----------------------------
# 2. Validate session
# ----------------------------
required_objects <- c("ovarian", "baseline_ovarian", "ontology_v04")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects) > 0L) {
  stop("Missing objects: ", paste(missing_objects, collapse = ", "))
}
if (!inherits(ovarian, "Seurat")) stop("`ovarian` is not a Seurat object.")
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("Package 'SeuratObject' is required.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required.")
}
if (!requireNamespace("MacrophageFramework", quietly = TRUE)) {
  stop("Package 'MacrophageFramework' is not installed.")
}

counts_original <- SeuratObject::LayerData(
  ovarian, assay = ASSAY, layer = COUNTS_LAYER
)
if (!inherits(counts_original, "sparseMatrix")) {
  stop("The RNA counts layer must be a sparse Matrix object.")
}
if (is.null(rownames(counts_original)) || is.null(colnames(counts_original))) {
  stop("The counts matrix must have gene and cell names.")
}
if (any(counts_original@x < 0, na.rm = TRUE)) {
  stop("Negative values detected in counts.")
}
if (any(abs(counts_original@x - round(counts_original@x)) > 1e-8)) {
  stop("Counts are not integer-valued; binomial thinning requires integer counts.")
}

clusters_all <- as.character(ovarian[[CLUSTER_COLUMN, drop = TRUE]])
names(clusters_all) <- colnames(ovarian)
clusters_all <- clusters_all[colnames(counts_original)]

if (anyNA(clusters_all)) stop("Missing cluster assignments detected.")
if (!identical(names(clusters_all), colnames(counts_original))) {
  stop("Cluster vector could not be aligned to counts columns.")
}

if (is.null(baseline_ovarian$scores) ||
    is.null(baseline_ovarian$annotations) ||
    is.null(baseline_ovarian$analyzed_cells) ||
    is.null(baseline_ovarian$analyzed_clusters)) {
  stop("baseline_ovarian is missing required result components.")
}

baseline_cells <- as.character(baseline_ovarian$analyzed_cells)
if (!all(baseline_cells %in% colnames(counts_original))) {
  stop("Some baseline analyzed_cells are absent from the original counts.")
}

baseline_clusters <- as.character(baseline_ovarian$analyzed_clusters)
if (!all(baseline_clusters %in% unique(clusters_all))) {
  stop("Some baseline analyzed_clusters are absent from the original metadata.")
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
LOG_FILE <- file.path(OUTPUT_DIR, "audit_log.txt")
writeLines(
  c(
    "MacrophageFramework sequencing-depth sensitivity audit",
    paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "Dataset: Human ovarian GSE224392",
    "Baseline: existing baseline_ovarian object",
    paste("Fractions:", paste(FRACTIONS, collapse = ", ")),
    paste("Replicates per fraction:", N_REPLICATES),
    paste("Base seed:", BASE_SEED),
    ""
  ),
  LOG_FILE
)

log_message <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(
    paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n"),
    file = LOG_FILE, append = TRUE
  )
}

# ----------------------------
# 3. Helper functions
# ----------------------------
downsample_counts <- function(count_matrix, fraction, seed) {
  set.seed(seed)
  out <- count_matrix
  out@x <- stats::rbinom(
    n = length(count_matrix@x),
    size = as.integer(round(count_matrix@x)),
    prob = fraction
  )
  storage.mode(out@x) <- "double"
  Matrix::drop0(out)
}

run_framework_matrix <- function(expression_matrix, cluster_vector) {
  MacrophageFramework::analyze_macrophage_profile(
    expression = expression_matrix,
    clusters = unname(cluster_vector),
    ontology = ontology_v04,
    detect_macrophages = TRUE,
    min_core_score = MIN_CORE_SCORE,
    species = SPECIES,
    save_outputs = FALSE
  )
}

score_to_long <- function(score_matrix, run_id, fraction, replicate, seed) {
  if (is.null(score_matrix) || length(score_matrix) == 0L) return(NULL)
  score_matrix <- as.matrix(score_matrix)
  if (is.null(rownames(score_matrix)) || is.null(colnames(score_matrix))) {
    stop("Score matrix must have program rownames and cell colnames.")
  }
  data.frame(
    run_id = run_id,
    fraction = fraction,
    replicate = replicate,
    seed = seed,
    cell = rep(colnames(score_matrix), each = nrow(score_matrix)),
    program = rep(rownames(score_matrix), times = ncol(score_matrix)),
    score = as.numeric(score_matrix),
    stringsAsFactors = FALSE
  )
}

annotation_to_long <- function(ann, run_id, fraction, replicate, seed) {
  if (is.null(ann) || nrow(ann) == 0L) return(NULL)
  ann <- as.data.frame(ann, stringsAsFactors = FALSE)
  if (!"Cluster" %in% names(ann)) {
    ann$Cluster <- rownames(ann)
  }
  if (!"Cluster" %in% names(ann)) {
    stop("Annotation table has no Cluster column or rownames.")
  }
  cbind(
    data.frame(
      run_id = run_id, fraction = fraction,
      replicate = replicate, seed = seed,
      stringsAsFactors = FALSE
    ),
    ann
  )
}

detector_to_long <- function(detector, run_id, fraction, replicate, seed) {
  if (is.null(detector) || is.null(detector$cluster_metrics)) return(NULL)
  metrics <- as.data.frame(detector$cluster_metrics, stringsAsFactors = FALSE)
  if (nrow(metrics) == 0L) return(NULL)
  if (!"Cluster" %in% names(metrics)) {
    metrics$Cluster <- rownames(metrics)
  }
  cbind(
    data.frame(
      run_id = run_id, fraction = fraction,
      replicate = replicate, seed = seed,
      stringsAsFactors = FALSE
    ),
    metrics
  )
}

compare_scores <- function(reference, current, run_id, fraction, replicate, seed) {
  ref <- as.matrix(reference$scores)
  cur <- as.matrix(current$scores)
  common_programs <- intersect(rownames(ref), rownames(cur))
  common_cells <- intersect(colnames(ref), colnames(cur))
  if (length(common_cells) < 3L || length(common_programs) == 0L) return(NULL)

  rows <- lapply(common_programs, function(program) {
    x <- as.numeric(ref[program, common_cells])
    y <- as.numeric(cur[program, common_cells])
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) == 0L) return(NULL)
    abs_diff <- abs(y - x)
    denom <- stats::median(abs(x), na.rm = TRUE)

    data.frame(
      run_id = run_id,
      fraction = fraction,
      replicate = replicate,
      seed = seed,
      program = program,
      n_cells_compared = length(x),
      pearson_r = if (
        length(x) > 2L && stats::sd(x) > 0 && stats::sd(y) > 0
      ) stats::cor(x, y, method = "pearson") else NA_real_,
      spearman_rho = if (length(x) > 2L) {
        suppressWarnings(stats::cor(x, y, method = "spearman"))
      } else {
        NA_real_
      },
      median_absolute_change = stats::median(abs_diff, na.rm = TRUE),
      median_relative_change_pct = if (is.finite(denom) && denom > 0) {
        100 * stats::median(abs_diff, na.rm = TRUE) / denom
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

compare_annotations <- function(reference, current, run_id, fraction, replicate, seed) {
  ref <- as.data.frame(reference$annotations, stringsAsFactors = FALSE)
  cur <- as.data.frame(current$annotations, stringsAsFactors = FALSE)
  fields <- c("Lineage", "Functional_State_1", "Full_annotation")

  if (is.null(ref) || is.null(cur) || nrow(ref) == 0L || nrow(cur) == 0L) {
    return(NULL)
  }
  if (!"Cluster" %in% names(ref)) ref$Cluster <- rownames(ref)
  if (!"Cluster" %in% names(cur)) cur$Cluster <- rownames(cur)
  if (!"Cluster" %in% names(ref) || !"Cluster" %in% names(cur)) {
    stop("Cannot compare annotations: Cluster identifier is missing.")
  }

  ref$Cluster <- as.character(ref$Cluster)
  cur$Cluster <- as.character(cur$Cluster)
  joined <- merge(ref, cur, by = "Cluster", suffixes = c("_100", "_current"), all = TRUE)

  detail <- data.frame(
    run_id = run_id,
    fraction = fraction,
    replicate = replicate,
    seed = seed,
    Cluster = joined$Cluster,
    stringsAsFactors = FALSE
  )

  summary_rows <- lapply(fields, function(field) {
    ref_col <- paste0(field, "_100")
    cur_col <- paste0(field, "_current")
    if (!all(c(ref_col, cur_col) %in% names(joined))) {
      return(data.frame(
        run_id = run_id, fraction = fraction, replicate = replicate,
        seed = seed, field = field, n_reference_clusters = nrow(ref),
        n_current_clusters = nrow(cur), n_clusters_compared = 0L,
        n_matches = 0L, concordance_pct = NA_real_,
        annotation_retention_pct = 0,
        stringsAsFactors = FALSE
      ))
    }
    a <- as.character(joined[[ref_col]])
    b <- as.character(joined[[cur_col]])
    comparable <- !is.na(a) & !is.na(b)
    matched <- comparable & a == b

    detail[[ref_col]] <<- a
    detail[[cur_col]] <<- b
    detail[[paste0(field, "_match")]] <<- matched

    data.frame(
      run_id = run_id,
      fraction = fraction,
      replicate = replicate,
      seed = seed,
      field = field,
      n_reference_clusters = nrow(ref),
      n_current_clusters = nrow(cur),
      n_clusters_compared = sum(comparable),
      n_matches = sum(matched),
      concordance_pct = if (sum(comparable) > 0L) {
        100 * sum(matched) / sum(comparable)
      } else {
        NA_real_
      },
      annotation_retention_pct = if (nrow(ref) > 0L) {
        100 * sum(matched) / nrow(ref)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })

  list(detail = detail, summary = do.call(rbind, summary_rows))
}

bind_nonempty <- function(x) {
  x <- Filter(function(z) !is.null(z) && NROW(z) > 0L, x)
  if (length(x) == 0L) return(data.frame())
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

# ----------------------------
# 4. Initialize with existing baseline
# ----------------------------
reference_result <- baseline_ovarian
log_message(
  "Using existing baseline_ovarian: ",
  length(reference_result$analyzed_cells), " cells; clusters ",
  paste(reference_result$analyzed_clusters, collapse = ", ")
)

all_scores <- list(
  score_to_long(reference_result$scores, "depth100", 1, 0, NA_integer_)
)
all_annotations <- list(
  annotation_to_long(reference_result$annotations, "depth100", 1, 0, NA_integer_)
)
all_detector <- list(
  detector_to_long(
    reference_result$macrophage_detection, "depth100", 1, 0, NA_integer_
  )
)
lib0 <- Matrix::colSums(counts_original)
all_depth <- list(data.frame(
  run_id = "depth100",
  fraction = 1,
  replicate = 0,
  seed = NA_integer_,
  n_genes = nrow(counts_original),
  n_cells = ncol(counts_original),
  total_counts = sum(counts_original),
  mean_library_size = mean(lib0),
  median_library_size = stats::median(lib0),
  min_library_size = min(lib0),
  max_library_size = max(lib0),
  stringsAsFactors = FALSE
))

all_score_comparisons <- list()
all_annotation_details <- list()
all_annotation_summaries <- list()
all_errors <- list()

# ----------------------------
# 5. Downsampling audit
# ----------------------------
seeds <- BASE_SEED + seq_len(N_REPLICATES) - 1L
total_runs <- length(FRACTIONS) * N_REPLICATES
run_counter <- 0L

for (fraction in FRACTIONS) {
  for (i in seq_along(seeds)) {
    seed <- seeds[i]
    run_counter <- run_counter + 1L
    run_id <- paste0(
      "depth", sprintf("%02d", round(100 * fraction)),
      "_rep", sprintf("%02d", i)
    )

    log_message(
      "[", run_counter, "/", total_runs, "] ",
      run_id, " — seed ", seed
    )

    tryCatch({
      down_counts <- downsample_counts(counts_original, fraction, seed)

      result <- run_framework_matrix(
        expression_matrix = down_counts,
        cluster_vector = clusters_all
      )

      all_scores[[length(all_scores) + 1L]] <- score_to_long(
        result$scores, run_id, fraction, i, seed
      )
      all_annotations[[length(all_annotations) + 1L]] <- annotation_to_long(
        result$annotations, run_id, fraction, i, seed
      )
      all_detector[[length(all_detector) + 1L]] <- detector_to_long(
        result$macrophage_detection, run_id, fraction, i, seed
      )

      lib <- Matrix::colSums(down_counts)
      all_depth[[length(all_depth) + 1L]] <- data.frame(
        run_id = run_id,
        fraction = fraction,
        replicate = i,
        seed = seed,
        n_genes = nrow(down_counts),
        n_cells = ncol(down_counts),
        total_counts = sum(down_counts),
        mean_library_size = mean(lib),
        median_library_size = stats::median(lib),
        min_library_size = min(lib),
        max_library_size = max(lib),
        stringsAsFactors = FALSE
      )

      score_cmp <- compare_scores(
        reference_result, result, run_id, fraction, i, seed
      )
      if (!is.null(score_cmp)) {
        all_score_comparisons[[length(all_score_comparisons) + 1L]] <- score_cmp
      }

      ann_cmp <- compare_annotations(
        reference_result, result, run_id, fraction, i, seed
      )
      if (!is.null(ann_cmp)) {
        all_annotation_details[[length(all_annotation_details) + 1L]] <- ann_cmp$detail
        all_annotation_summaries[[length(all_annotation_summaries) + 1L]] <- ann_cmp$summary
      }

      rm(down_counts, result)
      invisible(gc())
      log_message("Completed: ", run_id)

    }, error = function(e) {
      error_entry <- data.frame(
        run_id = run_id,
        fraction = fraction,
        replicate = i,
        seed = seed,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      all_errors[[length(all_errors) + 1L]] <<- error_entry
      log_message("ERROR in ", run_id, ": ", conditionMessage(e))
    })
  }
}

# ----------------------------
# 6. Combine and save outputs
# ----------------------------
scores_df <- bind_nonempty(all_scores)
annotations_df <- bind_nonempty(all_annotations)
detector_df <- bind_nonempty(all_detector)
depth_df <- bind_nonempty(all_depth)
score_comparison_df <- bind_nonempty(all_score_comparisons)
annotation_detail_df <- bind_nonempty(all_annotation_details)
annotation_summary_df <- bind_nonempty(all_annotation_summaries)
errors_df <- bind_nonempty(all_errors)

utils::write.csv(scores_df, file.path(OUTPUT_DIR, "scores_by_cell.csv"), row.names = FALSE)
utils::write.csv(annotations_df, file.path(OUTPUT_DIR, "annotations_by_cluster.csv"), row.names = FALSE)
utils::write.csv(detector_df, file.path(OUTPUT_DIR, "detector_metrics_by_cluster.csv"), row.names = FALSE)
utils::write.csv(depth_df, file.path(OUTPUT_DIR, "library_depth_by_run.csv"), row.names = FALSE)
utils::write.csv(score_comparison_df, file.path(OUTPUT_DIR, "score_stability_vs_100pct.csv"), row.names = FALSE)
utils::write.csv(annotation_detail_df, file.path(OUTPUT_DIR, "annotation_concordance_detail.csv"), row.names = FALSE)
utils::write.csv(annotation_summary_df, file.path(OUTPUT_DIR, "annotation_concordance_summary.csv"), row.names = FALSE)
utils::write.csv(errors_df, file.path(OUTPUT_DIR, "errors.csv"), row.names = FALSE)

# Replicate-level summary by depth and program.
if (nrow(score_comparison_df) > 0L) {
  grouping <- split(
    score_comparison_df,
    interaction(
      score_comparison_df$fraction,
      score_comparison_df$program,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  summary_rows <- lapply(grouping, function(d) {
    out <- data.frame(
      fraction = d$fraction[1],
      program = d$program[1],
      n_replicates = nrow(d),
      stringsAsFactors = FALSE
    )
    for (metric in c(
      "pearson_r", "spearman_rho",
      "median_absolute_change", "median_relative_change_pct"
    )) {
      x <- d[[metric]]
      out[[paste0(metric, "_mean")]] <- if (all(is.na(x))) {
        NA_real_
      } else {
        mean(x, na.rm = TRUE)
      }
      out[[paste0(metric, "_sd")]] <- if (sum(is.finite(x)) > 1L) {
        stats::sd(x, na.rm = TRUE)
      } else {
        NA_real_
      }
      out[[paste0(metric, "_median")]] <- if (all(is.na(x))) {
        NA_real_
      } else {
        stats::median(x, na.rm = TRUE)
      }
    }
    out
  })
  score_summary_df <- do.call(rbind, summary_rows)
  rownames(score_summary_df) <- NULL
} else {
  score_summary_df <- data.frame()
}
utils::write.csv(
  score_summary_df,
  file.path(OUTPUT_DIR, "score_stability_summary_by_depth.csv"),
  row.names = FALSE
)

# Configuration and reproducibility information.
config <- data.frame(
  parameter = c(
    "dataset", "assay", "counts_layer", "cluster_column", "ontology",
    "species", "min_core_score", "fractions", "replicates_per_fraction",
    "base_seed", "seed_sequence", "n_genes", "n_cells",
    "baseline_analyzed_cells", "baseline_analyzed_clusters"
  ),
  value = c(
    "Human ovarian GSE224392", ASSAY, COUNTS_LAYER, CLUSTER_COLUMN,
    "v0.4", SPECIES, as.character(MIN_CORE_SCORE),
    paste(FRACTIONS, collapse = ";"),
    as.character(N_REPLICATES), as.character(BASE_SEED),
    paste(seeds, collapse = ";"),
    as.character(nrow(counts_original)),
    as.character(ncol(counts_original)),
    as.character(length(baseline_ovarian$analyzed_cells)),
    paste(baseline_ovarian$analyzed_clusters, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  config, file.path(OUTPUT_DIR, "audit_configuration.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)

saveRDS(
  list(
    configuration = config,
    baseline_scores = reference_result$scores,
    baseline_annotations = reference_result$annotations,
    baseline_detector = reference_result$macrophage_detection$cluster_metrics,
    score_stability = score_comparison_df,
    score_stability_summary = score_summary_df,
    annotation_concordance = annotation_summary_df,
    errors = errors_df
  ),
  file.path(OUTPUT_DIR, "audit_summary_objects.rds")
)

# ----------------------------
# 7. Final report
# ----------------------------
successful_runs <- length(all_depth) - 1L
log_message("Audit finished.")
log_message("Successful downsampled runs: ", successful_runs, " / ", total_runs)
log_message("Errors recorded: ", nrow(errors_df))
log_message("Output directory: ", normalizePath(OUTPUT_DIR))

cat("\nSEQUENCING-DEPTH AUDIT FINISHED\n")
cat("Successful runs:", successful_runs, "/", total_runs, "\n")
cat("Errors:", nrow(errors_df), "\n")
cat("Output directory:", normalizePath(OUTPUT_DIR), "\n")

if (successful_runs != total_runs || nrow(errors_df) > 0L) {
  warning("Not all runs completed successfully. Inspect errors.csv and audit_log.txt.")
} else {
  message("All 30 downsampling runs completed successfully.")
}
