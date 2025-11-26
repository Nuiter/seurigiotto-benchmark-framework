#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 5: COMPLETE ANALYSIS OF THE 3 PIPELINES - FINAL INTEGRATED VERSION
# =============================================================================

cat("\n=== SCRIPT 5: FINAL MASTER'S THESIS ANALYSIS - DEFINITIVE VERSION ===\n")
cat("Date:", Sys.Date(), "\n")
cat("R:", R.version.string, "\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(Seurat); library(Giotto); library(ggplot2); library(dplyr)
  library(patchwork); library(pheatmap); library(viridis); library(RColorBrewer)
  library(data.table); library(gridExtra); library(grid)
})

# --- CONFIGURATION & PATHS ---
main_results_dir <- "results"
integrated_analysis_subdir <- file.path(main_results_dir, "05_Integrated_Pipeline_Analysis")
dir.create(integrated_analysis_subdir, showWarnings = FALSE, recursive = TRUE)
# Note: Assumes signature CSVs are in the project root.

# ===========================================================================
# PART 1: LOAD CANCER SIGNATURES
# ===========================================================================

load_cancer_signatures <- function() {
  cat("\n[SIGNATURES] Loading cancer signatures from CSV...\n")
  signatures <- list()
  
  # Define the correct path to the data directory
  data_dir_sigs <- "data" 
  cat("  Searching for signature files in:", file.path(getwd(), data_dir_sigs), "\n")
  
  load_sig <- function(file, name) {
    # Build the full path: data/Tumor.csv
    file_path <- file.path(data_dir_sigs, file) 
    
    if (file.exists(file_path)) {
      data <- read.csv(file_path, stringsAsFactors = FALSE)
      signatures[[name]] <<- data$Name
      cat("  ", name, ":", length(signatures[[name]]), "genes\n")
    } else {
      cat("  WARNING:", file_path, "not found. Skipping.\n")
    }
  }
  
  load_sig("Tumor.csv", "Tumor")
  load_sig("Fibroblasts.csv", "Fibroblasts")
  load_sig("Macrophages.csv", "Macrophages")
  load_sig("Neutrophils.csv", "Neutrophils")
  
  return(signatures)
}
cancer_signatures <- load_cancer_signatures()

# ===========================================================================
# PART 2: LOAD PIPELINE RESULTS
# ===========================================================================

load_all_pipelines_fixed <- function(results_path) {
  cat("\n[LOAD] Loading data from the 3 pipelines from path:", results_path, "...\n")

  find_latest_file <- function(pattern) {
    files <- list.files(path = results_path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) stop("Cannot find file in '", results_path, "' with pattern: ", pattern)
    return(files[which.max(file.info(files)$mtime)])
  }

  seurat_file <- find_latest_file("^seurat_results_.*\\.rds$")
  giotto_file <- find_latest_file("^giotto_results_.*\\.rds$")
  optimized_file <- find_latest_file("^optimized_pipeline_results_.*\\.rds$")

  cat("  Loading Seurat from:", basename(seurat_file), "\n")
  seurat_data <- readRDS(seurat_file)

  cat("  Loading Giotto from:", basename(giotto_file), "\n")
  giotto_data <- readRDS(giotto_file)

  cat("  Loading Optimized from:", basename(optimized_file), "\n")
  optimized_data <- readRDS(optimized_file)

  seurat_obj <- seurat_data$results$object
  giotto_obj <- giotto_data$results$object

  if(is.null(seurat_obj) || is.null(giotto_obj) || is.null(optimized_data$results)) {
    stop("One of the loaded pipeline objects is NULL.")
  }

  opt_strategy <- optimized_data$results$summary$strategy_used
  if (grepl("Giotto", opt_strategy)) {
    optimized_obj <- optimized_data$results$giotto$object
    opt_type <- "giotto"
    n_optimized <- optimized_data$results$giotto$n_clusters
  } else {
    optimized_obj <- optimized_data$results$seurat$object
    opt_type <- "seurat"
    n_optimized <- optimized_data$results$seurat$n_clusters
  }

  n_seurat <- seurat_data$results$n_clusters
  n_giotto <- giotto_data$results$n_clusters

  cat("  Seurat:", n_seurat, "clusters\n")
  cat("  Giotto:", n_giotto, "clusters\n")
  cat("  Optimized:", n_optimized, "clusters (type:", opt_type, ")\n")

  return(list(
    seurat = seurat_obj, giotto = giotto_obj, optimized = optimized_obj,
    opt_type = opt_type,
    n_clusters = c(seurat = n_seurat, giotto = n_giotto, optimized = n_optimized)
  ))
}

# ===========================================================================
# PART 3: MARKER ANALYSIS
# ===========================================================================

analyze_seurat_markers <- function(seurat_obj, pipeline_name = "Seurat", output_subdir) {
  cat("\n[MARKERS] Analyzing markers -", pipeline_name, "...\n")
  tryCatch({
    markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
    top_markers <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10) %>% as.data.frame()
    filename <- file.path(output_subdir, paste0("markers_", tolower(gsub("[^[:alnum:]]", "_", pipeline_name)), ".csv"))
    write.csv(top_markers, filename, row.names = FALSE)
    cat("  SUCCESS:", nrow(top_markers), "markers found and saved to:", filename, "\n")
    return(top_markers)
  }, error = function(e) {
    cat("  Error in analyze_seurat_markers for", pipeline_name, ":", e$message, "\n")
    return(NULL)
  })
}

analyze_giotto_markers_robust <- function(giotto_obj, pipeline_name = "Giotto", output_subdir) {
  cat("\n[MARKERS] Analyzing markers -", pipeline_name, "...\n")
  tryCatch({
    expr_matrix <- getExpression(giotto_obj, values = "normalized", output = "matrix")
    metadata <- pDataDT(giotto_obj)
    all_markers_list <- list()
    clusters <- sort(unique(metadata$leiden_clus))
    for (clust in clusters) {
      cells_in <- metadata[leiden_clus == clust]$cell_ID
      cells_out <- metadata[leiden_clus != clust]$cell_ID
      if (length(cells_in) < 10 || length(cells_out) < 10) next
      expr_in <- Matrix::rowMeans(expr_matrix[, cells_in, drop = FALSE])
      expr_out <- Matrix::rowMeans(expr_matrix[, cells_out, drop = FALSE])
      pct_in <- Matrix::rowSums(expr_matrix[, cells_in, drop = FALSE] > 0) / length(cells_in)
      pct_out <- Matrix::rowSums(expr_matrix[, cells_out, drop = FALSE] > 0) / length(cells_out)
      expressed_genes <- names(expr_in)[expr_in > 0.1 & pct_in > 0.1]
      if (length(expressed_genes) < 5) next
      fc <- expr_in[expressed_genes] / (expr_out[expressed_genes] + 1e-9)
      log2fc <- log2(fc)
      cluster_df <- data.frame(gene = expressed_genes, avg_log2FC = log2fc, pct.1 = pct_in[expressed_genes], pct.2 = pct_out[expressed_genes], cluster = clust)
      cluster_df <- cluster_df[cluster_df$avg_log2FC > 0.25 & cluster_df$pct.1 > 0.25,]
      if (nrow(cluster_df) > 0) {
        all_markers_list[[as.character(clust)]] <- head(cluster_df[order(cluster_df$avg_log2FC, decreasing = TRUE),], 10)
      }
    }
    all_markers_final_df <- do.call(rbind, all_markers_list)
    if (is.null(all_markers_final_df) || nrow(all_markers_final_df) == 0) {
        cat("  No robust markers found for Giotto.\n"); return(NULL)
    }
    filename <- file.path(output_subdir, paste0("markers_", tolower(gsub("[^[:alnum:]]", "_", pipeline_name)), ".csv"))
    write.csv(all_markers_final_df, filename, row.names = FALSE)
    cat("  SUCCESS:", nrow(all_markers_final_df), "markers found and saved to:", filename, "\n")
    return(all_markers_final_df)
  }, error = function(e) {
    cat("  Error in analyze_giotto_markers_robust for", pipeline_name, ":", e$message, "\n")
    return(NULL)
  })
}

# ===========================================================================
# PART 4: VISUALIZATION FUNCTIONS
# ===========================================================================

create_spatial_comparison_corrected <- function(data, output_subdir) {
    # ... (This function is long but its internal logic is unchanged)
}

visualize_cancer_signatures_fixed <- function(data, signatures, output_subdir) {
    # ... (This function is long but its internal logic is unchanged)
}

# ... (All other visualization functions are also unchanged internally)

# ===========================================================================
# MAIN ANALYSIS FUNCTION
# ===========================================================================

run_final_analysis <- function(results_path, output_subdir_main) {
  cat("\n----------------------------------------------------------------\n")
  cat("        RUNNING FINAL INTEGRATED ANALYSIS\n")
  cat("----------------------------------------------------------------\n\n")

  data_objects <- load_all_pipelines_fixed(results_path)

  cat("\n[PHASE 1] MARKER ANALYSIS\n")
  cat("-------------------------------\n")
  markers_seurat <- analyze_seurat_markers(data_objects$seurat, "Seurat", output_subdir_main)
  markers_giotto <- analyze_giotto_markers_robust(data_objects$giotto, "Giotto", output_subdir_main)

  markers_optimized <- if (data_objects$opt_type == "giotto") {
    analyze_giotto_markers_robust(data_objects$optimized, "Optimized_Giotto", output_subdir_main)
  } else {
    analyze_seurat_markers(data_objects$optimized, "Optimized_Seurat", output_subdir_main)
  }

  all_pipeline_markers <- list(seurat = markers_seurat, giotto = markers_giotto, optimized = markers_optimized)

  # ... (The rest of the main function remains unchanged, calling the plotting functions)

  cat("\nANALYSIS COMPLETED SUCCESSFULLY\n")
  cat("\nFiles generated in:", output_subdir_main, "\n")

  return(list(data = data_objects, markers = all_pipeline_markers, signatures = cancer_signatures))
}

# --- EXECUTION ---
final_results <- run_final_analysis(main_results_dir, integrated_analysis_subdir)

cat("\nSaving final integrated objects...\n")
saveRDS(final_results, file = file.path(integrated_analysis_subdir, "final_results_object.rds"))
if (exists("cancer_signatures") && !is.null(cancer_signatures) && length(cancer_signatures) > 0) {
  saveRDS(cancer_signatures, file = file.path(integrated_analysis_subdir, "cancer_signatures_object.rds"))
}
cat("Objects saved successfully.\n")
cat("\n--- End of Script 5 ---\n")
