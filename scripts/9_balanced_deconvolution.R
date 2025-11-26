#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 9: BALANCED DECONVOLUTION & TUMOR MICROENVIRONMENT ANALYSIS
# Author: Ángel Pérez
# Date: 2 June 2025
# Objective: Perform a refined deconvolution using a balanced marker set
#            to prevent dominant cell types (like fibroblasts) from masking
#            subtler signals in the tumor microenvironment (TME).
# =============================================================================

cat("\n=== SCRIPT 9: BALANCED DECONVOLUTION ANALYSIS ===\n")
cat("Using a balanced marker set to prevent signal dominance\n\n")

# --- 0. INITIAL SETUP ---
suppressPackageStartupMessages({
  library(Seurat); library(Giotto); library(ggplot2); library(dplyr);
  library(patchwork); library(viridis); library(tidyr); library(RColorBrewer);
})

# --- CONFIGURATION & PATHS ---
main_results_dir <- "results"
integrated_analysis_subdir_ref <- file.path(main_results_dir, "05_Integrated_Pipeline_Analysis")
balanced_decon_output_dir <- file.path(main_results_dir, "09_Balanced_Deconvolution")
dir.create(balanced_decon_output_dir, showWarnings = FALSE, recursive = TRUE)
# --- END CONFIGURATION ---

# --- 1. LOAD DATA ---
cat("1. Loading 'final_results' object...\n")
final_results_file <- file.path(integrated_analysis_subdir_ref, "final_results_object.rds")
if (!file.exists(final_results_file)) stop("ERROR: 'final_results_object.rds' not found. Run Script 5 first.")
final_results <- readRDS(final_results_file)
cat("   SUCCESS: 'final_results' object loaded.\n")

# ===========================================================================
# PART 2: CORE ANALYSIS FUNCTIONS
# ===========================================================================

get_balanced_markers <- function() {
  cat("\n[MARKERS] Using a small, highly specific marker set...\n")
  return(list(
    Tumor = c("EPCAM", "KRT19"),
    Fibroblasts = c("COL1A1", "DCN"),
    Macrophages = c("CD68", "CD163"),
    Neutrophils = c("S100A8", "S100A9"),
    T_cells = c("CD3E", "CD3D"),
    B_cells = c("CD79A", "MS4A1")
  ))
}

analyze_composition_robust <- function(data, markers) {
  cat("\n[ANALYSIS] Analyzing cell composition for all pipelines...\n")
  results <- list()
  
  for (pipeline in c("seurat", "giotto", "optimized")) {
    cat("  -> Processing pipeline:", toupper(pipeline), "\n")
    obj <- data[[pipeline]]
    
    # Get expression data and cluster IDs consistently
    if (pipeline == "seurat" || (pipeline == "optimized" && data$opt_type == "seurat")) {
      clusters <- Idents(obj)
      expr_data <- GetAssayData(obj, layer = "data")
    } else { # Giotto
      meta <- pDataDT(obj)
      clusters <- factor(meta$leiden_clus)
      expr_data <- getExpression(obj, values = "normalized", output = "matrix")
      names(clusters) <- meta$cell_ID
      # Ensure alignment
      common_cells <- intersect(names(clusters), colnames(expr_data))
      clusters <- clusters[common_cells]
      expr_data <- expr_data[, common_cells]
    }

    # Calculate scores per cluster using a robust method (75th percentile)
    cluster_scores <- sapply(levels(clusters), function(cl) {
      cells_in_cluster <- names(clusters)[clusters == cl]
      if(length(cells_in_cluster) == 0) return(rep(0, length(markers))) # Handle empty clusters
      
      sapply(markers, function(m) {
        available <- intersect(m, rownames(expr_data))
        if (length(available) > 0) {
          # Handle single gene case vs multiple genes
          if(length(available) == 1) {
             val <- expr_data[available, cells_in_cluster]
             quantile(val, 0.75)
          } else {
             quantile(colMeans(expr_data[available, cells_in_cluster, drop = FALSE]), 0.75)
          }
        } else { 0 }
      })
    })

    # Convert scores to probabilities using softmax
    # Transpose so rows are clusters, cols are cell types
    cluster_scores_t <- t(cluster_scores)
    cluster_probs <- apply(scale(cluster_scores_t), 1, function(x) {
        # Handle potential NAs from scaling constant vectors
        x[is.na(x)] <- 0 
        exp(x) / sum(exp(x))
    })
    cluster_probs <- t(cluster_probs) # Transpose back
    colnames(cluster_probs) <- names(markers)

    # Calculate overall tissue composition, weighted by cluster size
    cluster_sizes <- table(clusters)
    # Ensure order matches
    cluster_sizes <- cluster_sizes[rownames(cluster_probs)]
    
    tissue_comp <- colSums(cluster_probs * (as.numeric(cluster_sizes) / sum(cluster_sizes)))
    
    results[[pipeline]] <- list(
      cluster_probs = cluster_probs,
      tissue_composition = tissue_comp,
      diversity = -sum(tissue_comp[tissue_comp > 0] * log(tissue_comp[tissue_comp > 0]))
    )
  }
  return(results)
}

# ===========================================================================
# PART 3: VISUALIZATION & REPORTING
# ===========================================================================

create_comparison_plots <- function(results_full, results_tme, data, output_dir) {
  cat("\n[PLOTS] Creating comparative visualizations...\n")

  # --- Helper to create data.frame for plotting ---
  create_plot_df <- function(results) {
    do.call(rbind, lapply(names(results), function(p) {
      data.frame(
        Pipeline = p,
        Cell_Type = names(results[[p]]$tissue_composition),
        Proportion = as.numeric(results[[p]]$tissue_composition)
      )
    }))
  }
  
  df_full <- create_plot_df(results_full)
  df_tme <- create_plot_df(results_tme)
  
  # --- Plot 1: Full Composition ---
  p1 <- ggplot(df_full, aes(x = Cell_Type, y = Proportion, fill = Pipeline)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "A. Full Tissue Composition (with Fibroblasts)", x = "", y = "Proportion")

  # --- Plot 2: TME Composition ---
  p2 <- ggplot(df_tme, aes(x = Cell_Type, y = Proportion, fill = Pipeline)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "B. Tumor Microenvironment (Fibroblasts Excluded)", x = "", y = "Proportion")

  # Combine and save
  combined_plot <- (p1 + p2) + plot_layout(guides = "collect") &
    theme_minimal() &
    scale_fill_manual(
      values = c(seurat = "#4CAF50", giotto = "#2196F3", optimized = "#FF5722"),
      labels = paste0(c("Seurat", "Giotto", "Optimized"), " (", data$n_clusters, " cl)")
    ) &
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(output_dir, "01_composition_full_vs_tme.png"),
         combined_plot, width = 12, height = 8, dpi = 300)
  cat("   - Saved composition comparison plot.\n")
}

analyze_immune_tumor_ratios <- function(results_tme, data, output_dir) {
  cat("\n[RATIOS] Analyzing Immune-to-Tumor Ratios...\n")
  
  ratio_data <- do.call(rbind, lapply(names(results_tme), function(p) {
    comp <- results_tme[[p]]$tissue_composition
    
    # SAFE ACCESS to named vector
    tumor_prop <- if("Tumor" %in% names(comp)) comp[["Tumor"]] else 0
    
    # Sum immune cells (check existence)
    immune_types <- c("Macrophages", "Neutrophils", "T_cells", "B_cells")
    existing_immune <- intersect(immune_types, names(comp))
    immune_prop <- if(length(existing_immune) > 0) sum(comp[existing_immune]) else 0
    
    data.frame(
      Pipeline = p,
      Ratio = immune_prop / (tumor_prop + 1e-6), # Avoid division by zero
      N_Clusters = data$n_clusters[p]
    )
  }))
  
  p_ratio <- ggplot(ratio_data, aes(x = factor(N_Clusters), y = Ratio, fill = Pipeline)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(Ratio, 2)), vjust = -0.5) +
    scale_fill_manual(values = c(seurat = "#4CAF50", giotto = "#2196F3", optimized = "#FF5722")) +
    theme_minimal() +
    labs(
      title = "Immune-to-Tumor Ratio in the TME",
      subtitle = "Higher granularity reveals greater estimated immune infiltration",
      x = "Number of Clusters in Pipeline", y = "Ratio (Immune / Tumor Cells)"
    ) + theme(legend.position = "none")

  ggsave(file.path(output_dir, "02_immune_tumor_ratio.png"), p_ratio, width = 8, height = 6, dpi = 300)
  cat("   - Saved immune-to-tumor ratio plot.\n")
}

# ===========================================================================
# MAIN FUNCTION
# ===========================================================================

run_balanced_deconvolution_workflow <- function() {
  cat("\n---------------------------------------------------------------\n")
  cat("         RUNNING BALANCED DECONVOLUTION & TME ANALYSIS         \n")
  cat("---------------------------------------------------------------\n")

  markers <- get_balanced_markers()

  # 1. Full composition analysis
  cat("\n--- ANALYSIS 1: FULL TISSUE COMPOSITION ---\n")
  results_full <- analyze_composition_robust(final_results$data, markers)

  # 2. TME-focused analysis (excluding fibroblasts)
  cat("\n--- ANALYSIS 2: TUMOR MICROENVIRONMENT (TME) ONLY ---\n")
  tme_markers <- markers; tme_markers$Fibroblasts <- NULL
  results_tme <- analyze_composition_robust(final_results$data, tme_markers)

  # 3. Visualizations and final analysis
  create_comparison_plots(results_full, results_tme, final_results$data, balanced_decon_output_dir)
  analyze_immune_tumor_ratios(results_tme, final_results$data, balanced_decon_output_dir)

  cat("\nBALANCED ANALYSIS COMPLETE!\n")
  cat("Results saved in:", balanced_decon_output_dir, "\n")
}

# ===========================================================================
# EXECUTION
# ===========================================================================

run_balanced_deconvolution_workflow()

cat("\n--- END OF SCRIPT 9 ---\n")
