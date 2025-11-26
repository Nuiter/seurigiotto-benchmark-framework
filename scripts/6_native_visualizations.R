#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT 6: NATIVE PIPELINE VISUALIZATION COMPARISON
# Author: Ángel Pérez
# Date: 2 June 2025
# Objective: Compare how each pipeline visualizes its own results using
#            its native plotting functions (e.g., Seurat's SpatialDimPlot vs.
#            Giotto's spatPlot).
# =============================================================================

cat("\n=== SCRIPT 6: NATIVE VISUALIZATION OF EACH PIPELINE ===\n")
cat("Objective: Compare how each tool visualizes its own results\n\n")

# --- 0. INITIAL SETUP ---
# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Giotto)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
})

# --- CONFIGURATION & PATHS ---
# This script depends on the output of Script 5.
main_results_dir <- "results"
# Source directory for the integrated results object
integrated_analysis_subdir_ref <- file.path(main_results_dir, "05_Integrated_Pipeline_Analysis")
# Specific output directory for this script's plots
native_viz_output_dir <- file.path(main_results_dir, "06_Native_Pipeline_Visualizations")
dir.create(native_viz_output_dir, showWarnings = FALSE, recursive = TRUE)
# --- END CONFIGURATION ---

# --- 1. LOAD INTEGRATED RESULTS OBJECT ---
cat("1. Loading the 'final_results_object.rds'...\n")

final_results_file <- file.path(integrated_analysis_subdir_ref, "final_results_object.rds")

if (file.exists(final_results_file)) {
  final_results <- readRDS(final_results_file)
  cat("   SUCCESS: 'final_results' object loaded.\n")
} else {
  stop(paste(
    "ERROR: 'final_results_object.rds' not found in:", integrated_analysis_subdir_ref,
    "\nPlease ensure you have successfully run the '5_integrated_analysis.R' script first."
  ))
}

# Verify that the required data structure exists
if (!exists("final_results") || !("data" %in% names(final_results))) {
  stop("ERROR: The 'final_results' object is missing or does not contain the 'data' element.")
}


# =============================================================================
# FUNCTION 1: NATIVE SEURAT VISUALIZATION
# =============================================================================

visualize_seurat_native <- function(seurat_obj, save_name = "seurat_native", output_dir, custom_title = NULL) {
  cat("\n[SEURAT] Creating native visualization with SpatialDimPlot...\n")
  
  # Define Title
  plot_title <- if(!is.null(custom_title)) custom_title else "Seurat - Native Visualization"

  p_seurat <- SpatialDimPlot(
    seurat_obj,
    group.by = "seurat_clusters", # Assumes clusters are in Idents() or this column
    pt.size.factor = 1.6,
    image.alpha = 0, # Remove background image for clarity
    alpha = 1,       # Full opacity for points
    stroke = 0.1,
    label = TRUE,
    label.size = 3,
    repel = TRUE
  ) +
    ggtitle(plot_title,
            subtitle = paste0(length(unique(Idents(seurat_obj))), " clusters")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )

  ggsave(file.path(output_dir, paste0(save_name, ".png")),
         p_seurat, width = 8, height = 8, dpi = 300)

  cat("   Seurat visualization saved to:", output_dir, "\n")
  return(p_seurat)
}

# =============================================================================
# FUNCTION 2: NATIVE GIOTTO VISUALIZATION
# =============================================================================

visualize_giotto_native <- function(giotto_obj, save_name = "giotto_native", output_dir, custom_title = NULL) {
  cat("\n[GIOTTO] Creating native visualization with spatPlot...\n")

  cluster_column_name <- "leiden_clus" # Standard Giotto cluster column
  if (!cluster_column_name %in% colnames(pDataDT(giotto_obj))) {
    stop("ERROR: Cluster column '", cluster_column_name, "' not found in Giotto object.")
  }

  n_clusters <- length(unique(pDataDT(giotto_obj)[[cluster_column_name]]))
  cat("   Clusters detected:", n_clusters, "\n")
  
  # Define Title
  plot_title <- if(!is.null(custom_title)) custom_title else "Giotto - Native Visualization"

  # Create a consistent color palette
  if (n_clusters > 0) {
    if (n_clusters <= 12) {
      colors_giotto <- brewer.pal(max(3, n_clusters), "Set3")
      if (n_clusters < 3) colors_giotto <- colors_giotto[1:n_clusters]
    } else {
      colors_giotto <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)
    }
  } else {
    colors_giotto <- "blue" # Fallback
  }

  p_giotto <- spatPlot(
    gobject = giotto_obj,
    cell_color = cluster_column_name,
    point_size = 1.2,
    cell_color_code = colors_giotto,
    show_image = FALSE,
    show_plot = FALSE, # Required to modify the ggplot object
    return_plot = TRUE,
    axis_text = 12,
    axis_title = 12,
    legend_text = 10
  )

  if (!is.null(p_giotto)) {
    p_giotto <- p_giotto +
      ggtitle(plot_title,
              subtitle = paste0(n_clusters, " clusters")) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right"
      ) +
      labs(color = "Cluster")

    ggsave(file.path(output_dir, paste0(save_name, ".png")),
           p_giotto, width = 8, height = 8, dpi = 300)

    cat("   Giotto visualization saved to:", output_dir, "\n")
  }
  return(p_giotto)
}

# =============================================================================
# FUNCTION 3: SIDE-BY-SIDE NATIVE COMPARISON PLOT
# =============================================================================

create_native_comparison_plot <- function(data, output_dir) {
  cat("\n[COMPARISON] Creating side-by-side native visualization...\n")

  # 1. Process Seurat
  p_seurat <- SpatialDimPlot(
    data$seurat, pt.size.factor = 1.6, stroke = 0.1, label = FALSE,
    image.alpha = 0, alpha = 1
  ) +
    ggtitle("Seurat Pipeline", subtitle = paste0(data$n_clusters["seurat"], " clusters")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )

  # 2. Process Giotto
  p_giotto <- spatPlot(
    gobject = data$giotto, cell_color = "leiden_clus", point_size = 1.2,
    show_image = FALSE, show_plot = FALSE, return_plot = TRUE
  ) +
    ggtitle("Giotto Pipeline", subtitle = paste0(data$n_clusters["giotto"], " clusters")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )

  # 3. Process Optimized Pipeline (conditionally)
  if (data$opt_type == "giotto") {
    p_optimized <- spatPlot(
      gobject = data$optimized, cell_color = "leiden_clus", point_size = 1.2,
      show_image = FALSE, show_plot = FALSE, return_plot = TRUE
    )
  } else { # seurat
    p_optimized <- SpatialDimPlot(
      data$optimized, pt.size.factor = 1.6, stroke = 0.1, label = FALSE,
      image.alpha = 0, alpha = 1
    )
  }
  
  # Apply common theme and correct title
  p_optimized <- p_optimized +
    ggtitle("Optimized Pipeline", subtitle = paste0(data$n_clusters["optimized"], " clusters")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )

  # Combine all three
  combined_native <- p_seurat + p_giotto + p_optimized +
    plot_layout(ncol = 3) +
    plot_annotation(
      title = "Comparison of Native Visualizations",
      subtitle = "Each pipeline visualized with its own plotting function",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14)
      )
    )

  ggsave(file.path(output_dir, "native_comparison_plot.png"),
         combined_native, width = 18, height = 6, dpi = 300)

  cat("   Native comparison plot saved to:", output_dir, "\n")
  return(combined_native)
}


# =============================================================================
# MAIN FUNCTION: EXECUTE ALL NATIVE VISUALIZATIONS
# =============================================================================

run_native_visualizations <- function(data, output_dir) {
  cat("\n---------------------------------------------------------------\n")
  cat("       GENERATING NATIVE VISUALIZATIONS FOR EACH PIPELINE      \n")
  cat("---------------------------------------------------------------\n\n")

  # 1. Create individual, detailed native plots
  cat("[PHASE 1] Individual native visualizations\n")
  cat("-----------------------------------------------\n")
  
  # Standard Seurat and Giotto (Default titles are fine)
  p_seurat <- visualize_seurat_native(data$seurat, "seurat_native_detailed", output_dir)
  p_giotto <- visualize_giotto_native(data$giotto, "giotto_native_detailed", output_dir)

  # Optimized Pipeline - WE FORCE THE TITLE HERE
  if (data$opt_type == "giotto") {
    visualize_giotto_native(
      data$optimized, 
      "optimized_native_detailed", 
      output_dir, 
      custom_title = "Optimized Pipeline - Native Visualization"
    )
  } else {
    visualize_seurat_native(
      data$optimized, 
      "optimized_native_detailed", 
      output_dir,
      custom_title = "Optimized Pipeline - Native Visualization"
    )
  }

  # 2. Create side-by-side comparison plot
  cat("\n[PHASE 2] Side-by-side comparison\n")
  cat("-----------------------------------------------\n")
  create_native_comparison_plot(data, output_dir)

  # 3. Summary
  cat("\nNATIVE VISUALIZATIONS COMPLETE\n")
  cat("\nFiles generated in:", output_dir, "\n")
  cat("- seurat_native_detailed.png\n")
  cat("- giotto_native_detailed.png\n")
  cat("- optimized_native_detailed.png\n")
  cat("- native_comparison_plot.png\n")

  cat("\nNOTE ON ORIENTATION:\n")
  cat("- Giotto's spatPlot may invert the Y-axis by default.\n")
  cat("- Seurat's SpatialDimPlot maintains the original orientation.\n")
  cat("- The side-by-side comparison shows each tool's default rendering.\n")
}

# =============================================================================
# EXECUTE VISUALIZATION SCRIPT
# =============================================================================

run_native_visualizations(final_results$data, native_viz_output_dir)

cat("\nAll native visualizations have been successfully generated.\n")
cat("Please check the directory:", native_viz_output_dir, "\n")