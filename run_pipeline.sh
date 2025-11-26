#!/bin/bash

# =============================================================================
# MASTER BASH SCRIPT: SEURIGIOTTO BENCHMARK
# Objective: Run all R scripts sequentially in clean environments.
# =============================================================================

# Function to run a script and check for errors
run_script() {
    script_name=$1
    echo "----------------------------------------------------------------"
    echo ">>> STARTING: $script_name"
    echo "----------------------------------------------------------------"
    
    # Execute the script
    Rscript "$script_name"
    
    # Check the exit status ($?)
    if [ $? -eq 0 ]; then
        echo ">>> SUCCESS: $script_name completed."
        echo ""
    else
        echo "!!! ERROR: $script_name failed. Stopping pipeline."
        exit 1
    fi
}

# === SEQUENCE ===

# 1. Core Benchmarking
run_script "scripts/1_seurat_pipeline.R"
run_script "scripts/2_giotto_pipeline.R"
run_script "scripts/3_compare_pipelines.R"

# 2. Advanced & Integrated Analysis
run_script "scripts/4_optimized_pipeline.R"
run_script "scripts/5_integrated_analysis.R"
run_script "scripts/6_native_visualizations.R"
run_script "scripts/7_final_report_generation.R"
run_script "scripts/8_deconvolution_analysis.R"
run_script "scripts/9_balanced_deconvolution.R"

echo "================================================================"
echo "       🚀 ALL PIPELINES COMPLETED SUCCESSFULLY 🚀"
echo "================================================================"
