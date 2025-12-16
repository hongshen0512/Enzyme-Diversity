# ================================================================
# LEfSe-style Differential Abundance Analysis for Microbiome Data
# Description: This script performs LEfSe-style differential analysis
#              using the MicrobiotaProcess package
# ================================================================

# 1. Load required packages --------------------------------------------------
# Install packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages <- c("phyloseq", "MicrobiotaProcess", "ggplot2", "coin")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("phyloseq", "MicrobiotaProcess")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# 2. Set working directory and create output folder -------------------------
# Uncomment and modify the following line to set your working directory
# setwd("/path/to/your/project")

# Create output directory for results
output_dir <- "lefse_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 3. Define file paths -------------------------------------------------------
# Modify these paths according to your data structure
meta_file <- "data/META.txt"          # Metadata file
otu_file <- "data/ko_frequency.txt"   # OTU/Feature table
output_file <- file.path(output_dir, "species_lefse_results.csv")

# 4. Data loading and validation ---------------------------------------------
cat("Loading metadata from:", meta_file, "\n")
if (!file.exists(meta_file)) {
  stop(paste("Metadata file not found:", meta_file))
}
sampledatamre <- read.table(meta_file, header = TRUE, row.names = 1)
samplemre <- sample_data(sampledatamre)

cat("Loading OTU table from:", otu_file, "\n")
if (!file.exists(otu_file)) {
  stop(paste("OTU table file not found:", otu_file))
}
otutablemre <- read.table(otu_file, header = TRUE, row.names = 1)
otumatmre <- data.matrix(otutablemre)
OTUmre <- otu_table(otumatmre, taxa_are_rows = TRUE)

# 5. Create phyloseq object --------------------------------------------------
cat("Creating phyloseq object...\n")
physeqmre <- phyloseq(OTUmre, samplemre)

# Print summary of the data
cat("\n=== Data Summary ===\n")
cat("Number of samples:", nsamples(physeqmre), "\n")
cat("Number of features:", ntaxa(physeqmre), "\n")
cat("Metadata variables:", colnames(sample_data(physeqmre)), "\n")
cat("Grouping variable: Group\n")
cat("Groups:", unique(sample_data(physeqmre)$Group), "\n")

# 6. Perform LEfSe-style differential analysis --------------------------------
cat("\nPerforming differential analysis...\n")

# Set seed for reproducibility
set.seed(1024)

# Run differential analysis
tryCatch({
  deres <- diff_analysis(
    obj = physeqmre,
    classgroup = "Group",        # Grouping variable in metadata
    firstcomfun = "kruskal_test", # First-level test (group differences)
    padjust = "fdr",             # Multiple testing correction
    filtermod = "pvalue",        # Filter mode
    firstalpha = 0.5,            # First-level significance threshold
    strictmod = TRUE,            # Strict mode
    secondcomfun = "wilcox_test", # Second-level test (pairwise)
    secondalpha = 0.5,           # Second-level significance threshold
    mlfun = "lda",               # Effect size measure (LDA)
    ldascore = 1                 # LDA score cutoff
  )
  
  cat("Differential analysis completed successfully!\n")
  
  # 7. Save results ------------------------------------------------------------
  cat("\nSaving results to:", output_file, "\n")
  write.csv(as.data.frame(deres), file = output_file)
  
  # 8. Generate summary statistics ---------------------------------------------
  cat("\n=== Analysis Summary ===\n")
  cat("Number of significant features:", sum(deres$pvalue < 0.05, na.rm = TRUE), "\n")
  cat("Number of features with LDA score > 2:", sum(deres$LDA > 2, na.rm = TRUE), "\n")
  
  # 9. Optional: Create visualization -----------------------------------------
  cat("\nGenerating visualization...\n")
  
  # Check if plot method is available for deres object
  if (exists("deres") && !is.null(deres)) {
    # Try to generate LDA plot
    try({
      p <- plot(deres)  # This might vary depending on the object structure
      ggsave(
        filename = file.path(output_dir, "lefse_results_plot.pdf"),
        plot = p,
        width = 10,
        height = 8
      )
      cat("Plot saved to:", file.path(output_dir, "lefse_results_plot.pdf"), "\n")
    }, silent = TRUE)
  }
  
  # 10. Save session info for reproducibility ---------------------------------
  sink(file.path(output_dir, "session_info.txt"))
  cat("Analysis completed on:", date(), "\n\n")
  print(sessionInfo())
  sink()
  
  cat("\n=== Analysis Complete ===\n")
  cat("Results saved in:", output_dir, "\n")
  
}, error = function(e) {
  cat("Error in differential analysis:", e$message, "\n")
  stop("Analysis failed. Please check your input data and parameters.")
})

# ============================================================================
# USAGE INSTRUCTIONS:
# 1. Place your metadata in 'data/META.txt' with 'Group' column for comparison
# 2. Place your feature table in 'data/ko_frequency.txt' (features as rows)
# 3. Run this script in R or RStudio
# 4. Results will be saved in the 'lefse_results' folder
#
# PARAMETER EXPLANATION:
# - classgroup: Column name in metadata containing group labels
# - firstcomfun: First-level statistical test (kruskal_test for >2 groups)
# - secondcomfun: Second-level pairwise test
# - padjust: Multiple testing correction method (fdr = false discovery rate)
# - ldascore: Minimum LDA score cutoff for reporting
# ============================================================================