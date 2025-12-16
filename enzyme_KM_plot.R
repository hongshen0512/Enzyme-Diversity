# ================================================================================
# Enzyme Kinetics Parameter Analysis and Visualization
# Description: This script analyzes and visualizes enzyme kinetic parameters (KM, 
#              Vmax, Kcat, Kcat/KM) across different sample clusters
# ================================================================================

# 1. LOAD REQUIRED PACKAGES ----------------------------------------------------
# Install packages if not already installed
required_packages <- c("ggplot2", "dplyr", "ggpubr", "gridExtra", "rstatix")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 2. SETUP AND CONFIGURATION ---------------------------------------------------
# Create output directory
output_dir <- "enzyme_kinetics_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set file paths (modify these according to your data structure)
input_file <- "data/enzyme_data.txt"  # Input data file
output_plot <- file.path(output_dir, "enzyme_kinetics_plot")
log_file <- file.path(output_dir, "analysis_log.txt")

# Initialize log
sink(log_file)
cat("Enzyme Kinetics Analysis Log\n")
cat("Analysis started:", date(), "\n\n")

# 3. LOAD AND PREPARE DATA -----------------------------------------------------
cat("Loading data from:", input_file, "\n")
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Read enzyme data
tryCatch({
  enzyme_data <- read.delim(
    input_file, 
    header = TRUE, 
    stringsAsFactors = FALSE, 
    fileEncoding = "UTF-8",
    check.names = FALSE  # Preserve column names as they are
  )
  cat("Data loaded successfully. Dimensions:", dim(enzyme_data), "\n")
}, error = function(e) {
  stop(paste("Error loading data:", e$message))
})

# Transpose data for easier analysis
enzyme_t <- as.data.frame(t(enzyme_data[,-1]))
colnames(enzyme_t) <- enzyme_data[,1]

# Define sample grouping (modify this based on your experimental design)
# Assuming 4 replicates per cluster
cat("\nDefining sample groups...\n")
num_clusters <- 4
samples_per_cluster <- 4

# Create grouping vector
cluster_labels <- paste0("cluster_", 3:(3 + num_clusters - 1))
enzyme_t$sample <- rep(cluster_labels, each = samples_per_cluster)

cat("Cluster labels:", cluster_labels, "\n")
cat("Samples per cluster:", samples_per_cluster, "\n")
cat("Total samples:", nrow(enzyme_t), "\n")

# 4. DATA CLEANING AND VALIDATION ----------------------------------------------
cat("\nConverting data types and validating...\n")

# List of expected parameters (modify if your data has different columns)
expected_params <- c("KM", "Vmax", "Kcat", "Kcat/KM")

# Check if all expected parameters are present
missing_params <- setdiff(expected_params, colnames(enzyme_t))
if (length(missing_params) > 0) {
  warning(paste("Missing parameters:", paste(missing_params, collapse = ", ")))
}

# Convert numeric columns
for (param in expected_params) {
  if (param %in% colnames(enzyme_t)) {
    enzyme_t[[param]] <- as.numeric(enzyme_t[[param]])
    cat(param, ": Converted to numeric. NAs:", sum(is.na(enzyme_t[[param]])), "\n")
  }
}

# Factorize sample groups with specified order
enzyme_t$sample <- factor(enzyme_t$sample, levels = cluster_labels)

# Remove rows with all NA values for parameters
complete_cases <- complete.cases(enzyme_t[, expected_params[expected_params %in% colnames(enzyme_t)]])
if (sum(!complete_cases) > 0) {
  cat("Removing", sum(!complete_cases), "rows with missing values\n")
  enzyme_t <- enzyme_t[complete_cases, ]
}

# 5. STATISTICAL ANALYSIS FUNCTION ---------------------------------------------
create_kinetics_plot <- function(data, y_var, y_label, title, 
                                 show_stats = TRUE, 
                                 p_adjust_method = "none") {
  """
  Creates a bar plot with error bars and statistical significance annotations
  
  Args:
    data: Data frame containing the data
    y_var: Column name of the variable to plot
    y_label: Label for y-axis
    title: Plot title
    show_stats: Whether to show statistical comparisons
    p_adjust_method: Method for p-value adjustment (none, bonferroni, fdr, etc.)
  
  Returns:
    ggplot object
  """
  
  # Calculate summary statistics
  summary_data <- data %>%
    group_by(sample) %>%
    summarise(
      mean = mean(!!sym(y_var), na.rm = TRUE),
      sd = sd(!!sym(y_var), na.rm = TRUE),
      se = sd / sqrt(n()),
      n = n(),
      .groups = 'drop'
    )
  
  # Set color palette
  cluster_colors <- c(
    "cluster_3" = "#1f77b4", 
    "cluster_4" = "#ff7f0e", 
    "cluster_5" = "#2ca02c", 
    "cluster_6" = "#d62728"
  )
  
  # Adjust colors based on available clusters
  available_clusters <- levels(data$sample)
  cluster_colors <- cluster_colors[names(cluster_colors) %in% available_clusters]
  
  # Create base plot
  p <- ggplot(summary_data, aes(x = sample, y = mean)) +
    geom_col(aes(fill = sample), width = 0.7, alpha = 0.8) +
    geom_errorbar(
      aes(ymin = mean - sd, ymax = mean + sd), 
      width = 0.2, 
      size = 0.7,
      color = "black"
    ) +
    geom_point(
      data = data,
      aes(x = sample, y = !!sym(y_var)),
      position = position_jitter(width = 0.2),
      size = 2,
      alpha = 0.6
    ) +
    scale_fill_manual(values = cluster_colors) +
    labs(
      title = title,
      y = y_label,
      x = "Cluster",
      fill = "Cluster"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
  
  # Add statistical comparisons if requested
  if (show_stats && length(available_clusters) > 1) {
    
    # Perform pairwise t-tests with p-value adjustment
    pairwise_tests <- data %>%
      t_test(
        as.formula(paste(y_var, "~ sample")),
        p.adjust.method = p_adjust_method,
        detailed = TRUE
      )
    
    # Filter significant comparisons
    sig_comparisons <- pairwise_tests %>%
      filter(p < 0.05) %>%
      mutate(
        p.adj.signif = case_when(
          p < 0.001 ~ "***",
          p < 0.01 ~ "**",
          p < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
    
    # Calculate y-axis position for significance bars
    y_max <- max(summary_data$mean + summary_data$sd, na.rm = TRUE)
    y_step <- y_max * 0.1
    
    # Add significance annotations
    if (nrow(sig_comparisons) > 0) {
      for (i in 1:nrow(sig_comparisons)) {
        group1 <- sig_comparisons$group1[i]
        group2 <- sig_comparisons$group2[i]
        p_val <- sig_comparisons$p.adj.signif[i]
        
        x1 <- which(available_clusters == group1)
        x2 <- which(available_clusters == group2)
        y_pos <- y_max + (i * y_step)
        
        p <- p +
          # Significance bracket
          annotate("segment", 
                   x = x1, xend = x2, 
                   y = y_pos, yend = y_pos,
                   size = 0.8) +
          annotate("segment",
                   x = x1, xend = x1,
                   y = y_pos - 0.02 * y_max, yend = y_pos,
                   size = 0.8) +
          annotate("segment",
                   x = x2, xend = x2,
                   y = y_pos - 0.02 * y_max, yend = y_pos,
                   size = 0.8) +
          # Significance stars
          annotate("text",
                   x = (x1 + x2)/2,
                   y = y_pos + 0.02 * y_max,
                   label = p_val,
                   size = 5,
                   fontface = "bold")
      }
      
      # Adjust y-axis limit to accommodate significance markers
      p <- p + expand_limits(y = y_max + (nrow(sig_comparisons) + 1) * y_step)
    } else {
      cat("No significant differences found for", y_var, "\n")
    }
    
    # Save statistical results
    stats_file <- file.path(output_dir, paste0("stats_", y_var, ".csv"))
    write.csv(pairwise_tests, stats_file, row.names = FALSE)
    cat("Statistical results saved to:", stats_file, "\n")
  }
  
  return(p)
}

# 6. CREATE VISUALIZATIONS -----------------------------------------------------
cat("\nCreating visualizations...\n")

# Define parameters to plot
parameters <- list(
  KM = list(var = "KM", label = "KM (¦ÌM)", title = "Michaelis Constant (KM)"),
  Vmax = list(var = "Vmax", label = "Vmax (¦ÌM/min)", title = "Maximum Velocity (Vmax)"),
  Kcat = list(var = "Kcat", label = "Kcat (min?1)", title = "Turnover Number (Kcat)"),
  Kcat_KM = list(var = "Kcat/KM", label = "Kcat/KM (¦ÌM?1 min?1)", 
                 title = "Catalytic Efficiency (Kcat/KM)")
)

# Create individual plots
plots <- list()
for (param_name in names(parameters)) {
  if (parameters[[param_name]]$var %in% colnames(enzyme_t)) {
    cat("Creating plot for:", param_name, "\n")
    
    plots[[param_name]] <- create_kinetics_plot(
      data = enzyme_t,
      y_var = parameters[[param_name]]$var,
      y_label = parameters[[param_name]]$label,
      title = parameters[[param_name]]$title,
      show_stats = TRUE,
      p_adjust_method = "bonferroni"  # Adjust for multiple comparisons
    )
    
    # Save individual plot
    ggsave(
      filename = file.path(output_dir, paste0("enzyme_", param_name, ".pdf")),
      plot = plots[[param_name]],
      width = 6,
      height = 5,
      device = "pdf"
    )
  } else {
    warning(paste("Parameter not found in data:", parameters[[param_name]]$var))
  }
}

# 7. CREATE COMBINED PLOT ------------------------------------------------------
cat("\nCreating combined plot...\n")

# Arrange plots in grid
if (length(plots) > 0) {
  combined_plot <- grid.arrange(
    grobs = plots,
    ncol = min(4, length(plots)),
    top = textGrob(
      "Enzyme Kinetic Parameters by Cluster",
      gp = gpar(fontsize = 16, fontface = "bold")
    )
  )
  
  # Save combined plot in multiple formats
  save_plot <- function(filename, width, height) {
    ggsave(
      filename = file.path(output_dir, filename),
      plot = combined_plot,
      width = width,
      height = height,
      dpi = 300
    )
    cat("Saved:", filename, "\n")
  }
  
  # PDF for publications
  save_plot("enzyme_kinetics_combined.pdf", width = 16, height = 5)
  
  # EPS for publications (requires extra packages)
  if (require("Cairo")) {
    ggsave(
      filename = file.path(output_dir, "enzyme_kinetics_combined.eps"),
      plot = combined_plot,
      width = 16,
      height = 5,
      device = cairo_ps,
      dpi = 300
    )
    cat("Saved: enzyme_kinetics_combined.eps\n")
  }
  
  # PNG for presentations
  save_plot("enzyme_kinetics_combined.png", width = 16, height = 5)
}

# 8. GENERATE SUMMARY STATISTICS -----------------------------------------------
cat("\nGenerating summary statistics...\n")

# Summary table by cluster
summary_stats <- enzyme_t %>%
  group_by(sample) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        se = ~sd(., na.rm = TRUE)/sqrt(n()),
        n = ~sum(!is.na(.))
      ),
      .names = "{.col}_{.fn}"
    )
  )

# Save summary statistics
write.csv(
  summary_stats,
  file.path(output_dir, "enzyme_kinetics_summary.csv"),
  row.names = FALSE
)

# 9. SESSION INFO AND LOGGING --------------------------------------------------
cat("\nSaving session information...\n")

# Save comprehensive session info
sink(file.path(output_dir, "session_info.txt"))
cat("Analysis completed on:", date(), "\n\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n\n")
print(sessionInfo())
sink()

# Close log file
cat("\nAnalysis completed successfully!\n")
cat("Results saved in:", output_dir, "\n")
sink()

# ==============================================================================
# USAGE INSTRUCTIONS:
# 1. Prepare your data file (tab-delimited) with the following structure:
#    - First column: Parameter names (KM, Vmax, Kcat, Kcat/KM)
#    - Subsequent columns: Sample measurements
# 2. Modify the 'input_file' path to point to your data
# 3. Adjust 'cluster_labels' and 'samples_per_cluster' according to your design
# 4. Run the script: source("enzyme_kinetics_analysis.R")
#
# OUTPUT FILES:
# - enzyme_kinetics_results/enzyme_KM.pdf: Individual KM plot
# - enzyme_kinetics_results/enzyme_Vmax.pdf: Individual Vmax plot
# - enzyme_kinetics_results/enzyme_Kcat.pdf: Individual Kcat plot
# - enzyme_kinetics_results/enzyme_Kcat_KM.pdf: Individual Kcat/KM plot
# - enzyme_kinetics_results/enzyme_kinetics_combined.pdf: Combined plot
# - enzyme_kinetics_results/enzyme_kinetics_combined.eps: EPS format
# - enzyme_kinetics_results/enzyme_kinetics_combined.png: PNG format
# - enzyme_kinetics_results/stats_*.csv: Statistical test results
# - enzyme_kinetics_results/enzyme_kinetics_summary.csv: Summary statistics
# - enzyme_kinetics_results/analysis_log.txt: Analysis log
# - enzyme_kinetics_results/session_info.txt: Session information
# ==============================================================================