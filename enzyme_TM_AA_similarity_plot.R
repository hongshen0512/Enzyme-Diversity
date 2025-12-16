# ================================================================================
# Amino Acid Similarity and TM Score Correlation Analysis
# Description: This script analyzes and visualizes the correlation between 
#              amino acid similarity and TM (Template Modeling) scores
#              for protein structural comparisons
# ================================================================================

# 1. LOAD REQUIRED PACKAGES ----------------------------------------------------
# Install packages if not already installed
required_packages <- c("ggplot2", "patchwork", "dplyr", "tidyr", "ggpubr", 
                       "gridExtra", "cowplot", "viridis")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 2. SETUP AND CONFIGURATION ---------------------------------------------------
# Create output directory
output_dir <- "aa_tm_correlation_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set file paths (uncomment the file you want to analyze)
data_dir <- "data"
input_files <- list(
  "All" = file.path(data_dir, "AA_Similarity_all.txt"),
  "Surface" = file.path(data_dir, "AA_Similarity_surface.txt"),
  "Core" = file.path(data_dir, "AA_Similarity_core.txt")
)

# Choose which file to analyze (default: "All")
selected_region <- "All"
input_file <- input_files[[selected_region]]

# Create log file
log_file <- file.path(output_dir, "correlation_analysis_log.txt")
sink(log_file)
cat("AA Similarity and TM Score Correlation Analysis Log\n")
cat("Analysis started:", date(), "\n\n")
cat("Selected region:", selected_region, "\n")
cat("Input file:", input_file, "\n\n")

# 3. DATA LOADING AND VALIDATION ------------------------------------------------
cat("Loading data...\n")

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Read data
tryCatch({
  df <- read.table(
    input_file, 
    header = TRUE, 
    sep = "\t", 
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cat("Data loaded successfully. Dimensions:", dim(df), "\n")
  cat("Columns:", colnames(df), "\n")
}, error = function(e) {
  stop(paste("Error loading data:", e$message))
})

# Check required columns
required_cols <- c("AA_Similarity", "TM_Score", "source")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# 4. DATA CLEANING AND PREPROCESSING -------------------------------------------
cat("\nCleaning and preprocessing data...\n")

# Convert source to factor with proper levels
df$source <- factor(df$source, levels = c("Non_Rumen", "Rumen"))

# Initial data summary
initial_rows <- nrow(df)
cat("Initial number of observations:", initial_rows, "\n")

# Remove rows with missing values
df_clean <- df[complete.cases(df$AA_Similarity, df$TM_Score), ]
removed_na <- initial_rows - nrow(df_clean)
cat("Rows removed due to missing values:", removed_na, "\n")

# Validate value ranges
df_clean <- df_clean %>%
  filter(
    AA_Similarity >= 0 & AA_Similarity <= 1,
    TM_Score >= 0 & TM_Score <= 1
  )

removed_range <- nrow(df_clean) - (initial_rows - removed_na)
cat("Rows removed due to values outside [0,1] range:", abs(removed_range), "\n")
cat("Final number of observations:", nrow(df_clean), "\n")

# Update df for analysis
df <- df_clean

# 5. DESCRIPTIVE STATISTICS ----------------------------------------------------
cat("\nCalculating descriptive statistics...\n")

# Summary statistics by source
summary_stats <- df %>%
  group_by(source) %>%
  summarise(
    n = n(),
    AA_Similarity_mean = mean(AA_Similarity, na.rm = TRUE),
    AA_Similarity_sd = sd(AA_Similarity, na.rm = TRUE),
    AA_Similarity_median = median(AA_Similarity, na.rm = TRUE),
    TM_Score_mean = mean(TM_Score, na.rm = TRUE),
    TM_Score_sd = sd(TM_Score, na.rm = TRUE),
    TM_Score_median = median(TM_Score, na.rm = TRUE),
    .groups = "drop"
  )

# Print summary statistics
cat("\nSummary statistics by source:\n")
print(summary_stats)

# Save statistics
write.csv(summary_stats, 
          file.path(output_dir, "descriptive_statistics.csv"),
          row.names = FALSE)

# 6. CORRELATION ANALYSIS ------------------------------------------------------
cat("\nPerforming correlation analysis...\n")

# Spearman correlation (non-parametric, robust to outliers)
cor_test <- cor.test(df$AA_Similarity, df$TM_Score, method = "spearman")
rho <- round(cor_test$estimate, 3)
pval <- cor_test$p.value

# Linear regression for visualization
lm_model <- lm(TM_Score ~ AA_Similarity, data = df)
r_squared <- round(summary(lm_model)$r.squared, 3)
lm_pval <- summary(lm_model)$coefficients[2, 4]

# Create prediction data for regression line
pred_df <- data.frame(AA_Similarity = seq(0, 1, length.out = 100))
pred_df$TM_Score <- predict(lm_model, newdata = pred_df)
# Ensure predictions stay within [0,1] range
pred_df$TM_Score <- pmax(0, pmin(1, pred_df$TM_Score))

# Print correlation results
cat("\nCorrelation Analysis Results:\n")
cat("Spearman's rho:", rho, "\n")
cat("p-value:", format(pval, scientific = TRUE, digits = 3), "\n")
cat("Linear regression R-squared:", r_squared, "\n")
cat("Linear regression p-value:", format(lm_pval, scientific = TRUE, digits = 3), "\n")

# Save correlation results
corr_results <- data.frame(
  Region = selected_region,
  Spearman_rho = rho,
  Spearman_p_value = pval,
  Linear_R_squared = r_squared,
  Linear_p_value = lm_pval,
  N_observations = nrow(df),
  Analysis_date = date()
)

write.csv(corr_results, 
          file.path(output_dir, "correlation_results.csv"),
          row.names = FALSE)

# 7. CREATE MAIN SCATTER PLOT --------------------------------------------------
cat("\nCreating main scatter plot...\n")

# Define color palette
source_colors <- c(
  "Non_Rumen" = "#E4DDC0",  # Light beige
  "Rumen" = "#AEC5DA"       # Light blue
)

source_border_colors <- c(
  "Non_Rumen" = "#C4BAA0",  # Darker beige
  "Rumen" = "#8FA8C9"       # Darker blue
)

# Create main scatter plot
main_plot <- ggplot(df, aes(x = AA_Similarity, y = TM_Score)) +
  # Add points with fill and border
  geom_point(
    aes(fill = source, color = source),
    shape = 21, 
    size = 2.5, 
    stroke = 0.5,
    alpha = 0.8
  ) + 
  # Add regression line with confidence interval
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#2D3E50",
    fill = "#95A5A6",
    alpha = 0.2,
    linewidth = 1
  ) +
  # Add custom regression line (for exact control)
  geom_line(
    data = pred_df,
    aes(x = AA_Similarity, y = TM_Score),
    color = "#E74C3C",
    linewidth = 1,
    linetype = "solid"
  ) +
  # Color scales
  scale_fill_manual(values = source_colors, name = "Source") +
  scale_color_manual(values = source_border_colors, name = "Source") +
  # Annotation for correlation statistics
  annotate(
    "text",
    x = 0.02, 
    y = 0.98,
    hjust = 0,
    vjust = 1,
    label = sprintf(
      "Spearman's ¦Ñ = %.3f\np = %s\nn = %d",
      rho,
      ifelse(pval < 0.001, "< 0.001", 
             format(pval, scientific = TRUE, digits = 2)),
      nrow(df)
    ),
    size = 4,
    color = "#2C3E50",
    fontface = "bold"
  ) +
  # Axis limits and breaks
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  # Labels
  labs(
    x = "Amino Acid Similarity",
    y = "TM Score",
    title = paste("AA Similarity vs TM Score Correlation -", selected_region),
    subtitle = "Protein Structure Comparison Analysis"
  ) +
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    plot.title = element_text(
      hjust = 0.5, 
      size = 16, 
      face = "bold",
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, 
      size = 12,
      margin = margin(b = 20)
    ),
    legend.position = c(0.85, 0.15),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key = element_rect(fill = "white"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # Square aspect ratio (1:1)
  coord_fixed(ratio = 1)

# 8. CREATE MARGINAL DENSITY PLOTS ---------------------------------------------
cat("Creating marginal density plots...\n")

# Create top density plot (AA Similarity distribution)
top_density <- ggplot(df, aes(x = AA_Similarity, fill = source)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(
    xintercept = median(df$AA_Similarity),
    linetype = "dashed",
    color = "#E74C3C",
    linewidth = 1
  ) +
  scale_fill_manual(values = source_colors) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Density") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    axis.title.y = element_text(
      size = 9,
      angle = 90,
      vjust = 2,
      hjust = 0.5
    )
  )

# Create right density plot (TM Score distribution)
right_density <- ggplot(df, aes(x = TM_Score, fill = source)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(
    xintercept = median(df$TM_Score),
    linetype = "dashed",
    color = "#E74C3C",
    linewidth = 1
  ) +
  scale_fill_manual(values = source_colors) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Density") +
  coord_flip() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_text(
      size = 9,
      vjust = -1,
      hjust = 0.5
    )
  )

# 9. CREATE EMPTY PLOT FOR LEGEND (IF NEEDED) ----------------------------------
empty_plot <- ggplot() + 
  theme_void()

# 10. COMBINE PLOTS USING PATCHWORK --------------------------------------------
cat("Combining plots...\n")

# Define layout
final_plot <- (
  (top_density + empty_plot + plot_layout(widths = c(1, 0.05))) / 
  (main_plot + right_density + plot_layout(widths = c(1, 0.2))) / 
  plot_layout(heights = c(0.2, 1))
) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '(',
    tag_suffix = ')'
  ) &
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )

# 11. SAVE PLOTS ---------------------------------------------------------------
cat("\nSaving plots...\n")

# Function to save plots in multiple formats
save_plot <- function(plot_obj, filename_base, width = 10, height = 10) {
  # PDF for publications
  ggsave(
    file.path(output_dir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    device = "pdf",
    dpi = 300
  )
  
  # PNG for presentations
  ggsave(
    file.path(output_dir, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    device = "png",
    dpi = 300,
    bg = "white"
  )
  
  # EPS for publications (if Cairo is available)
  if (require("Cairo")) {
    ggsave(
      file.path(output_dir, paste0(filename_base, ".eps")),
      plot = plot_obj,
      width = width,
      height = height,
      device = cairo_ps,
      dpi = 300
    )
  }
  
  cat("Saved:", filename_base, "(.pdf, .png, .eps)\n")
}

# Save combined plot
save_plot(final_plot, paste0("aa_tm_correlation_", tolower(selected_region)), 10, 10)

# Save main plot separately
save_plot(main_plot, paste0("main_scatter_", tolower(selected_region)), 8, 8)

# Save density plots separately
save_plot(top_density, paste0("aa_similarity_density_", tolower(selected_region)), 8, 2)
save_plot(right_density, paste0("tm_score_density_", tolower(selected_region)), 2, 8)

# 12. ADDITIONAL ANALYSES (OPTIONAL) -------------------------------------------
cat("\nPerforming additional analyses...\n")

# Correlation by source
corr_by_source <- df %>%
  group_by(source) %>%
  summarise(
    Spearman_rho = cor(AA_Similarity, TM_Score, method = "spearman"),
    Pearson_r = cor(AA_Similarity, TM_Score, method = "pearson"),
    n = n(),
    .groups = "drop"
  )

cat("\nCorrelation by source:\n")
print(corr_by_source)

# Save source-specific correlations
write.csv(corr_by_source, 
          file.path(output_dir, "correlation_by_source.csv"),
          row.names = FALSE)

# Create faceted plot by source
faceted_plot <- ggplot(df, aes(x = AA_Similarity, y = TM_Score)) +
  geom_point(alpha = 0.6, size = 2, color = "#3498DB") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#F1948A") +
  facet_wrap(~ source, ncol = 2) +
  labs(
    x = "Amino Acid Similarity",
    y = "TM Score",
    title = paste("AA Similarity vs TM Score by Source -", selected_region)
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(color = "white", face = "bold")
  )

save_plot(faceted_plot, paste0("faceted_by_source_", tolower(selected_region)), 10, 5)

# 13. SAVE SESSION INFO --------------------------------------------------------
cat("\nSaving session information...\n")
sink(file.path(output_dir, "session_info.txt"))
cat("Correlation Analysis completed on:", date(), "\n\n")
cat("Input file:", input_file, "\n")
cat("Selected region:", selected_region, "\n")
cat("Number of observations:", nrow(df), "\n")
cat("Spearman correlation (rho):", rho, "\n")
cat("Spearman p-value:", pval, "\n\n")
print(sessionInfo())
sink()

# 14. COMPLETION MESSAGE -------------------------------------------------------
cat("\n" + strrep("=", 60) + "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 60) + "\n")
cat("Results saved in:", output_dir, "\n")
cat("\nFiles generated:\n")
cat("  1. aa_tm_correlation_[region].pdf/.png/.eps - Combined plot\n")
cat("  2. main_scatter_[region].pdf/.png/.eps - Main scatter plot\n")
cat("  3. *_density_[region].pdf/.png/.eps - Marginal density plots\n")
cat("  4. faceted_by_source_[region].pdf/.png/.eps - Faceted plot\n")
cat("  5. descriptive_statistics.csv - Summary statistics\n")
cat("  6. correlation_results.csv - Correlation statistics\n")
cat("  7. correlation_by_source.csv - Source-specific correlations\n")
cat("  8. session_info.txt - Reproducibility information\n")
cat("  9. correlation_analysis_log.txt - Complete analysis log\n")
cat(strrep("=", 60) + "\n")

# Close log file
sink()

# ==============================================================================
# USAGE INSTRUCTIONS:
# 1. Prepare input data file in 'data' directory with columns:
#    - AA_Similarity: Amino acid similarity scores (0-1)
#    - TM_Score: Template Modeling scores (0-1)
#    - source: Source category (e.g., "Rumen", "Non_Rumen")
#
# 2. Set the analysis region (line 28):
#    selected_region <- "All"  # Options: "All", "Surface", "Core"
#
# 3. Run the script:
#    source("aa_tm_correlation_analysis.R")
#
# 4. Results will be saved in 'aa_tm_correlation_results' directory
#
# CUSTOMIZATION OPTIONS:
# - Modify source_colors for different color schemes
# - Adjust plot dimensions in save_plot() function calls
# - Change correlation method (Spearman vs Pearson)
# - Modify statistical significance thresholds
#
# INTERPRETATION:
# - AA Similarity: Measures sequence similarity (0 = no similarity, 1 = identical)
# - TM Score: Measures structural similarity (0 = no similarity, 1 = identical)
# - Positive correlation indicates that sequence similarity predicts structural similarity
# ==============================================================================