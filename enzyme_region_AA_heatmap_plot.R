# ================================================================================
# Combined Amino Acid Enrichment Heatmap and Preference Analysis
# Description: This script creates combined visualizations of amino acid 
#              enrichment analysis (heatmap) and regional preference analysis
#              (line plot) for protein structural regions
# ================================================================================

# 1. LOAD REQUIRED PACKAGES ----------------------------------------------------
# Install packages if not already installed
required_packages <- c("tidyverse", "cowplot", "patchwork", "ggtext", "RColorBrewer")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 2. SETUP AND CONFIGURATION ---------------------------------------------------
# Create output directory
output_dir <- "aa_combined_visualization_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set file paths
enrichment_file <- "data/AA_region_enrichment_results.txt"
preference_file <- "data/region_aa_preference.txt"

# Create log file
log_file <- file.path(output_dir, "visualization_log.txt")
sink(log_file)
cat("Amino Acid Combined Visualization Log\n")
cat("Analysis started:", date(), "\n\n")

# 3. DEFINE AMINO ACID ORDER AND PROPERTIES ------------------------------------
# Standard amino acid order (3-letter codes)
amino_acids <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", 
                 "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", 
                 "Thr", "Val", "Trp", "Tyr")

# Amino acid classification for grouping
aa_classification <- list(
  "Aliphatic" = c("Ala", "Gly", "Ile", "Leu", "Val"),
  "Aromatic" = c("Phe", "Trp", "Tyr"),
  "Acidic" = c("Asp", "Glu"),
  "Basic" = c("Arg", "His", "Lys"),
  "Hydroxylic" = c("Ser", "Thr"),
  "Sulfur" = c("Cys", "Met"),
  "Amide" = c("Asn", "Gln"),
  "Cyclic" = c("Pro")
)

# Create a data frame for amino acid grouping
aa_groups <- data.frame(
  AA = character(),
  Group = character(),
  stringsAsFactors = FALSE
)

for (group_name in names(aa_classification)) {
  group_aas <- aa_classification[[group_name]]
  aa_groups <- rbind(aa_groups, 
                     data.frame(AA = group_aas, Group = group_name))
}

# 4. LOAD AND PREPARE DATA -----------------------------------------------------
cat("Loading enrichment data from:", enrichment_file, "\n")

if (!file.exists(enrichment_file)) {
  stop(paste("Enrichment file not found:", enrichment_file))
}

# Read enrichment data
enrich_data <- read.table(enrichment_file, header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)

# Ensure proper factor levels
enrich_data$AA <- factor(enrich_data$AA, levels = amino_acids)
enrich_data$Region <- factor(enrich_data$Region, 
                             levels = c("Core", "Surface", "Key"))

cat("Enrichment data loaded. Dimensions:", dim(enrich_data), "\n")
cat("Regions:", paste(levels(enrich_data$Region), collapse = ", "), "\n")
cat("Amino acids:", length(unique(enrich_data$AA)), "\n")

# 5. PREPARE HEATMAP DATA ------------------------------------------------------
cat("\nPreparing heatmap data...\n")

# Calculate log2 enrichment factor
heatmap_data <- enrich_data %>%
  mutate(
    log2_EF = log2(Enrichment_Factor),
    # Define significance levels
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Add asterisk for FDR < 0.1 as dot for visualization
    Trend = ifelse(FDR < 0.1 & FDR >= 0.05, "¡¤", "")
  )

# Set limits for color scale
color_limits <- c(-2, 2)
heatmap_data$log2_EF <- pmax(pmin(heatmap_data$log2_EF, color_limits[2]), 
                             color_limits[1])

# 6. CREATE ENRICHMENT HEATMAP -------------------------------------------------
cat("Creating enrichment heatmap...\n")

# Define color palette for heatmap
heatmap_colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                    "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")

# Create heatmap plot
heatmap_plot <- ggplot(heatmap_data, 
                       aes(x = Region, y = AA, fill = log2_EF)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.9, height = 0.9) +
  # Add significance asterisks
  geom_text(aes(label = Significance), 
            size = 4, 
            vjust = 0.8, 
            color = "black",
            fontface = "bold") +
  # Add trend dots for borderline significance
  geom_text(aes(label = Trend), 
            size = 5, 
            vjust = 0.8, 
            color = "gray30",
            fontface = "bold") +
  # Color gradient
  scale_fill_gradientn(
    colors = heatmap_colors,
    limits = color_limits,
    na.value = "gray90",
    name = expression(Log[2]("Enrichment")),
    breaks = seq(-2, 2, by = 1),
    labels = c("¡Ü-2", "-1", "0", "1", "¡Ý2")
  ) +
  # Axis labels
  labs(
    x = "Protein Region",
    y = "Amino Acid"
  ) +
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 0, 
      hjust = 0.5, 
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),
    axis.title.x = element_text(
      face = "bold", 
      size = 12,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      face = "bold", 
      size = 12,
      margin = margin(r = 10)
    ),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = 1
    ),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10, "pt")
  ) +
  # Reverse amino acid order (alphabetical from top to bottom)
  scale_y_discrete(limits = rev(amino_acids)) +
  # Add background shading for amino acid groups
  annotate("rect", 
           xmin = 0.5, xmax = 3.5,
           ymin = c(15.5, 12.5, 10.5, 7.5, 5.5, 3.5, 2.5, 0.5),
           ymax = c(20.5, 15.5, 12.5, 10.5, 7.5, 5.5, 3.5, 2.5),
           fill = c("#F0F0F0", "#E0E0E0", "#F0F0F0", "#E0E0E0", 
                    "#F0F0F0", "#E0E0E0", "#F0F0F0", "#E0E0E0"),
           alpha = 0.3)

# 7. CREATE HEATMAP LEGEND -----------------------------------------------------
cat("Creating heatmap legend...\n")

heatmap_legend <- cowplot::get_legend(
  ggplot(heatmap_data, aes(x = Region, y = AA, fill = log2_EF)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = heatmap_colors,
      limits = color_limits,
      na.value = "gray90",
      name = expression(Log[2]("Enrichment")),
      breaks = seq(-2, 2, by = 1),
      labels = c("¡Ü-2", "-1", "0", "1", "¡Ý2")
    ) +
    theme(
      legend.position = "right",
      legend.key.height = unit(2, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.title = element_text(
        face = "bold", 
        size = 10,
        margin = margin(b = 5)
      ),
      legend.text = element_text(size = 9),
      legend.background = element_rect(
        fill = "white", 
        color = "black",
        size = 0.5
      )
    )
)

# 8. LOAD AND PREPARE PREFERENCE DATA ------------------------------------------
cat("\nLoading preference data from:", preference_file, "\n")

if (!file.exists(preference_file)) {
  stop(paste("Preference file not found:", preference_file))
}

# Read preference data
pref_data <- read.table(preference_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

# Ensure proper factor levels
pref_data$AA <- factor(pref_data$AA, levels = amino_acids)
pref_data$Region <- factor(pref_data$Region, 
                           levels = c("Core", "Surface", "Key"))

# Join with amino acid groups
pref_data <- pref_data %>%
  left_join(aa_groups, by = "AA") %>%
  mutate(
    Group = factor(Group, levels = rev(names(aa_classification)))
  )

cat("Preference data loaded. Dimensions:", dim(pref_data), "\n")

# 9. CREATE PREFERENCE LINE PLOT -----------------------------------------------
cat("Creating preference line plot...\n")

# Define region colors (matching heatmap)
region_colors <- c(
  "Core" = "#B2182B",    # Red
  "Surface" = "#2166AC", # Blue
  "Key" = "#4DAF4A"      # Green
)

# Create preference plot
pref_plot <- ggplot(pref_data, 
                    aes(x = AA, y = Preference, group = Region, color = Region)) +
  # Add horizontal grid lines
  geom_hline(yintercept = c(0.1, 0.2, 0.5, 1, 2, 5, 10), 
             color = "gray90", 
             size = 0.3,
             linetype = "dashed") +
  # Add reference line at y = 1
  geom_hline(yintercept = 1, 
             linetype = "solid", 
             color = "red", 
             size = 0.8,
             alpha = 0.7) +
  # Add lines and points
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  # Color scheme
  scale_color_manual(values = region_colors) +
  # Axis labels
  labs(
    x = "Amino Acid",
    y = "Preference Ratio (Log Scale)",
    color = "Region"
  ) +
  # Log scale for y-axis
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10")
  ) +
  # Flip coordinates
  coord_flip() +
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 9
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 12,
      margin = margin(t = 10)
    ),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_line(
      color = "gray90",
      size = 0.3
    ),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      size = 1
    ),
    legend.position = "top",
    legend.title = element_text(
      face = "bold", 
      size = 11,
      margin = margin(b = 5)
    ),
    legend.text = element_text(
      face = "bold", 
      size = 10
    ),
    legend.key = element_rect(
      fill = "white", 
      color = "black"
    ),
    legend.box.spacing = unit(0.5, "cm"),
    plot.margin = margin(10, 10, 10, 10, "pt")
  ) +
  # Amino acid order
  scale_x_discrete(limits = rev(amino_acids)) +
  # Add amino acid group labels
  annotate("text",
           x = c(18.5, 15.5, 13.5, 11.5, 9.5, 7.5, 5.5, 3.5),
           y = 0.07,
           label = rev(names(aa_classification)),
           angle = 0,
           hjust = 0,
           size = 3.5,
           fontface = "italic",
           color = "gray40")

# 10. CREATE COMBINED FIGURE ---------------------------------------------------
cat("\nCreating combined figure...\n")

# Create title
title_plot <- ggplot() +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = "Amino Acid Regional Enrichment and Preference Analysis",
           size = 6, 
           fontface = "bold",
           hjust = 0.5) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 5, "pt")
  )

# Create legend for the figure
figure_legend <- ggplot() +
  annotate("text",
           x = 0.1,
           y = 0.5,
           label = "Significance: *FDR < 0.05, **FDR < 0.01, ***FDR < 0.001",
           size = 3.5,
           hjust = 0) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 5, "pt")
  )

# Arrange plots using patchwork
combined_plot <- (title_plot / 
  (heatmap_plot + 
     plot_spacer() + 
     heatmap_legend + 
     plot_layout(widths = c(0.8, 0.05, 0.2))) / 
  pref_plot / 
  figure_legend) +
  plot_layout(
    heights = c(0.05, 0.6, 0.3, 0.05),
    ncol = 1
  ) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.01, 0.98)
  )

# 11. SAVE THE COMBINED FIGURE -------------------------------------------------
cat("Saving figures...\n")

# Save in multiple formats
save_plot <- function(filename, width, height, device = NULL) {
  filepath <- file.path(output_dir, filename)
  
  if (is.null(device)) {
    device <- tools::file_ext(filename)
  }
  
  ggsave(
    filename = filepath,
    plot = combined_plot,
    width = width,
    height = height,
    device = device,
    dpi = 300,
    bg = "white"
  )
  
  cat("Saved:", filename, "\n")
}

# High-resolution PDF for publications
save_plot("aa_combined_analysis.pdf", width = 12, height = 14)

# EPS for publications (requires Cairo)
if (require("Cairo")) {
  save_plot("aa_combined_analysis.eps", width = 12, height = 14, 
            device = cairo_ps)
}

# PNG for presentations
save_plot("aa_combined_analysis.png", width = 12, height = 14)

# Save individual components as well
ggsave(
  file.path(output_dir, "aa_enrichment_heatmap.pdf"),
  plot = heatmap_plot,
  width = 6,
  height = 8
)

ggsave(
  file.path(output_dir, "aa_preference_plot.pdf"),
  plot = pref_plot,
  width = 8,
  height = 6
)

# 12. SAVE VISUALIZATION DATA --------------------------------------------------
cat("\nSaving visualization data...\n")

# Save processed data for reproducibility
write.table(heatmap_data, 
            file.path(output_dir, "heatmap_processed_data.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(pref_data, 
            file.path(output_dir, "preference_processed_data.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 13. SAVE SESSION INFO --------------------------------------------------------
cat("\nSaving session information...\n")
sink(file.path(output_dir, "session_info.txt"))
cat("Visualization completed on:", date(), "\n\n")
cat("Input files:\n")
cat("  Enrichment data:", enrichment_file, "\n")
cat("  Preference data:", preference_file, "\n\n")
print(sessionInfo())
sink()

# 14. DISPLAY THE PLOT (OPTIONAL) ----------------------------------------------
# Uncomment the following lines to display the plot in R
# cat("\nDisplaying combined plot...\n")
# print(combined_plot)

# 15. COMPLETION MESSAGE -------------------------------------------------------
cat("\n" + strrep("=", 60) + "\n")
cat("VISUALIZATION COMPLETE!\n")
cat(strrep("=", 60) + "\n")
cat("Results saved in:", output_dir, "\n")
cat("\nFiles generated:\n")
cat("  1. aa_combined_analysis.pdf - Main combined figure (PDF)\n")
cat("  2. aa_combined_analysis.eps - Main combined figure (EPS)\n")
cat("  3. aa_combined_analysis.png - Main combined figure (PNG)\n")
cat("  4. aa_enrichment_heatmap.pdf - Individual heatmap\n")
cat("  5. aa_preference_plot.pdf - Individual preference plot\n")
cat("  6. heatmap_processed_data.txt - Processed heatmap data\n")
cat("  7. preference_processed_data.txt - Processed preference data\n")
cat("  8. session_info.txt - Reproducibility information\n")
cat("  9. visualization_log.txt - Complete analysis log\n")
cat(strrep("=", 60) + "\n")

# Close log file
sink()

# ==============================================================================
# USAGE INSTRUCTIONS:
# 1. Prepare input files:
#    - AA_region_enrichment_results.txt: Output from enrichment analysis
#    - region_aa_preference.txt: Output from preference analysis
#
# 2. Expected data formats:
#    - Enrichment file: Tab-delimited with columns: AA, Region, Enrichment_Factor, FDR
#    - Preference file: Tab-delimited with columns: AA, Region, Preference
#
# 3. Run the script:
#    source("aa_combined_visualization.R")
#
# 4. Results will be saved in 'aa_combined_visualization_results' directory
#
# CUSTOMIZATION OPTIONS:
# - Modify 'amino_acids' vector for different amino acid order
# - Adjust 'color_limits' in heatmap for different data ranges
# - Change 'region_colors' for different color schemes
# - Modify significance thresholds in 'Significance' calculation
#
# DEPENDENCIES:
# - tidyverse: Data manipulation and plotting
# - cowplot: Plot arrangement and legend extraction
# - patchwork: Advanced plot composition
# - ggtext: Enhanced text rendering (optional)
# - RColorBrewer: Color palettes
# ==============================================================================