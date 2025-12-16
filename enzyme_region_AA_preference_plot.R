# ================================================================================
# Regional Amino Acid Preference Analysis for Protein Structures
# Description: This script analyzes amino acid preferences in different protein
#              structural regions (Core, Surface, Key residues) including:
#              1. Amino acid composition distribution
#              2. Regional amino acid preferences
#              3. Enrichment analysis using Fisher's exact test
# ================================================================================

# 1. LOAD REQUIRED PACKAGES ----------------------------------------------------
# Install packages if not already installed
required_packages <- c("tidyverse", "corrplot", "ggpubr", "viridis", "gridExtra", "cowplot")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 2. SETUP AND CONFIGURATION ---------------------------------------------------
# Create output directory
output_dir <- "regional_aa_analysis_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set file paths (modify these according to your data structure)
data_dir <- "data"  # Directory containing input files
core_file <- file.path(data_dir, "core_aa90.txt")
surface_file <- file.path(data_dir, "surface_aa90.txt")
key_file <- file.path(data_dir, "key_aa.txt")

# Create log file
log_file <- file.path(output_dir, "analysis_log.txt")
sink(log_file)
cat("Regional Amino Acid Preference Analysis Log\n")
cat("Analysis started:", date(), "\n\n")

# 3. DEFINE AMINO ACID PROPERTIES ----------------------------------------------
# Standard 20 amino acids
aa_list <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# Amino acid properties for visualization
aa_properties <- data.frame(
  AA = aa_list,
  FullName = c("Alanine", "Cysteine", "Aspartic acid", "Glutamic acid", 
               "Phenylalanine", "Glycine", "Histidine", "Isoleucine", 
               "Lysine", "Leucine", "Methionine", "Asparagine", 
               "Proline", "Glutamine", "Arginine", "Serine", 
               "Threonine", "Valine", "Tryptophan", "Tyrosine"),
  Type = c("Aliphatic", "Sulfur", "Acidic", "Acidic", "Aromatic",
           "Aliphatic", "Basic", "Aliphatic", "Basic", "Aliphatic",
           "Sulfur", "Amide", "Cyclic", "Amide", "Basic", "Hydroxylic",
           "Hydroxylic", "Aliphatic", "Aromatic", "Aromatic"),
  Hydrophobicity = c(1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5,
                     -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8,
                     -0.7, 4.2, -0.9, -1.3)
)

# 4. DEFINE DATA READING FUNCTIONS --------------------------------------------
read_aa_count <- function(file_path, region_name) {
  """
  Reads amino acid count data from formatted text files
  
  Args:
    file_path: Path to input file
    region_name: Name of the protein region (e.g., "Core", "Surface")
  
  Returns:
    Data frame with amino acid counts
  """
  
  cat("Reading data from:", file_path, "\n")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  tryCatch({
    lines <- readLines(file_path)
    # Remove comment lines
    lines <- lines[!grepl("^#", lines)]  
    
    data_list <- list()
    
    for (line in lines) {
      parts <- str_split(line, "\\s+")[[1]]
      protein <- parts[1]
      
      # Initialize counts for all amino acids
      counts <- setNames(rep(0, length(aa_list)), aa_list) 
      
      # Parse amino acid counts
      for (item in parts[-1]) {
        if (grepl(":", item)) {
          aa_count <- str_split(item, ":")[[1]]
          aa <- aa_count[1]
          count <- as.integer(aa_count[2])
          
          if (aa %in% aa_list) {
            counts[aa] <- count
          }
        }
      }
      data_list[[protein]] <- counts
    }
    
    # Convert to tidy format
    result <- bind_rows(data_list, .id = "Protein") %>% 
      pivot_longer(-Protein, names_to = "AA", values_to = "Count") %>%
      mutate(Region = region_name)
    
    cat("Successfully read", length(data_list), "proteins for", region_name, "\n")
    return(result)
    
  }, error = function(e) {
    warning(paste("Error reading", file_path, ":", e$message))
    return(NULL)
  })
}

read_aa_sequence <- function(file_path, region_name) {
  """
  Reads amino acid sequences from tab-delimited files
  
  Args:
    file_path: Path to input file
    region_name: Name of the protein region
  
  Returns:
    Data frame with amino acid counts from sequences
  """
  
  cat("Reading sequence data from:", file_path, "\n")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  tryCatch({
    # Read sequence data
    data <- read_delim(
      file_path, 
      delim = "\t", 
      col_names = FALSE, 
      comment = "#", 
      trim_ws = TRUE, 
      show_col_types = FALSE
    )
    
    names(data)[1] <- "Protein"
    
    # Convert to long format and count amino acids
    long_data <- data %>% 
      pivot_longer(
        cols = -Protein, 
        names_to = "Position", 
        values_to = "AA", 
        values_drop_na = TRUE
      ) %>% 
      filter(
        AA != "" & !is.na(AA) & nchar(AA) == 1
      ) 
    
    # Count amino acids per protein
    count_data <- long_data %>% 
      group_by(Protein, AA) %>% 
      summarise(Count = n(), .groups = "drop")
    
    # Ensure all amino acids are represented
    complete_data <- complete(
      count_data, 
      Protein, 
      AA = aa_list, 
      fill = list(Count = 0)
    ) %>% 
      mutate(Region = region_name)
    
    cat("Successfully read", n_distinct(complete_data$Protein), 
        "proteins for", region_name, "\n")
    return(complete_data)
    
  }, error = function(e) {
    warning(paste("Error reading", file_path, ":", e$message))
    return(NULL)
  })
}

# 5. LOAD AND COMBINE DATA -----------------------------------------------------
cat("\nLoading data...\n")

# Read data for each region
core_df <- read_aa_count(core_file, "Core")
surface_df <- read_aa_count(surface_file, "Surface")
key_df <- read_aa_sequence(key_file, "Key")

# Check if all data loaded successfully
if (is.null(core_df) || is.null(surface_df) || is.null(key_df)) {
  stop("Failed to load one or more input files. Check file paths and formats.")
}

# Combine all data
combined <- bind_rows(core_df, surface_df, key_df) %>%
  mutate(Region = factor(Region, levels = c("Core", "Surface", "Key")))

cat("\nData Summary:\n")
cat("Total proteins:", n_distinct(combined$Protein), "\n")
cat("Total observations:", nrow(combined), "\n")
cat("Regions:", paste(unique(combined$Region), collapse = ", "), "\n")

# 6. CALCULATE AMINO ACID COMPOSITION ------------------------------------------
cat("\nCalculating amino acid composition...\n")

# Calculate percentages for each protein and region
combined <- combined %>%
  group_by(Protein, Region) %>%
  mutate(
    Total = sum(Count),
    Percentage = ifelse(Total > 0, Count / Total * 100, 0)
  ) %>%
  ungroup()

# Save composition data
composition_file <- file.path(output_dir, "region_aa_composition.txt")
write_tsv(combined, composition_file)
cat("Composition data saved to:", composition_file, "\n")

# 7. CREATE AMINO ACID DISTRIBUTION PLOT ---------------------------------------
cat("\nCreating amino acid distribution plot...\n")

# Calculate regional averages
region_avg <- combined %>%
  group_by(Region, AA) %>%
  summarise(
    AvgCount = mean(Count, na.rm = TRUE),
    AvgPercentage = mean(Percentage, na.rm = TRUE),
    SEM = sd(Percentage, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(aa_properties, by = "AA") %>%
  mutate(AA = factor(AA, levels = aa_list))

# Save distribution data
distribution_file <- file.path(output_dir, "region_aa_distribution.txt")
write_tsv(region_avg, distribution_file)
cat("Distribution data saved to:", distribution_file, "\n")

# Create distribution plot
p_distribution <- ggplot(region_avg, aes(x = AA, y = AvgPercentage, fill = Region)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = AvgPercentage - SEM, ymax = AvgPercentage + SEM),
    position = position_dodge(width = 0.8),
    width = 0.2,
    size = 0.5
  ) +
  labs(
    title = "Amino Acid Distribution by Protein Region",
    subtitle = "Average percentage composition with standard error",
    y = "Average Percentage (%)", 
    x = "Amino Acid",
    fill = "Region"
  ) +
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      vjust = 1, 
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# 8. CREATE AMINO ACID PREFERENCE PLOT ----------------------------------------
cat("\nCalculating amino acid preferences...\n")

# Calculate global average for each amino acid
global_avg <- region_avg %>%
  group_by(AA) %>%
  summarise(
    GlobalAvg = mean(AvgPercentage),
    .groups = "drop"
  )

# Calculate preferences
preference <- region_avg %>%
  left_join(global_avg, by = "AA") %>%
  mutate(
    Preference = AvgPercentage / GlobalAvg,
    Log2Preference = log2(Preference)
  )

# Save preference data
preference_file <- file.path(output_dir, "region_aa_preference.txt")
write_tsv(preference, preference_file)
cat("Preference data saved to:", preference_file, "\n")

# Create preference plot
p_preference <- ggplot(preference, aes(x = AA, y = Log2Preference, color = Region, group = Region)) +
  geom_line(size = 1.2, alpha = 0.7) + 
  geom_point(size = 3, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  labs(
    title = "Regional Amino Acid Preferences",
    subtitle = "Log2 fold-change compared to global average",
    y = expression(Log[2]("Region/Global")),
    x = "Amino Acid",
    color = "Region"
  ) +
  scale_color_viridis(discrete = TRUE, begin = 0.2, end = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      vjust = 1, 
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# 9. PERFORM ENRICHMENT ANALYSIS ----------------------------------------------
cat("\nPerforming enrichment analysis...\n")

# Amino acid full name mapping
aa_mapping <- c(
  "A" = "Ala", "C" = "Cys", "D" = "Asp", "E" = "Glu", "F" = "Phe", 
  "G" = "Gly", "H" = "His", "I" = "Ile", "K" = "Lys", "L" = "Leu", 
  "M" = "Met", "N" = "Asn", "P" = "Pro", "Q" = "Gln", "R" = "Arg", 
  "S" = "Ser", "T" = "Thr", "V" = "Val", "W" = "Trp", "Y" = "Tyr"
)

# Convert amino acid symbols to full names
all_data <- combined %>%
  mutate(AA_full = aa_mapping[as.character(AA)])

# Calculate background frequencies
total_background <- sum(all_data$Count)
background <- all_data %>%
  group_by(AA_full) %>%
  summarise(
    Background_count = sum(Count),
    .groups = "drop"
  ) %>%
  mutate(
    Total_background = total_background,
    Percentage_background = Background_count / total_background * 100
  )

# Calculate regional statistics
region_stats <- all_data %>%
  group_by(Region) %>%
  mutate(Total_region = sum(Count)) %>%
  group_by(Region, AA_full) %>%
  summarise(
    Count_in_region = sum(Count),
    Total_region = first(Total_region),
    .groups = "drop"
  ) %>%
  mutate(
    Percentage_in_region = Count_in_region / Total_region * 100
  )

# Define function for Fisher's exact test
run_fisher_test <- function(count_in_region, background_count, total_region, total_background) {
  target_in_other <- background_count - count_in_region
  other_in_region <- total_region - count_in_region
  other_in_other <- total_background - total_region - target_in_other
  
  if (any(c(target_in_other, other_in_region, other_in_other) < 0)) {
    return(NA)
  }
  
  contingency_table <- matrix(
    c(count_in_region, target_in_other, other_in_region, other_in_other),
    nrow = 2,
    byrow = TRUE
  )
  
  fisher.test(contingency_table, alternative = "greater")$p.value
}

# Perform enrichment analysis
enrich_data <- region_stats %>%
  left_join(background, by = "AA_full") %>%
  mutate(
    Enrichment_Factor = (Count_in_region / Total_region) / (Background_count / Total_background),
    Log2_Enrichment = log2(Enrichment_Factor)
  ) %>%
  mutate(
    P_value = pmap_dbl(
      list(Count_in_region, Background_count, Total_region, Total_background),
      run_fisher_test
    )
  ) %>%
  mutate(
    FDR = p.adjust(P_value, method = "fdr"),
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(
    AA = AA_full, Region, Count_in_region, Total_region, Percentage_in_region,
    Background_count, Total_background, Percentage_background,
    Enrichment_Factor, Log2_Enrichment, P_value, FDR, Significance
  ) %>%
  arrange(Region, AA)

# Save enrichment results
enrichment_file <- file.path(output_dir, "aa_region_enrichment_results.txt")
write.table(enrich_data, enrichment_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Enrichment analysis results saved to:", enrichment_file, "\n")

# 10. CREATE ENRICHMENT HEATMAP -----------------------------------------------
cat("\nCreating enrichment heatmap...\n")

# Prepare data for heatmap
heatmap_data <- enrich_data %>%
  select(AA, Region, Log2_Enrichment) %>%
  pivot_wider(names_from = Region, values_from = Log2_Enrichment) %>%
  column_to_rownames("AA") %>%
  as.matrix()

# Create heatmap
pdf(file.path(output_dir, "enrichment_heatmap.pdf"), width = 8, height = 10)
corrplot(
  heatmap_data,
  method = "color",
  is.corr = FALSE,
  col = colorRampPalette(c("blue", "white", "red"))(100),
  tl.col = "black",
  tl.srt = 45,
  cl.pos = "r",
  cl.ratio = 0.1,
  cl.length = 5,
  title = "Log2 Enrichment of Amino Acids by Region",
  mar = c(0, 0, 2, 0)
)
dev.off()

# 11. COMBINE AND SAVE PLOTS --------------------------------------------------
cat("\nSaving plots...\n")

# Save individual plots
ggsave(
  file.path(output_dir, "aa_distribution_plot.pdf"),
  plot = p_distribution,
  width = 14,
  height = 8,
  device = "pdf"
)

ggsave(
  file.path(output_dir, "aa_preference_plot.pdf"),
  plot = p_preference,
  width = 14,
  height = 8,
  device = "pdf"
)

# Create combined plot
combined_plot <- plot_grid(
  p_distribution + theme(legend.position = "none"),
  p_preference + theme(legend.position = "none"),
  ncol = 2,
  align = "h",
  labels = c("A", "B"),
  label_size = 16
)

# Extract legend
legend <- get_legend(p_distribution + theme(legend.box.margin = margin(0, 0, 0, 12)))

# Add legend to combined plot
final_plot <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(10, 1)
)

# Save combined plot
ggsave(
  file.path(output_dir, "combined_aa_analysis_plot.pdf"),
  plot = final_plot,
  width = 18,
  height = 9,
  device = "pdf"
)

# 12. GENERATE SUMMARY REPORT -------------------------------------------------
cat("\nGenerating summary report...\n")

summary_stats <- combined %>%
  group_by(Region) %>%
  summarise(
    Total_Proteins = n_distinct(Protein),
    Total_Residues = sum(Count),
    Unique_Amino_Acids = n_distinct(AA[Count > 0]),
    .groups = "drop"
  )

# Save summary
write_tsv(summary_stats, file.path(output_dir, "analysis_summary.txt"))

# 13. SAVE SESSION INFO -------------------------------------------------------
cat("\nSaving session information...\n")
sink(file.path(output_dir, "session_info.txt"))
cat("Analysis completed on:", date(), "\n\n")
cat("Input files:\n")
cat("  Core data:", core_file, "\n")
cat("  Surface data:", surface_file, "\n")
cat("  Key residue data:", key_file, "\n\n")
print(sessionInfo())
sink()

# 14. COMPLETION MESSAGE ------------------------------------------------------
cat("\n" + strrep("=", 60) + "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 60) + "\n")
cat("Results saved in:", output_dir, "\n")
cat("Files generated:\n")
cat("  1. region_aa_composition.txt - Raw composition data\n")
cat("  2. region_aa_distribution.txt - Distribution statistics\n")
cat("  3. region_aa_preference.txt - Preference calculations\n")
cat("  4. aa_region_enrichment_results.txt - Enrichment analysis\n")
cat("  5. aa_distribution_plot.pdf - Distribution bar plot\n")
cat("  6. aa_preference_plot.pdf - Preference line plot\n")
cat("  7. enrichment_heatmap.pdf - Enrichment heatmap\n")
cat("  8. combined_aa_analysis_plot.pdf - Combined figure\n")
cat("  9. analysis_summary.txt - Summary statistics\n")
cat("  10. session_info.txt - Reproducibility information\n")
cat(strrep("=", 60) + "\n")

# Close log file
sink()

# ==============================================================================
# USAGE INSTRUCTIONS:
# 1. Prepare input files in the 'data' directory:
#    - core_aa90.txt: Core residue amino acid counts
#    - surface_aa90.txt: Surface residue amino acid counts  
#    - key_aa.txt: Key residue sequences
#
# 2. File formats:
#    - For core_aa90.txt and surface_aa90.txt:
#      ProteinID A:5 C:3 D:2 ...
#    - For key_aa.txt (tab-delimited):
#      ProteinID  A  C  D  E  ...
#      Protein1   M  L  V  I  ...
#
# 3. Run the analysis:
#    source("regional_amino_acid_analysis.R")
#
# 4. Results will be saved in 'regional_aa_analysis_results' directory
#
# CUSTOMIZATION:
# - Modify aa_list to include non-standard amino acids
# - Adjust color schemes in scale_fill_viridis() and scale_color_viridis()
# - Change statistical thresholds in enrichment analysis
# ==============================================================================