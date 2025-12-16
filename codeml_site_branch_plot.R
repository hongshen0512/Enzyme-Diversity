# ================================================================================
# PAML CodeML Branch-Site Model Positive Selection Visualization
# Author: [Your Name]
# Date: [Date]
# Description: This script visualizes positive selection sites identified by
#              PAML CodeML branch-site model analysis. Creates Nature-style
#              lollipop plots for posterior probability distributions across
#              evolutionary branches.
# ================================================================================

# 1. LOAD REQUIRED PACKAGES ----------------------------------------------------
# Install packages if not already installed
required_packages <- c("ggplot2", "grid", "gridExtra", "dplyr", 
                       "stringr", "scales", "RColorBrewer", "cowplot", "viridis")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 2. SETUP AND CONFIGURATION ---------------------------------------------------
# Create output directory
output_dir <- "codeml_positive_selection_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set file paths
input_file <- "data/results.txt"  # PAML CodeML output file

# Create log file
log_file <- file.path(output_dir, "analysis_log.txt")
sink(log_file)
cat("PAML CodeML Positive Selection Analysis Log\n")
cat("Analysis started:", date(), "\n\n")

# 3. DATA PARSING FUNCTIONS ----------------------------------------------------
parse_codeml_results <- function(file_path) {
  """
  Parse PAML CodeML branch-site model results from output file
  
  Args:
    file_path: Path to CodeML output file
  
  Returns:
    List containing site data for each branch/cluster
  """
  
  cat("Reading PAML CodeML results from:", file_path, "\n")
  
  if (!file.exists(file_path)) {
    stop(paste("Input file not found:", file_path))
  }
  
  tryCatch({
    data <- readLines(file_path)
    
    sites <- list()
    current_branch <- NULL
    dnds_value <- NULL
    sites_df <- data.frame()
    
    # Parse file line by line
    for (line in data) {
      if (line == "") next
      
      parts <- str_split(line, "\t")[[1]]
      if (length(parts) < 3) next
      
      branch_name <- parts[1]
      
      # Check if this is a dN/dS line
      if (parts[2] == "dN/dS") {
        # Save current branch data if exists
        if (!is.null(current_branch) && nrow(sites_df) > 0) {
          sites[[current_branch]] <- list(
            dnds = as.numeric(dnds_value),
            sites = sites_df
          )
        }
        
        # Start new branch
        current_branch <- branch_name
        dnds_value <- parts[3]
        sites_df <- data.frame()
        
      } else if (grepl("NEB_P", parts[2])) {
        # Parse site information line
        site_info <- parts[2]
        
        # Extract position number
        position_str <- str_extract(site_info, "\\d+")
        if (is.na(position_str)) next
        position <- as.numeric(position_str)
        
        # Extract amino acid residue
        aa_part <- str_remove_all(site_info, "\\d+|\\s+|-|NEB_P|\\*")
        aa <- str_sub(aa_part, 1, 1)
        if (is.na(aa) || aa == "") aa <- "X"
        
        # Check for significance asterisk
        significant <- grepl("\\*", site_info)
        
        # Posterior probability value
        prob <- as.numeric(parts[3])
        
        # Add to current branch data frame
        new_row <- data.frame(
          position = position,
          amino_acid = aa,
          posterior_prob = prob,
          significant = significant,
          stringsAsFactors = FALSE
        )
        
        sites_df <- rbind(sites_df, new_row)
      }
    }
    
    # Save the last branch data
    if (!is.null(current_branch) && nrow(sites_df) > 0) {
      sites[[current_branch]] <- list(
        dnds = as.numeric(dnds_value),
        sites = sites_df
      )
    }
    
    cat("Successfully parsed", length(sites), "branches/clusters\n")
    return(sites)
    
  }, error = function(e) {
    stop(paste("Error parsing CodeML results:", e$message))
  })
}

# 4. AMINO ACID COLOR SCHEME ---------------------------------------------------
create_amino_acid_color_scheme <- function() {
  """
  Create a consistent, publication-ready color scheme for amino acids
  
  Returns:
    Named vector of colors for each amino acid
  """
  
  # 20 standard amino acids (one-letter codes)
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Define amino acid groups for color grouping
  aa_groups <- list(
    "Aliphatic" = c("A", "G", "I", "L", "V"),
    "Aromatic" = c("F", "W", "Y"),
    "Acidic" = c("D", "E"),
    "Basic" = c("H", "K", "R"),
    "Hydroxylic" = c("S", "T"),
    "Sulfur" = c("C", "M"),
    "Amide" = c("N", "Q"),
    "Cyclic" = c("P")
  )
  
  # Color palette by amino acid group
  group_colors <- c(
    "#1B9E77",  # Aliphatic - Teal
    "#D95F02",  # Aromatic - Orange
    "#7570B3",  # Acidic - Purple
    "#E7298A",  # Basic - Magenta
    "#66A61E",  # Hydroxylic - Green
    "#E6AB02",  # Sulfur - Gold
    "#A6761D",  # Amide - Brown
    "#666666"   # Cyclic - Gray
  )
  
  # Create color mapping
  aa_color_mapping <- setNames(rep(NA, length(amino_acids)), amino_acids)
  
  # Assign colors based on groups
  for (i in seq_along(aa_groups)) {
    group_name <- names(aa_groups)[i]
    group_aas <- aa_groups[[i]]
    aa_color_mapping[group_aas] <- group_colors[i]
  }
  
  # Verify all amino acids have colors
  missing_aas <- amino_acids[is.na(aa_color_mapping)]
  if (length(missing_aas) > 0) {
    warning("Missing colors for amino acids: ", paste(missing_aas, collapse = ", "))
    # Assign gray to missing amino acids
    aa_color_mapping[missing_aas] <- "#999999"
  }
  
  return(aa_color_mapping)
}

# 5. VISUALIZATION FUNCTIONS ---------------------------------------------------
create_branch_lollipop_plot <- function(branch_name, branch_data, aa_colors, 
                                        total_sites = 463, 
                                        prob_threshold = 0.95) {
  """
  Create Nature-style lollipop plot for a single branch
  
  Args:
    branch_name: Name of the evolutionary branch/cluster
    branch_data: List containing dnds and site data
    aa_colors: Amino acid color mapping
    total_sites: Total number of sites analyzed
    prob_threshold: Posterior probability threshold for significance
  
  Returns:
    ggplot object
  """
  
  # Check if data is available
  if (is.null(branch_data) || nrow(branch_data$sites) == 0) {
    # Create empty plot for branches with no data
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0(branch_name, "\nNo significant sites\n(", 
                             total_sites, " sites analyzed)"),
               size = 4, color = "#666666") +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", color = NA))
    return(p)
  }
  
  # Extract significant sites
  sig_sites <- branch_data$sites %>% 
    filter(significant) %>% 
    arrange(position)
  
  # Handle cases with no significant sites
  if (nrow(sig_sites) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(branch_name, "\nNo sites with posterior probability ¡Ý", 
                             prob_threshold, "\n(", total_sites, " sites analyzed)"),
               size = 4, color = "#666666") +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", color = NA))
    return(p)
  }
  
  # Filter for amino acids with valid colors
  sig_sites <- sig_sites %>%
    filter(amino_acid %in% names(aa_colors))
  
  if (nrow(sig_sites) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(branch_name, "\nNo sites with valid amino acid codes\n(", 
                             total_sites, " sites analyzed)"),
               size = 4, color = "#666666") +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", color = NA))
    return(p)
  }
  
  # Limit the number of sites shown for readability
  max_sites_to_show <- 50
  if (nrow(sig_sites) > max_sites_to_show) {
    sig_sites <- sig_sites %>%
      arrange(desc(posterior_prob)) %>%
      head(max_sites_to_show) %>%
      arrange(position)
  }
  
  # Calculate statistics
  sig_count <- nrow(sig_sites)
  percent_sig <- round(sig_count / total_sites * 100, 1)
  dnds_value <- round(branch_data$dnds, 3)
  
  # Create the lollipop plot
  p <- ggplot(sig_sites, aes(x = position, y = posterior_prob)) +
    # Add significance threshold line
    geom_hline(yintercept = prob_threshold, 
               linetype = "dashed", 
               color = "#E74C3C", 
               alpha = 0.7, 
               linewidth = 0.5) +
    # Add vertical segments (lollipop stems)
    geom_segment(aes(x = position, xend = position, y = 0, yend = posterior_prob),
                 color = "#2C3E50", 
                 alpha = 0.6, 
                 linewidth = 0.4) +
    # Add points colored by amino acid
    geom_point(aes(fill = amino_acid, color = amino_acid), 
               shape = 21, 
               size = 3, 
               stroke = 0.5,
               alpha = 0.8) +
    # Amino acid color scale
    scale_fill_manual(values = aa_colors) +
    scale_color_manual(values = aa_colors) +
    # Axis labels and plot title
    labs(
      title = paste0(branch_name, 
                    ": ", sig_count, " sites with posterior probability ¡Ý", prob_threshold,
                    "\n¦Ø = ", dnds_value, " (", sig_count, "/", total_sites, 
                    " sites, ", percent_sig, "%)"),
      x = "Position",
      y = "Posterior probability"
    ) +
    # Nature-style theme
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 11, 
        face = "bold", 
        hjust = 0.5,
        margin = margin(b = 8),
        lineheight = 1.2
      ),
      axis.title = element_text(
        size = 10, 
        face = "plain", 
        color = "black"
      ),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.text = element_text(
        size = 8, 
        color = "black"
      ),
      axis.line = element_line(
        color = "black", 
        linewidth = 0.35
      ),
      axis.ticks = element_line(
        color = "black", 
        linewidth = 0.35
      ),
      axis.ticks.length = unit(1.5, "pt"),
      panel.grid.major.y = element_line(
        color = "gray92", 
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black", 
        fill = NA, 
        linewidth = 0.5
      ),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 8, 8, 8),
      legend.position = "none"  # Legend will be separate
    ) +
    # Y-axis settings
    scale_y_continuous(
      limits = c(0, 1.05), 
      breaks = seq(0, 1, by = 0.2),
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    # X-axis settings
    scale_x_continuous(
      expand = expansion(mult = 0.03)
    )
  
  return(p)
}

create_amino_acid_legend <- function(aa_colors) {
  """
  Create a compact, publication-ready amino acid color legend
  
  Args:
    aa_colors: Named vector of amino acid colors
  
  Returns:
    ggplot object for legend
  """
  
  # Create data frame for legend
  aa_data <- data.frame(
    amino_acid = names(aa_colors),
    color = aa_colors,
    group = NA,
    stringsAsFactors = FALSE
  )
  
  # Define amino acid groups for organization
  group_info <- list(
    "Aliphatic" = c("A", "G", "I", "L", "V"),
    "Aromatic" = c("F", "W", "Y"),
    "Acidic" = c("D", "E"),
    "Basic" = c("H", "K", "R"),
    "Hydroxylic" = c("S", "T"),
    "Sulfur" = c("C", "M"),
    "Amide" = c("N", "Q"),
    "Cyclic" = c("P")
  )
  
  # Assign groups
  for (group_name in names(group_info)) {
    group_aas <- group_info[[group_name]]
    aa_data$group[aa_data$amino_acid %in% group_aas] <- group_name
  }
  
  # Order by group
  aa_data <- aa_data[order(match(aa_data$group, names(group_info))), ]
  
  # Create legend plot
  p <- ggplot(aa_data, aes(x = 1, y = seq_len(nrow(aa_data)))) +
    geom_point(aes(color = amino_acid, fill = amino_acid),
               shape = 21, 
               size = 4, 
               stroke = 0.5) +
    geom_text(aes(label = amino_acid),
              hjust = 0, 
              nudge_x = 0.2, 
              size = 3.5,
              color = "black") +
    scale_color_manual(values = aa_colors) +
    scale_fill_manual(values = aa_colors) +
    labs(title = "Amino Acid Key") +
    xlim(0.8, 2.5) +
    ylim(0, nrow(aa_data) + 1) +
    theme_void() +
    theme(
      plot.title = element_text(
        size = 10, 
        face = "bold", 
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.margin = margin(10, 10, 10, 10),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# 6. STATISTICAL SUMMARIES -----------------------------------------------------
generate_statistical_summary <- function(sites_data, branches_to_analyze) {
  """
  Generate comprehensive statistical summary of positive selection analysis
  
  Args:
    sites_data: Parsed CodeML results
    branches_to_analyze: Vector of branch names to include
  
  Returns:
    Data frame with summary statistics
  """
  
  summary_list <- list()
  
  for (branch in branches_to_analyze) {
    if (branch %in% names(sites_data)) {
      data <- sites_data[[branch]]
      sig_sites <- data$sites %>% filter(significant)
      sig_count <- nrow(sig_sites)
      dnds_value <- data$dnds
      
      # Calculate additional statistics if there are significant sites
      if (sig_count > 0) {
        avg_prob <- mean(sig_sites$posterior_prob, na.rm = TRUE)
        max_prob <- max(sig_sites$posterior_prob, na.rm = TRUE)
        min_position <- min(sig_sites$position, na.rm = TRUE)
        max_position <- max(sig_sites$position, na.rm = TRUE)
        
        # Amino acid composition of significant sites
        aa_composition <- table(sig_sites$amino_acid)
        most_common_aa <- names(which.max(aa_composition))
        aa_diversity <- length(unique(sig_sites$amino_acid))
        
      } else {
        avg_prob <- NA
        max_prob <- NA
        min_position <- NA
        max_position <- NA
        most_common_aa <- NA
        aa_diversity <- 0
      }
      
      summary_list[[branch]] <- data.frame(
        Branch = branch,
        dN_dS = dnds_value,
        Total_Sites_Analyzed = 463,  # Fixed as per original analysis
        Significant_Sites = sig_count,
        Percent_Significant = round(sig_count / 463 * 100, 2),
        Average_Posterior_Probability = round(avg_prob, 3),
        Maximum_Posterior_Probability = round(max_prob, 3),
        Position_Range = ifelse(!is.na(min_position) && !is.na(max_position),
                               paste(min_position, "-", max_position),
                               "N/A"),
        Most_Common_Amino_Acid = most_common_aa,
        Amino_Acid_Diversity = aa_diversity,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(summary_list) > 0) {
    return(do.call(rbind, summary_list))
  } else {
    return(NULL)
  }
}

# 7. MAIN ANALYSIS FUNCTION ----------------------------------------------------
run_codeml_visualization <- function() {
  """
  Main function to run the complete PAML CodeML visualization pipeline
  """
  
  cat("Starting PAML CodeML positive selection visualization...\n")
  
  # Parse PAML results
  sites_data <- parse_codeml_results(input_file)
  
  # Display parsed data summary
  cat("\nParsed Data Summary:\n")
  cat("=====================\n")
  for (branch in names(sites_data)) {
    sig_count <- sum(sites_data[[branch]]$sites$significant)
    dnds_value <- sites_data[[branch]]$dnds
    cat(sprintf("%s: %d significant sites, ¦Ø = %.3f\n", 
                branch, sig_count, dnds_value))
  }
  
  # Create amino acid color scheme
  cat("\nCreating amino acid color scheme...\n")
  aa_colors <- create_amino_acid_color_scheme()
  
  # Define branches to visualize (Cluster 3, 4, 6 only - excluding Cluster 5)
  branches_to_plot <- c("Cluster3", "Cluster4", "Cluster6")
  
  # Filter for available branches
  available_branches <- intersect(branches_to_plot, names(sites_data))
  
  if (length(available_branches) == 0) {
    cat("Warning: None of the specified branches were found in the data.\n")
    available_branches <- names(sites_data)
  }
  
  cat(sprintf("\nVisualizing branches: %s\n", paste(available_branches, collapse = ", ")))
  cat("Note: Cluster 5 is excluded from visualization as per original analysis\n")
  
  # Create individual lollipop plots
  branch_plots <- list()
  
  for (branch in available_branches) {
    cat(sprintf("Creating plot for %s...\n", branch))
    
    plot <- create_branch_lollipop_plot(
      branch_name = branch,
      branch_data = sites_data[[branch]],
      aa_colors = aa_colors,
      total_sites = 463,  # Fixed total sites as per original
      prob_threshold = 0.95
    )
    
    branch_plots[[branch]] <- plot
    
    # Save individual plot
    ggsave(
      file.path(output_dir, paste0("positive_selection_", tolower(branch), ".pdf")),
      plot = plot,
      width = 6,
      height = 5,
      device = "pdf"
    )
  }
  
  # Create amino acid legend
  cat("Creating amino acid legend...\n")
  aa_legend <- create_amino_acid_legend(aa_colors)
  
  # Save legend separately
  ggsave(
    file.path(output_dir, "amino_acid_legend.pdf"),
    plot = aa_legend,
    width = 3,
    height = 8,
    device = "pdf"
  )
  
  # Create combined figure
  cat("Creating combined figure...\n")
  
  # Create title and footer
  title_grob <- textGrob(
    "Positive Selection Sites Identified by Branch-Site Model Analysis",
    gp = gpar(fontsize = 16, fontface = "bold"),
    vjust = 1
  )
  
  footer_grob <- textGrob(
    "Analysis based on PAML CodeML branch-site model with posterior probability ¡Ý0.95\nTotal sites analyzed per branch: 463",
    gp = gpar(fontsize = 10),
    vjust = -0.5
  )
  
  # Arrange plots in grid (1x3 horizontal layout)
  plot_list <- branch_plots[available_branches]
  
  # Create combined plot using cowplot
  combined_plot <- plot_grid(
    plotlist = plot_list,
    ncol = 3,
    align = "h",
    labels = available_branches,
    label_size = 12,
    label_fontface = "bold"
  )
  
  # Add title and footer
  final_plot <- plot_grid(
    title_grob,
    combined_plot,
    footer_grob,
    ncol = 1,
    rel_heights = c(0.1, 1, 0.1)
  )
  
  # Save combined plot in multiple formats
  cat("Saving plots...\n")
  
  save_plot_formats <- function(plot_obj, filename_base, width = 14, height = 6) {
    # PDF
    ggsave(
      file.path(output_dir, paste0(filename_base, ".pdf")),
      plot = plot_obj,
      width = width,
      height = height,
      device = "pdf",
      dpi = 300
    )
    
    # PNG
    ggsave(
      file.path(output_dir, paste0(filename_base, ".png")),
      plot = plot_obj,
      width = width,
      height = height,
      device = "png",
      dpi = 300,
      bg = "white"
    )
    
    # EPS (if Cairo available)
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
  }
  
  # Save combined plot
  save_plot_formats(final_plot, "combined_positive_selection_plots", 14, 6)
  
  # Generate statistical summary
  cat("Generating statistical summary...\n")
  stats_summary <- generate_statistical_summary(sites_data, available_branches)
  
  if (!is.null(stats_summary)) {
    # Save summary to CSV
    write.csv(
      stats_summary,
      file.path(output_dir, "positive_selection_summary.csv"),
      row.names = FALSE
    )
    
    # Print summary to console
    cat("\nStatistical Summary:\n")
    cat("====================\n")
    print(stats_summary)
  }
  
  # Generate detailed site list
  cat("\nGenerating detailed site lists...\n")
  
  for (branch in available_branches) {
    if (branch %in% names(sites_data)) {
      sig_sites <- sites_data[[branch]]$sites %>% 
        filter(significant) %>%
        arrange(desc(posterior_prob))
      
      if (nrow(sig_sites) > 0) {
        # Save detailed site list
        write.csv(
          sig_sites,
          file.path(output_dir, paste0("significant_sites_", tolower(branch), ".csv")),
          row.names = FALSE
        )
      }
    }
  }
  
  # Save session information
  cat("\nSaving session information...\n")
  sink(file.path(output_dir, "session_info.txt"))
  cat("PAML CodeML Positive Selection Visualization\n")
  cat("Analysis completed on:", date(), "\n\n")
  cat("Input file:", input_file, "\n")
  cat("Branches analyzed:", paste(available_branches, collapse = ", "), "\n")
  cat("Total sites per branch: 463\n")
  cat("Posterior probability threshold: ¡Ý0.95\n\n")
  print(sessionInfo())
  sink()
  
  # Display completion message
  cat("\n" + strrep("=", 70) + "\n")
  cat("ANALYSIS COMPLETE!\n")
  cat(strrep("=", 70) + "\n")
  cat("Results saved in:", output_dir, "\n\n")
  cat("Files generated:\n")
  cat("  1. combined_positive_selection_plots.pdf/png/eps - Main combined figure\n")
  cat("  2. positive_selection_[branch].pdf - Individual branch plots\n")
  cat("  3. amino_acid_legend.pdf - Amino acid color legend\n")
  cat("  4. positive_selection_summary.csv - Statistical summary\n")
  cat("  5. significant_sites_[branch].csv - Detailed site lists\n")
  cat("  6. session_info.txt - Reproducibility information\n")
  cat("  7. analysis_log.txt - Complete analysis log\n")
  cat(strrep("=", 70) + "\n")
  
  # Close log file
  sink()
}

# 8. EXECUTE ANALYSIS ----------------------------------------------------------
# Uncomment the following line to run the analysis
# run_codeml_visualization()

# Or run with error handling
tryCatch({
  run_codeml_visualization()
}, error = function(e) {
  cat("Error during analysis:", e$message, "\n")
})

# ==============================================================================
# USAGE INSTRUCTIONS:
# 1. Prepare PAML CodeML output file:
#    - Run CodeML with branch-site model
#    - Save results in a text file (e.g., results.txt)
#    - Expected format: Tab-delimited with site-specific probabilities
#
# 2. Set the input file path (line 32):
#    input_file <- "data/results.txt"
#
# 3. Run the analysis:
#    source("codeml_positive_selection_analysis.R")
#    run_codeml_visualization()
#
# 4. Results will be saved in 'codeml_positive_selection_results' directory
#
# CUSTOMIZATION OPTIONS:
# - Modify branches_to_plot to include different branches
# - Adjust total_sites parameter if different from 463
# - Change prob_threshold for different significance levels
# - Modify color scheme in create_amino_acid_color_scheme()
#
# DATA FORMAT:
# The input file should contain:
# - Branch names (e.g., Cluster3, Cluster4, Cluster6)
# - dN/dS values for each branch
# - Site-specific posterior probabilities (NEB_P values)
# - Asterisks indicating statistical significance
#
# OUTPUT INTERPRETATION:
# - ¦Ø (dN/dS): >1 indicates positive selection
# - Posterior probability: Probability that site is under positive selection
# - Sites with probability ¡Ý0.95 are considered statistically significant
# ==============================================================================

# 9. HELPER FUNCTION FOR BATCH PROCESSING --------------------------------------
# (Optional) For processing multiple CodeML result files

process_multiple_codeml_files <- function(file_pattern = "codeml_results_*.txt") {
  """
  Process multiple PAML CodeML output files
  
  Args:
    file_pattern: Pattern to match CodeML result files
  """
  
  result_files <- list.files("data", pattern = file_pattern, full.names = TRUE)
  
  if (length(result_files) == 0) {
    cat("No files found matching pattern:", file_pattern, "\n")
    return()
  }
  
  for (file in result_files) {
    cat("\nProcessing file:", basename(file), "\n")
    
    # Create subdirectory for this file's results
    file_output_dir <- file.path(output_dir, tools::file_path_sans_ext(basename(file)))
    if (!dir.exists(file_output_dir)) {
      dir.create(file_output_dir)
    }
    
    # Temporarily modify input file
    original_input_file <- input_file
    input_file <<- file
    
    # Run visualization
    tryCatch({
      run_codeml_visualization()
    }, error = function(e) {
      cat("Error processing", basename(file), ":", e$message, "\n")
    })
    
    # Restore original input file
    input_file <<- original_input_file
  }
}