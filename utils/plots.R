##################################################
#################### PCA Plot ####################
##################################################

plot_pca <- function(dds, intgroup="time_point", ntop=500) {
  vsd <- vst(dds, blind=FALSE) # Get variance stabilized data
  pca_data <- plotPCA(vsd, intgroup=intgroup, returnData=TRUE, ntop=ntop) # Run PCA
  percentVar <- round(100 * attr(pca_data, "percentVar")) # Percentage of variance explained
  # Plot
  ggplot(pca_data, aes(x=PC1, y=PC2, color=time_point, shape=replicate)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of RNA-seq samples") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_color_brewer(palette = "Set1")
}

##################################################
############# Counts Histogram Plot ##############
##################################################
plot_count_distributions <- function(results_df, 
                                     metrics = c("baseMean", "lfcSE", "padj", "log2FoldChange"),
                                     labels = c("Base Mean", "LFC Standard Error", 
                                                "Adjusted P-value", "Log2 Fold Change"),
                                     colors = c("#3498db", "#2ecc71", "#e74c3c", "#f1c40f"),
                                     output_file = NULL,
                                     width = 10, 
                                     height = 8) {
  
  # Prepare the data for plotting
  histogram_data <- results_df %>%
    select(all_of(metrics)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "metric",
      values_to = "value"
    ) %>%
    # Convert metric to a factor with specific order and labels
    mutate(metric = factor(metric, levels = metrics, labels = labels))
  
  # Create color palette
  histogram_colors <- setNames(colors, labels)
  
  # Create the plot
  faceted_plot <- ggplot(histogram_data, aes(x = value, fill = metric)) +
    geom_histogram(bins = 50, color = "white", alpha = 0.85) +
    facet_wrap(~ metric, scales = "free", ncol = 2) +
    scale_fill_manual(values = histogram_colors) +
    labs(
      title = "Distribution of DESeq2 Statistics",
      x = NULL,
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  # Save the plot if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, faceted_plot, width = width, height = height, dpi = 300)
  }
}

##################################################
############# Tpm Histogram Plot ##############
##################################################
plot_tpm_distributions <- function(avg_sd_df, 
                                   log_transform = TRUE,
                                   bins = 50,
                                   theme_style = "minimal",
                                   title = "Distribution of TPM Values Across Time Points",
                                   output_file = NULL,
                                   width = 10, 
                                   height = 8) {
  
  # Find columns containing average values
  avg_columns <- grep("avg", names(avg_sd_df), value = TRUE)
  
  # Create a long-format dataframe for plotting
  plot_data <- avg_sd_df %>%
    select(all_of(c("gene_id", "gene_name", avg_columns))) %>%
    pivot_longer(
      cols = all_of(avg_columns),
      names_to = "time_point",
      values_to = "tpm"
    ) %>%
    # Create cleaner labels for time points
    mutate(
      time_point = factor(
        gsub("_avg", "", time_point),
        levels = sort(as.numeric(gsub("_avg", "", avg_columns)))
      )
    )
  
  # Determine x-axis values based on log_transform parameter
  x_values <- if(log_transform) {
    plot_data %>% mutate(plot_value = log2(tpm + 1)) %>% pull(plot_value)
  } else {
    plot_data %>% pull(tpm)
  }
  
  x_label <- if(log_transform) "log2(TPM + 1)" else "TPM"
  
  # Create the plot
  tpm_plot <- ggplot(plot_data, aes(x = if(log_transform) log2(tpm + 1) else tpm)) +
    geom_histogram(bins = bins, fill = "#3498db", color = "white", alpha = 0.8) +
    facet_wrap(~ time_point, ncol = 2) +
    labs(
      title = title,
      x = x_label,
      y = "Number of genes"
    )
  
  # Apply selected theme
  if (theme_style == "paperwhite") {
    tpm_plot <- tpm_plot + 
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold")
      )
  } else {
    tpm_plot <- tpm_plot + get(paste0("theme_", theme_style))()
  }
  
  # Save the plot if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, tpm_plot, width = width, height = height, dpi = 300)
  }
}

##################################################
################## Volcano Plot ##################
##################################################
create_volcano_plot <- function(results_df, 
                                p_value_filter = 0.05,
                                log2_fold_filter = 1,
                                comparisons = NULL,
                                plot_title = "Volcano Plot of Differential Expression",
                                color_scheme = c("red", "blue", "grey"),
                                point_size = 1.0,
                                point_alpha = 0.5,
                                ncol = 2,
                                output_file = NULL,
                                width = 20, 
                                height = 12) {
  
  # Pre-process results for volcano plot
  volcano_data <- results_df %>% 
    mutate(
      negLog10Padj = -log10(padj),
      significance = case_when(
        padj < p_value_filter & log2FoldChange > log2_fold_filter  ~ "Upregulated",
        padj < p_value_filter & log2FoldChange < -log2_fold_filter ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      significance = factor(significance, levels = c("Upregulated", "Downregulated", "Not Significant")),
      gene_label = if_else(significance != "Not Significant", gene_name, "")
    ) %>%
    filter(!is.na(padj) & !is.na(log2FoldChange))
  
  # Filter by comparisons if specified
  if (!is.null(comparisons)) {
    volcano_data <- volcano_data %>%
      filter(comparison %in% comparisons) %>%
      mutate(comparison = factor(comparison, levels = comparisons))
  }
  
  # Define color map
  color_map <- setNames(color_scheme, c("Upregulated", "Downregulated", "Not Significant"))
  
  # Common theme
  common_volcano_theme <- theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.2)),
      strip.text = element_text(face = "bold", size = rel(0.9)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  # Create the plot
  volcano_plot <- ggplot(
    volcano_data,
    aes(x = log2FoldChange, y = negLog10Padj, color = significance, label = gene_label)
  ) +
    geom_point(alpha = point_alpha, size = point_size) +
    scale_color_manual(values = color_map, name = "Significance") +
    geom_vline(xintercept = c(-log2_fold_filter, log2_fold_filter), 
               lty = "dashed", color = "grey40", alpha = 0.6) +
    geom_hline(yintercept = -log10(p_value_filter), 
               lty = "dashed", color = "grey40", alpha = 0.6) +
    facet_wrap(~ comparison, ncol = ncol, scales = "free_y") +
    labs(
      title = plot_title,
      x = expression(Log[2]("Fold Change")),
      y = expression(-Log[10]("Adjusted P-value"))
    ) +
    common_volcano_theme +
    coord_cartesian(ylim = c(0, NA))
  
  # Save the plot if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, volcano_plot, width = width, height = height, dpi = 300)
  }
  
  return(volcano_plot)
}

##################################################
######### LFC Plot - Prolongued Genes ############
##################################################
plot_lfc_prolonged <- function(results, gene_groups, output_file = NULL, 
                               width = 20, height = 12, dpi = 500) {
  
  # Merge results with genes of interest
  res_df <- merge(results, gene_groups)
  
  # Split data into baseline and intermediate comparisons
  df_time0 <- res_df %>% 
    filter(grepl("_vs_0$", comparison)) %>%
    mutate(Time = as.numeric(gsub("time_([0-9]+)_vs_0", "\\1", comparison)))
  
  df_intermediate <- res_df %>% 
    filter(!grepl("_vs_0$", comparison)) %>%
    mutate(
      Time1 = as.numeric(gsub("time_([0-9]+)_vs_([0-9]+)", "\\1", comparison)),
      Time2 = as.numeric(gsub("time_([0-9]+)_vs_([0-9]+)", "\\2", comparison))
    ) %>%
    group_by(gene_name) %>%
    mutate(Offset_Index = row_number() - 1) %>%
    ungroup()
  
  # Calculate y-axis limits for each gene
  y_axis_limits <- res_df %>%
    group_by(gene_name) %>%
    summarize(
      min_log2FC = min(log2FoldChange, na.rm = TRUE),
      max_log2FC = max(log2FoldChange, na.rm = TRUE),
      y_axis_limit = max_log2FC + (max_log2FC - min_log2FC) * 0.05
    )
  
  # Add significance markers
  df_intermediate <- df_intermediate %>%
    left_join(y_axis_limits, by = "gene_name") %>%
    mutate(
      Y_Pos = y_axis_limit + (Offset_Index * 0.3),
      Significance = case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01  ~ "**",
        pvalue < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  # Create the plot
  prolonged_lfc <- ggplot() +
    # Plot log2FoldChange for time 0 comparisons
    geom_point(data = df_time0, 
               aes(x = Time, y = log2FoldChange, color = padj < 0.05), 
               size = 1.5) +
    geom_line(data = df_time0, 
              aes(x = Time, y = log2FoldChange, group = gene_name)) +
    
    # Add p-value markers
    geom_segment(data = df_intermediate %>% filter(Significance != ""), 
                 aes(x = Time1, xend = Time2, y = Y_Pos, yend = Y_Pos), 
                 color = "black", size = 0.5) +
    
    # Add significance markers
    geom_text(data = df_intermediate %>% filter(Significance != ""), 
              aes(x = (Time1 + Time2) / 2, y = Y_Pos + 0.025, label = Significance), 
              size = 4, color = "black") +
    
    # Facet and formatting
    facet_wrap(~gene_name, scales = "free") +
    coord_cartesian() +
    scale_x_continuous(breaks = c(12, 24, 48, 96), 
                       labels = c("12", "24", "48", "96")) +
    labs(title = "Genes regulated by dox over time",
         caption = "*: p < 0.05\n**: p < 0.01\n***: p < 0.001",
         x = "Time after dox induction (hours)",
         y = "log2fold change",
         color = "Significant (padj < 0.05)\ncompared to the zero timepoint") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, prolonged_lfc, width = width, height = height, dpi = dpi)
  }
}

##################################################
########## TPM plot - Prolonged Genes ############
##################################################
plot_tpm_prolonged <- function(avg_and_sd_values, genes_of_interest, g2s, filtered_results,
                               output_file = NULL, width = 20, height = 12, dpi = 500) {
  # Get significant genes
  sig_gene <- data.frame(gene_id = unique(filtered_results$gene_id))
  
  # Merge with TPM values
  sig_gene_tpm <- inner_join(avg_and_sd_values, sig_gene, by = "gene_id") %>%
    merge(g2s)
  
  # Prepare data for genes of interest
  df_plot_genes <- merge(sig_gene_tpm, genes_of_interest)
  
  # Process average values
  time_avg_values <- grep("avg", names(df_plot_genes), value = TRUE)
  df_plot_genes_avg <- df_plot_genes %>%
    pivot_longer(cols = all_of(time_avg_values), names_to = "variable", values_to = "avg") %>%
    mutate(
      Time = as.numeric(gsub("(.*)_avg", "\\1", variable))
    ) %>%
    select(gene_name, Time, avg)
  
  # Process standard deviation values
  time_sd_values <- grep("sd", names(df_plot_genes), value = TRUE)
  df_plot_genes_sd <- df_plot_genes %>%
    pivot_longer(cols = all_of(time_sd_values), names_to = "variable", values_to = "sd") %>%
    mutate(
      Time = as.numeric(gsub("(.*)_sd", "\\1", variable))
    ) %>%
    select(gene_name, Time, sd)
  
  # Combine average and standard deviation
  df_plot_genes_avg_sd <- inner_join(df_plot_genes_avg, df_plot_genes_sd, 
                                     by = c("gene_name", "Time"))
  
  # Create the plot
  prolonged_tpm <- ggplot(df_plot_genes_avg_sd, aes(x = Time, y = avg)) +
    geom_line() +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), width = 3, size = 0.5) +
    labs(x = 'Time after dox induction (hours)', y = 'Average TPM') +
    ggtitle('Average TPM across Time Points') +
    facet_wrap(~ gene_name, scales = "free") +
    scale_x_continuous(breaks = c(0, 12, 24, 48, 96), 
                       labels = c("0", "12", "24", "48", "96")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, prolonged_tpm, width = width, height = height, dpi = dpi)
  }
}

##################################################
################# Heatmap Plot ###################
##################################################
plot_gene_heatmap <- function(data_matrix, 
                              title = "Heatmap of Gene Expression",
                              show_gene_names = FALSE,
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                              scale_method = "row",
                              color_palette = c("red", "white", "blue"),
                              colors_n = 50,
                              distance_method = "euclidean",
                              clustering_method = "ward.D2",
                              break_min = -2,
                              break_max = 2,
                              breaks_n = 51,
                              output_file = NULL,
                              width = 10,
                              height = 12) {
  
  # Define distance function
  custom_dist <- function(x) dist(x, method = distance_method)
  
  # Define clustering function
  custom_hclust <- function(x) hclust(x, method = clustering_method)
  
  # Set breaks for color scale
  breaks <- seq(break_min, break_max, length.out = breaks_n)
  
  # Create color palette
  color_ramp <- colorRampPalette(color_palette)(colors_n)
  
  # Plot the heatmap
  heatmap_plot <- pheatmap(
    data_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    scale = scale_method,
    show_rownames = show_gene_names,
    show_colnames = TRUE,
    main = title,
    distfun = custom_dist,
    color = color_ramp,
    hclustfun = custom_hclust,
    breaks = breaks,
    filename = output_file,
    silent = !is.null(output_file)
  )
  
  return(heatmap_plot)
}

##################################################
################# Circos Plot ####################
##################################################
plot_circos_differential_expression <- function(annotated_results,
                                                output_file = "results/circos_differential_expression.png",
                                                comparisons = NULL,
                                                width = 3000,
                                                height = 3000,
                                                resolution = 300,
                                                color_cap = 2,
                                                point_size = 0.6,
                                                track_height = 0.15,
                                                ylim = c(-20, 20)) {
  # Initialize output file
  png(output_file, width = width, height = height, res = resolution)
  # Clear existing circos plots
  circos.clear()
  # Determine species based on chromosomes
  species <- if(any(grepl("chr22", unique(annotated_results$chr)))) {
    "Human (hg38)"
  } else {
    "Mouse (mm10)"
  }
  # Create a correctly ordered factor for chromosomes based on species
  if(species == "Human (hg38)") {
    sorted_chrs <- factor(
      paste0("chr", c(1:22, "X", "Y")),
      levels = paste0("chr", c(1:22, "X", "Y"))
    )
    species_code <- "hg38"
  } else {
    sorted_chrs <- factor(
      paste0("chr", c(1:19, "X", "Y")),
      levels = paste0("chr", c(1:19, "X", "Y"))
    )
    species_code <- "mm10"
  }
  # Set up circos parameters
  circos.par(start.degree = 90, gap.degree = 2)
  # Initialize with ideogram
  circos.initializeWithIdeogram(
    species = species_code,
    chromosome.index = levels(sorted_chrs)
  )
  # Filter comparisons if specified
  if (!is.null(comparisons)) {
    annotated_results <- annotated_results %>%
      filter(comparison %in% comparisons)
  }
  # Get unique comparisons
  time_points <- unique(annotated_results$comparison)
  # Find global max log2FC value with margin
  global_max_fc <- max(abs(annotated_results$log2FoldChange), na.rm = TRUE) * 1.2
  # Define color mapping function
  get_color <- function(value) {
    if (is.na(value)) return("gray80")
    capped_value <- max(min(value, color_cap), -color_cap)
    intensity <- abs(capped_value) / color_cap
    intensity <- sqrt(intensity) * 0.8
    if (capped_value < 0) {
      return(rgb(1 - intensity * 0.7, 1 - intensity * 0.7, 1))  # Blue
    } else if (capped_value > 0) {
      return(rgb(1, 1 - intensity * 0.7, 1 - intensity * 0.7))  # Red
    } else {
      return("white")  # Zero value
    }
  }
  # Function to plot each comparison
  plot_timepoint_comparison <- function(comparison_data, track_title) {
    if (nrow(comparison_data) > 0) {
      circos.genomicTrack(
        comparison_data,
        ylim = ylim,
        track.height = track_height,
        bg.border = "black",
        panel.fun = function(region, value, ...) {
          # For each genomic region in this cell/sector
          for (i in 1:nrow(region)) {
            log2fc <- value$log2FoldChange[i]
            point_color <- get_color(log2fc)
            
            circos.points(
              mean(c(region$start[i], region$end[i])),
              log2fc,
              col = point_color,
              pch = 16,
              cex = point_size
            )
          }
          # Add horizontal line at y=0
          circos.lines(CELL_META$xlim, c(0, 0), col = "black", lty = 2)
          # Add track title on the first sector
          if (CELL_META$sector.index == levels(sorted_chrs)[1]) {
            circos.text(
              CELL_META$xcenter, 
              max(ylim) * 0.8,
              track_title,
              facing = "inside", 
              niceFacing = TRUE,
              adj = c(0.5, 0),
              cex = 0.8
            )
          }
        }
      )
    } else {
      message(paste("No data for comparison:", track_title))
    }
  }
  # Plot each timepoint comparison
  for (time_point in time_points) {
    comparison_data <- annotated_results %>%
      filter(comparison == time_point) %>%
      dplyr::select(chr, start, end, log2FoldChange)
    plot_timepoint_comparison(comparison_data, time_point)
  }
  # Add title and legend
  title(main = paste0("Differential Expression Across Time Points (", species, ")"), cex.main = 1.2)
  par(xpd = TRUE)
  # Add color legend
  legend_values <- seq(-color_cap, color_cap, length.out = 7)
  legend_x <- 0.9
  legend_y <- -1.3
  color_legend_width <- 0.2
  color_legend_height <- 0.02
  # Draw color bar segments
  for (i in 1:(length(legend_values)-1)) {
    mid_val <- (legend_values[i] + legend_values[i+1]) / 2
    segment_color <- get_color(mid_val)
    rect(
      legend_x + (i-1) * color_legend_width/(length(legend_values)-1),
      legend_y - color_legend_height,
      legend_x + i * color_legend_width/(length(legend_values)-1),
      legend_y,
      col = segment_color,
      border = NA
    )
  }
  # Add tick marks and labels
  for (i in 1:length(legend_values)) {
    text(
      legend_x + (i-1) * color_legend_width/(length(legend_values)-1),
      legend_y - color_legend_height*1.5,
      labels = sprintf("%.1f", legend_values[i]),
      cex = 0.7,
      adj = c(0.5, 0.5)
    )
  }
  # Add legend title
  text(
    legend_x + color_legend_width/2,
    legend_y + color_legend_height*1.5,
    "log2FoldChange",
    cex = 0.8
  )
  # Close device
  dev.off()
}

##################################################
########## TPM Plot with Peak Group ##############
##################################################
plot_tpm_atac <- function(avg_tpm_up_prolonged_genes,
                          avg_tpm_up_prolonged_genes_nopeak,
                          title = "Mean TPM ± SEM at Each Timepoint",
                          output_file = NULL, width = 10, height = 6, 
                          dpi = 500) {
# Number of genes
n_peak <- nrow(avg_tpm_up_prolonged_genes)
n_nopeak <- nrow(avg_tpm_up_prolonged_genes_nopeak)

# Means
means_peak <- colMeans(avg_tpm_up_prolonged_genes, na.rm = TRUE)
means_nopeak <- colMeans(avg_tpm_up_prolonged_genes_nopeak, na.rm = TRUE)

# Standard deviations
sd_peak <- apply(avg_tpm_up_prolonged_genes, 2, sd, na.rm = TRUE)
sd_nopeak <- apply(avg_tpm_up_prolonged_genes_nopeak, 2, sd, na.rm = TRUE)

# Calculate SEM
sem_peak <- sd_peak / sqrt(n_peak)
sem_nopeak <- sd_nopeak / sqrt(n_nopeak)

# Make combined data frame
df_combined <- data.frame(
  timepoint = rep(names(means_peak), 2),
  mean_tpm = c(means_peak, means_nopeak),
  sem_tpm = c(sem_peak, sem_nopeak),
  group = rep(c(
    paste0("Peak (n = ", n_peak, ")"),
    paste0("No Peak (n = ", n_nopeak, ")")
  ), each = length(means_peak))
)

prolonged_tpm <- ggplot(df_combined, aes(x = timepoint, y = mean_tpm, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean_tpm - sem_tpm, ymax = mean_tpm + sem_tpm),
                width = 0.2, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  labs(title = title,
       x = "Timepoint",
       y = "Mean TPM",
       fill = "Group") +
  scale_fill_manual(values = c("grey", "orange"))

# Save if output file specified
if (!is.null(output_file)) {
  ggsave(output_file, prolonged_tpm, width = width, height = height, dpi = dpi)}
}
