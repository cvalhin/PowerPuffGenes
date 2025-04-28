PowerPuffGenes Report
================
2025-04-28

# SECTION 1: RNAseq analysis

Understanding the temporal dynamics of gene expression is crucial for
deciphering complex biological processes like stem cell differentiation
or response to stimuli. This project utilizes a doxycycline
(dox)-inducible system in stem cells, where dox triggers gene expression
via the rtTA-TRE mechanism. While powerful, concerns exist that dox
itself might exert off-target effects on the transcriptome. Leveraging
Dr. Rinn’s lab RNA-seq dataset spanning five time points (0, 12, 24, 48,
96 hours) with replicates, this analysis aims to:

1)  Identify genes significantly regulated following dox administration.
2)  Characterize the temporal profiles of these regulated genes,
    distinguishing between transient and sustained responses.
3)  Investigate the role, significance, and chromosomal location of
    these genes.

We employ the nf-core/rnaseq pipeline for data processing and DESeq2 for
differential expression analysis, providing a practical application of
bioinformatics skills to time-course data exploration.

## 1.1) Data preprocessing

To filter genes that significantly change from the zero time point (no
dox)

### Parameters and paths definition

``` r
# Set species parameter - set to either "mouse" or "human"
SPECIES          <- "human"  # Change to "human" for human data

# Analysis parameters
FILTER_TRESHOLD  <- 0
P_VALUE_FILTER   <- 0.01
LOG2_FOLD_FILTER <- 1
```

### Loading and preproccesing count

``` r
# Read raw data and make dictionary mapping gene id with gene name
counts_orig  <- read.table(COUNTS_PATH, header=TRUE, row.names=1)
g2s          <- data.frame(gene_id = rownames(counts_orig), 
                           gene_name = counts_orig[ , 1]) 

# Remove gene name column from counts and convert counts to numerical rounded matrix
counts           <- counts_orig %>% select(-gene_name)

# Filter genes with no counts for any of the samples
counts_filtered <- counts %>% 
  as.matrix() %>% 
  round() %>% 
  .[rowSums(.) > FILTER_TRESHOLD, ]
```

### DESeq2 analysis for counts

``` r
# Make a column from the titles of the columns of the counts_matrix table
deseq_samples  <- data.frame(sample_id = colnames(counts))

# Get the time and replicate values from the sample names in deseq_samples
split_values     <- strsplit(deseq_samples$sample_id, "_")
time_values      <- sapply(split_values, function(x) x[[2]])
replicate_values <- sapply(split_values, function(x) x[[3]])

# Add time and replicate values as columns to deseq_samples and factor them
deseq_samples$time_point <- time_values
deseq_samples$replicate  <- replicate_values
deseq_samples$time_point <- factor(deseq_samples$time_point, levels=c("0","12","24","48","96"))
deseq_samples$replicate  <- factor(deseq_samples$replicate, levels =c("1","2","3"))
deseq_samples <- drop_na(deseq_samples)

# Remove extra timepoints and replicates (counts)
keep_cols <- deseq_samples$sample_id
counts_filtered  <- counts_filtered[, colnames(counts_filtered) %in% keep_cols, drop = FALSE]

# Testing whether sample sheet and counts are arranged properly 
stopifnot(all(colnames(counts) == rownames(deseq_samples$sample_id)))

# Prepare DESeq dataset and run DESea2
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = deseq_samples,
                              design = ~ time_point) 
dds <- DESeq(dds)
```

``` r
pca_plot <- plot_pca(dds, intgroup=c("time_point", "replicate"))
ggsave(file.path(RESULTS_PATH, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)
```

### Obtaining and filtering count results

``` r
# Define the comparisons
comparison_list <- list(
  # Time points against baseline
  time_12_vs_0 = c("time_point", "12", "0"),
  time_24_vs_0 = c("time_point", "24", "0"),
  time_48_vs_0 = c("time_point", "48", "0"),
  time_96_vs_0 = c("time_point", "96", "0"),
  # Time series comparisons
  time_24_vs_12 = c("time_point", "24", "12"),
  time_48_vs_24 = c("time_point", "48", "24"),
  time_96_vs_48 = c("time_point", "96", "48")
)

# Create empty df to store results values
results <- tibble(
  gene_id = character(0),
  baseMean = numeric(0),
  log2FoldChange = numeric(0),
  lfcSE = numeric(0),
  stat = numeric(0),
  pvalue = numeric(0),
  padj = numeric(0),
  gene_name = character(0),
  comparison = character(0)
)

# Process all comparisons
for (comp_name in names(comparison_list)) {
  comp_result <- process_comparison(
    dds_obj = dds,
    comparison_params = comparison_list[[comp_name]],
    comparison_name = comp_name,
    gene_map = g2s
  )
  if (!is.null(comp_result)) {
    results <- bind_rows(results, comp_result)
  }
}

results_unique <- unique(results$gene_name)
results_unique <- as.data.frame(results_unique)
print(c("This are the total genes for all our comparisons before filtering", length(results_unique$results_unique)))
```

    ## [1] "This are the total genes for all our comparisons before filtering"
    ## [2] "42871"

``` r
# Filter based on p-value and log2fold changes defined at top
filtered_results <- results %>%
  filter(padj < P_VALUE_FILTER, abs(log2FoldChange) > LOG2_FOLD_FILTER)

print(c("This are the total genes our filtered results with the se pvalue and log2 fold filter",  length(filtered_results$gene_name)))
```

    ## [1] "This are the total genes our filtered results with the se pvalue and log2 fold filter"
    ## [2] "3092"

``` r
# Get unique gene names that are significant
filtered_genes <- filtered_results %>%
  select(gene_name) %>%
  distinct()

save(results, filtered_results, filtered_genes, 
     file=file.path(RESULTS_PATH, "results_counts.RData"))
```

### Obtaining and filtering tmp results

``` r
# Read raw tpm data and delete gene name column
tpm_orig  <- read.table(TPM_PATH, header=TRUE, row.names=1)
tpm       <- tpm_orig %>% select(-gene_name)

# Filter genes with no counts for any of the samples
tpm_filtered <- tpm[rowSums(tpm) > FILTER_TRESHOLD, ]

# Remove extra timepoints and replicates (tpm)
tpm_filtered <- tpm_filtered[, colnames(tpm_filtered) %in% keep_cols, drop = FALSE]

# Calculate average and standard deviation for replicates at each time point
time_points <- c("0", "12", "24", "48", "96")

# Use map to process each time point and bind the results
avg_and_sd_values <- time_points %>%
  map_dfc(function(tp) {
    # Find columns for this time point - adapt for different column naming patterns
    cols <- grep(paste0("^", SAMPLE_PREFIX, "_", tp, "_"), colnames(tpm_filtered), value = TRUE)
    
    # Calculate mean and sd across replicates
    data.frame(
      avg = rowMeans(tpm_filtered[, cols]),
      sd = apply(tpm_filtered[, cols], 1, sd)
    ) %>% 
    # Rename columns to include time point
    setNames(paste0(tp, c("_avg", "_sd")))
  }) %>% 
  # Add gene identifiers and merge with gene symbols
  rownames_to_column("gene_id") %>% 
  merge(g2s)

save(tpm_filtered, avg_and_sd_values, file = file.path(RESULTS_PATH, "results_tpm.RData"))
```

### Distribution of filtered results

``` r
plot_count_distributions(
  filtered_results,
  output_file = file.path(RESULTS_PATH, "counts_distributions.png")
)

# Use the function with your avg_and_sd_values
plot_tpm_distributions(
  avg_sd_df = avg_and_sd_values,
  output_file = file.path(RESULTS_PATH, "tpm_distributions.png")
)
```

### Volcano plots

``` r
# Baseline comparisons
baseline_comparisons <- c("time_12_vs_0", "time_24_vs_0", "time_48_vs_0", "time_96_vs_0")
plot_baseline <- create_volcano_plot(
  results,
  p_value_filter = P_VALUE_FILTER,
  log2_fold_filter = LOG2_FOLD_FILTER,
  comparisons = baseline_comparisons,
  plot_title = paste0("Volcano Plots (", toupper(SPECIES), "): Comparisons vs Baseline"),
  output_file = file.path(RESULTS_PATH, "volcano_baseline.png")
)

# Time series comparisons
timeseries_comparisons <- c("time_12_vs_0", "time_24_vs_12", "time_48_vs_24", "time_96_vs_48")
plot_timeseries <- create_volcano_plot(
  results,
  p_value_filter = P_VALUE_FILTER,
  log2_fold_filter = LOG2_FOLD_FILTER,
  comparisons = timeseries_comparisons,
  plot_title = paste0("Volcano Plots (", toupper(SPECIES), "): Sequential Time Comparisons"),
  output_file = file.path(RESULTS_PATH, "volcano_timeseries.png")
)
```

### Circos Plot

``` r
# Define chromosomes based on species
valid_chromosomes <- if(SPECIES == "mouse") {
  paste0("chr", c(1:19, "X", "Y"))
} else if(SPECIES == "human") {
  paste0("chr", c(1:22, "X", "Y"))
} else {
  stop("Invalid species: ", SPECIES)
}

# Prepare gene annotations
gene_locations <- prepare_gene_annotations(
  gtf_file = GTF_PATH,
  valid_chromosomes = valid_chromosomes
)

# Merge DE results with annotations
annotated_results <- merge_results_with_annotations(
  de_results = filtered_results,
  gene_locations = gene_locations
)

# Generate circos plot
# Using only baseline comparisons (time vs 0)
baseline_comparisons <- grep("_vs_0$", unique(annotated_results$comparison), value = TRUE)

plot_circos_differential_expression(
  annotated_results = annotated_results,
  output_file = file.path(RESULTS_PATH, "circos.png"),
  comparisons = baseline_comparisons,
  color_cap = 2,       # Cap color intensity at log2FC = 2
  point_size = 0.6,    # Size of points
  track_height = 0.15  # Height of each track
)
```

    ## png 
    ##   2

## 1.2) Gene selection

To identify and visualize prolonged and transient expression patterns
from significant gene expression dynamics after dox induction

### Classifying genes by their differential expression dynamics

Find genes that hold differential expression up to the end time point
(96 hours)

``` r
# Extract unique gene IDs from the filtered results and create a data frame
sig_gene <- as.data.frame(unique(filtered_results$gene_id))  # Get unique gene IDs
names(sig_gene) <- "gene_id"  # Rename the column to "gene_id"

# Subset the data frame containing TPM values
tpm_filtered$gene_id <- rownames(tpm_filtered)
rownames(tpm_filtered) <- NULL
sig_gene_tpm_all <- inner_join(tpm_filtered, sig_gene, by = "gene_id")
sig_gene_tpm_all <- merge(sig_gene_tpm_all, g2s)  # Merge with gene-to-symbol mapping (g2s)

# Categorize genes into expression groups
gene_groups <- categorize_gene_expression(filtered_results)
gene_groups_df <- data.frame(
  gene_name = c(gene_groups$change_at_12_to_end, 
               gene_groups$change_at_24_to_end, 
               gene_groups$change_at_48_to_end),
  change_from = c(rep("12", length(gene_groups$change_at_12_to_end)),
                 rep("24", length(gene_groups$change_at_24_to_end)),
                 rep("48", length(gene_groups$change_at_48_to_end)))
)

# Prepare TPM data for plotting genes of interest by merging with genes of interest
selected_tpm <- merge(sig_gene_tpm_all, gene_groups_df)

# Dynamically generate column names for 0h and 96h triplicates
t0_cols <- paste0(SAMPLE_PREFIX, "_0_", 1:3)
t96_cols <- paste0(SAMPLE_PREFIX, "_96_", 1:3)

# Apply function to each row/gene
stat_df <- do.call(rbind, apply(selected_tpm, 1, calculate_stats))
selected_tpm <- cbind(selected_tpm, stat_df)

# Adjust p-values for multiple testing (Benjamini-Hochberg)
selected_tpm$adjusted_pvalue <- p.adjust(selected_tpm$pvalue, method = "BH")

# Get significant results
significant_genes <- selected_tpm[selected_tpm$adjusted_pvalue < P_VALUE_FILTER, ]
genes_of_interest <- significant_genes %>% select(gene_name, change_from, tpm_log2fc)
genes_of_interest <- genes_of_interest[order(genes_of_interest$change_from),]

# Get transient genes
transient_genes <- filtered_genes[!(filtered_genes$gene_name %in%
                                      genes_of_interest$gene_name),]

# Save to csv file
write_csv(genes_of_interest, file.path(RESULTS_PATH, "prolonged_gene_list.csv"))

print(paste("Number of genes with expression sustained until the end:",
            length(genes_of_interest$gene_name)))
```

    ## [1] "Number of genes with expression sustained until the end: 134"

``` r
print(paste("Number of genes exhibiting transient expression:",
            length(transient_genes$gene_name)))
```

    ## [1] "Number of genes exhibiting transient expression: 1494"

### Heatmaps

``` r
# Prepare data for the heatmap of all significant genes
rownames(tpm_filtered) <- tpm_filtered$gene_id
tpm_filtered$gene_id <- NULL
all_sig_genes_data <- prepare_heatmap_data(
  tpm_matrix = tpm_filtered,
  gene_symbol_map = g2s,
  gene_filter = filtered_genes$gene_name
)

# Generate heatmap for all significant genes
all_genes_heatmap <- plot_gene_heatmap(
  data_matrix = all_sig_genes_data$numeric_matrix,
  title = paste0("Heatmap of Log2-Transformed TPM Values for All Significant Genes (", 
                toupper(SPECIES), ")"),
  show_gene_names = FALSE,
  output_file = file.path(RESULTS_PATH, "heatmap_sig_genes.png")
)
# Prepare data for the heatmap of prolonged genes
prolonged_genes_data <- prepare_heatmap_data(
  tpm_matrix = tpm_filtered,
  gene_symbol_map = g2s,
  gene_filter = genes_of_interest$gene_name
)

# Generate heatmap for prolonged genes
prolonged_genes_heatmap <- plot_gene_heatmap(
  data_matrix = prolonged_genes_data$numeric_matrix,
  title = paste0("Heatmap of Log2-Transformed TPM Values of Prolonged Genes (", 
                toupper(SPECIES), ")"),
  show_gene_names = FALSE,
  output_file = file.path(RESULTS_PATH, "heatmap_prolonged_genes.png")
)
```

## SECTION 2: ATACseq analysis

To determine if chromatin accessibility changes due to dox exposure in
stem cells. We have performed a time course series of experiments
measuring chromatin accessibility (ATACseq peaks) upon exposure to dox.
These are 0, 30, 60, 90, 120, 150 minutes. Each time point has a
replicate. The fastq sequencing files were processed by the NF_CORE
ATACseq pipeline v-2.12 (<https://nf-co.re/atacseq/2.1.2/>).

### Loading in ATACseq peak files with custom function import_peaks (list of GRanges output)

``` r
# creating a file list also needed for import_peaks function to get sample name associated with file
fl <- list.files(PEAK_PATH, full.names = TRUE, pattern = ".broadPeak")

# running import_peaks
my_peaks <- import_peaks(PEAK_PATH)

# finding consensus peaks among replicates
my_consensus_peaks <- find_consensus_peaks(my_peaks)

# run find_common_peaks function
common_peaks <-  suppressWarnings(find_common_peaks(my_consensus_peaks))

print("here are the number of peaks for each sample")
```

    ## [1] "here are the number of peaks for each sample"

``` r
print(num_consensus_peaks <- sapply(my_consensus_peaks, length) %>% as.data.frame)
```

    ##                      .
    ## novageneGFP_0   101526
    ## novageneGFP_30   96730
    ## novageneGFP_60   69841
    ## novageneGFP_90   74112
    ## novageneGFP_120  99919
    ## novageneGFP_150 104256

### Peaks that are unique to dox and non-dox conditions

Using find_common_overlaps peaks that are specific to dox or non-dox
will be identified.

Non-dox samples

``` r
# common peaks in non-dox (0 time point)
non_dox_samples <- names(my_consensus_peaks)[c(grep("_0", names(my_consensus_peaks)))]
non_dox_peaks <- my_consensus_peaks[non_dox_samples]
if (length(non_dox_samples) > 1) {
  non_dox_common_peaks <- suppressWarnings(find_common_peaks(non_dox_peaks))
} else {
  non_dox_common_peaks <- non_dox_peaks[[1]]
  }
print(c("This is how many peaks are common in non-dox:", length(non_dox_common_peaks)))
```

    ## [1] "This is how many peaks are common in non-dox:"
    ## [2] "101526"

Dox samples

``` r
# common peaks in dox time points (!not time 0)
dox_samples <- names(my_consensus_peaks)[!grepl("_0$", names(my_consensus_peaks))]
dox_peaks <- my_consensus_peaks[dox_samples]
dox_common_peaks <- suppressWarnings(find_common_peaks(dox_peaks))

print(c("This is how many peaks are common in dox:", length(dox_common_peaks)))
```

    ## [1] "This is how many peaks are common in dox:"
    ## [2] "58378"

Overlap of dox and non-dox (to get unique to each condition)

``` r
# Now overlap between dox and non-dox common peaks
dox_compare_list <- list(non_dox = non_dox_common_peaks, dox = dox_common_peaks)
dox_non_dox_ov <- suppressWarnings(find_common_peaks(dox_compare_list))
print(c("This is how many peaks are common in both non- and dox treatments", length(dox_non_dox_ov)))
```

    ## [1] "This is how many peaks are common in both non- and dox treatments"
    ## [2] "58317"

``` r
# extracting peaks unique to each condition (dox non-dox)
# Peaks unique to non_dox
unique_to_non_dox <- suppressWarnings(find_my_peaks(dox_non_dox_ov, non_dox_common_peaks))
print(c("This is how many peaks are unique to non-dox condition", length(unique_to_non_dox)))
```

    ## [1] "This is how many peaks are unique to non-dox condition"
    ## [2] "43209"

``` r
# Peaks unique to dox
unique_to_dox <- suppressWarnings(find_my_peaks(dox_non_dox_ov, dox_common_peaks))
print(c("This is how many peaks are unique to dox condition", length(unique_to_dox)))
```

    ## [1] "This is how many peaks are unique to dox condition"
    ## [2] "328"

### Creating gene, lincrna, mRNA annotation GRange objects

``` r
# Loading gencode genome annotation as GRanges 
gencode_gr <- rtracklayer::import(GTF_PATH)

# all genes
gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
gene_promoters <- promoters(gencode_genes, upstream = 2000, downstream = 2000)
```

### Comparing overlaps of dox with gene annotaitons

Now we will overlap our dox and non-dox unique peaks with genome
annotations (gene promoters) First we will find number of overlaps with
gene promoters and then genes that had changed in RNAseq

``` r
# peaks common to dox condition overlapped with gene promoters
gr_list_gene_promoter_dox_ov <- list(gene_promoters = gene_promoters, dox_peaks = dox_common_peaks)
dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_dox_ov))

print(c("This is how many dox common peaks overlapped gene promoters", length(dox_gene_promoter_ov)))
```

    ## [1] "This is how many dox common peaks overlapped gene promoters"
    ## [2] "22896"

``` r
gr_list_gene_promoter_dox_ov <-  list(gene_promoters = gene_promoters, dox_peaks = unique_to_dox)
unique_to_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_dox_ov))

print(c("This is how many dox unique peaks overlapped gene promoters", length(unique_to_dox_gene_promoter_ov)))
```

    ## [1] "This is how many dox unique peaks overlapped gene promoters"
    ## [2] "40"

``` r
# peaks common to non-dox condition overlapped with gene promoters
gr_list_gene_promoter_dox_ov <- list(gene_promoters = gene_promoters, dox_peaks = non_dox_common_peaks)
non_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_dox_ov))

print(c("This is how many non_dox common peaks overlapped gene promoters", length(non_dox_gene_promoter_ov)))
```

    ## [1] "This is how many non_dox common peaks overlapped gene promoters"
    ## [2] "26274"

``` r
gene_promoter_lostindox <- setdiff(non_dox_gene_promoter_ov$gene_name,
                                   dox_gene_promoter_ov$gene_name)
print(c("This is how many dox lost peaks overlapped gene promoters", length(gene_promoter_lostindox)))
```

    ## [1] "This is how many dox lost peaks overlapped gene promoters"
    ## [2] "3425"

### Matching ATAC peaks to prolonged genes

``` r
# filter upregulated-prolonged genes
up_genes <- genes_of_interest %>%
  filter(tpm_log2fc > 0)
print(c("this is how many prolonged genes that are upregulated with dox", length(up_genes$gene_name)))
```

    ## [1] "this is how many prolonged genes that are upregulated with dox"
    ## [2] "54"

``` r
# get TPM for all prolonged genes that are upregulated
tpm_up_genes <- avg_and_sd_values[avg_and_sd_values$gene_name %in% 
                                    up_genes$gene_name, ]
avg_tpm_up_genes <- tpm_up_genes[, grep("_avg", colnames(tpm_up_genes))]
print("TPM summary of prolonged genes that are upregulated with dox")
```

    ## [1] "TPM summary of prolonged genes that are upregulated with dox"

``` r
print(summary(avg_tpm_up_genes))
```

    ##      0_avg             12_avg              24_avg              48_avg         
    ##  Min.   : 0.0000   Min.   :  0.06089   Min.   :  0.05253   Min.   :  0.06119  
    ##  1st Qu.: 0.1932   1st Qu.:  0.59910   1st Qu.:  0.85075   1st Qu.:  1.10122  
    ##  Median : 0.7968   Median :  1.44334   Median :  2.56292   Median :  2.96418  
    ##  Mean   : 3.2043   Mean   :  6.81834   Mean   :  9.26534   Mean   : 13.02380  
    ##  3rd Qu.: 1.9194   3rd Qu.:  3.00736   3rd Qu.:  4.19604   3rd Qu.:  6.56835  
    ##  Max.   :75.4711   Max.   :100.04414   Max.   :138.47797   Max.   :227.25070  
    ##      96_avg        
    ##  Min.   :  0.0566  
    ##  1st Qu.:  1.0854  
    ##  Median :  3.2383  
    ##  Mean   : 15.3117  
    ##  3rd Qu.:  5.9180  
    ##  Max.   :398.2863

``` r
unique_gene_promoters <- intersect(unique_to_dox_gene_promoter_ov$gene_name, up_genes$gene_name)
print(c("this is how many upregulated-prolonged genes overlapped with unique dox gene promoters", length(unique_gene_promoters)))
```

    ## [1] "this is how many upregulated-prolonged genes overlapped with unique dox gene promoters"
    ## [2] "0"

``` r
if (length(unique_gene_promoters) > 0) {
  print(unique_gene_promoters)
}

# filter upregulated-prolonged genes to dox promoter overlaps
up_prolonged_rnaseq_atac_dox <- dox_gene_promoter_ov[dox_gene_promoter_ov$gene_name %in% up_genes$gene_name]
up_atac_genes <- up_genes %>%
  filter(gene_name %in% up_prolonged_rnaseq_atac_dox$gene_name)

print(c("this is how many upregulated genes overlap between prolonged expression and ATACseq dox peaks", length(up_prolonged_rnaseq_atac_dox)))
```

    ## [1] "this is how many upregulated genes overlap between prolonged expression and ATACseq dox peaks"
    ## [2] "37"

``` r
# get TPM for prolonged genes that are upregulated with ATAC peaks
tpm_up_prolonged_genes <- avg_and_sd_values[avg_and_sd_values$gene_name %in%
                                              up_atac_genes$gene_name, ]
avg_tpm_up_prolonged_genes <- tpm_up_prolonged_genes[, grep("_avg",
                                                            colnames(tpm_up_prolonged_genes))]
print("TPM summary of prolonged genes that are upregulated with ATAC peaks in dox")
```

    ## [1] "TPM summary of prolonged genes that are upregulated with ATAC peaks in dox"

``` r
print(summary(avg_tpm_up_prolonged_genes))
```

    ##      0_avg              12_avg              24_avg              48_avg        
    ##  Min.   : 0.02471   Min.   :  0.06753   Min.   :  0.05615   Min.   :  0.3217  
    ##  1st Qu.: 0.36194   1st Qu.:  0.62470   1st Qu.:  0.86018   1st Qu.:  1.6252  
    ##  Median : 1.37009   Median :  1.86908   Median :  3.16803   Median :  4.6035  
    ##  Mean   : 4.48831   Mean   :  8.48570   Mean   : 11.45683   Mean   : 16.9570  
    ##  3rd Qu.: 2.37833   3rd Qu.:  3.42231   3rd Qu.:  4.92961   3rd Qu.:  6.8530  
    ##  Max.   :75.47109   Max.   :100.04414   Max.   :138.47797   Max.   :227.2507  
    ##      96_avg        
    ##  Min.   :  0.2395  
    ##  1st Qu.:  1.6551  
    ##  Median :  4.0641  
    ##  Mean   : 20.3072  
    ##  3rd Qu.:  6.7870  
    ##  Max.   :398.2863

``` r
# get TPM for prolonged genes that are upregulated without ATAC peaks
up_no_atac_genes <- setdiff(up_genes$gene_name, up_atac_genes$gene_name)
tpm_up_prolonged_genes_nopeak <- avg_and_sd_values[avg_and_sd_values$gene_name %in%
                                              up_no_atac_genes, ]
avg_tpm_up_prolonged_genes_nopeak <- tpm_up_prolonged_genes_nopeak[, grep("_avg",
                                                        colnames(tpm_up_prolonged_genes_nopeak))]
print("TPM summary of prolonged genes that are upregulated without ATAC peaks in dox")
```

    ## [1] "TPM summary of prolonged genes that are upregulated without ATAC peaks in dox"

``` r
print(summary(avg_tpm_up_prolonged_genes_nopeak))
```

    ##      0_avg             12_avg             24_avg             48_avg        
    ##  Min.   :0.00000   Min.   : 0.06089   Min.   : 0.05253   Min.   : 0.06119  
    ##  1st Qu.:0.08112   1st Qu.: 0.40149   1st Qu.: 0.50739   1st Qu.: 0.62973  
    ##  Median :0.15826   Median : 0.74447   Median : 1.22931   Median : 1.64629  
    ##  Mean   :0.40954   Mean   : 3.18939   Mean   : 4.49565   Mean   : 4.46338  
    ##  3rd Qu.:0.61811   3rd Qu.: 2.04511   3rd Qu.: 2.83060   3rd Qu.: 2.74017  
    ##  Max.   :1.93732   Max.   :24.23190   Max.   :33.58028   Max.   :30.67189  
    ##      96_avg       
    ##  Min.   : 0.0566  
    ##  1st Qu.: 0.5530  
    ##  Median : 1.1374  
    ##  Mean   : 4.4393  
    ##  3rd Qu.: 2.9951  
    ##  Max.   :31.3043

``` r
# plot avg TPM for upregulated-prolonged genes with peaks and no peak
plot_tpm_atac(avg_tpm_up_prolonged_genes, avg_tpm_up_prolonged_genes_nopeak,
              title="Mean TPM ± SEM of Upregulated-Prolonged Genes",
              output_file=file.path(RESULTS_PATH, "bar_tpm_prolonged_up_genes.png"))

# filter downregulated-prolonged genes
down_genes <- genes_of_interest %>%
  filter(tpm_log2fc < 0)
print(c("this is how many prolonged genes that are downregulated with dox", length(down_genes$gene_name)))
```

    ## [1] "this is how many prolonged genes that are downregulated with dox"
    ## [2] "80"

``` r
# get TPM for all prolonged genes that are downregulated
tpm_down_genes <- avg_and_sd_values[avg_and_sd_values$gene_name %in% 
                                    down_genes$gene_name, ]
avg_tpm_down_genes <- tpm_down_genes[, grep("_avg", colnames(tpm_down_genes))]
print("TPM summary of prolonged genes that are downregulated with dox")
```

    ## [1] "TPM summary of prolonged genes that are downregulated with dox"

``` r
print(summary(avg_tpm_down_genes))
```

    ##      0_avg              12_avg              24_avg              48_avg        
    ##  Min.   :  0.0875   Min.   :  0.04292   Min.   :  0.01916   Min.   : 0.01915  
    ##  1st Qu.:  1.7509   1st Qu.:  1.08741   1st Qu.:  0.79561   1st Qu.: 0.60420  
    ##  Median :  4.8924   Median :  2.75020   Median :  2.04085   Median : 1.65176  
    ##  Mean   : 20.5582   Mean   : 12.41719   Mean   :  9.62737   Mean   : 6.64077  
    ##  3rd Qu.: 12.4985   3rd Qu.:  7.51060   3rd Qu.:  6.35270   3rd Qu.: 6.01086  
    ##  Max.   :308.0928   Max.   :188.11414   Max.   :132.48833   Max.   :94.47125  
    ##      96_avg       
    ##  Min.   : 0.0125  
    ##  1st Qu.: 0.5536  
    ##  Median : 1.6455  
    ##  Mean   : 5.9127  
    ##  3rd Qu.: 5.1778  
    ##  Max.   :92.8702

``` r
lost_gene_promoters <- intersect(gene_promoter_lostindox, down_genes$gene_name)
print(c("this is how many downregulated-prolonged genes overlapped with the lost gene promoters in dox", length(lost_gene_promoters)))
```

    ## [1] "this is how many downregulated-prolonged genes overlapped with the lost gene promoters in dox"
    ## [2] "8"

``` r
if (length(lost_gene_promoters) > 0) {
  print(lost_gene_promoters)
}
```

    ## [1] "MIR34AHG"    "TRDN"        "RP11-69I8.2" "LMO3"        "MAG"        
    ## [6] "NLRP12"      "SERTM2"      "SRY"

``` r
# filter downregulated-prolonged genes to non dox promoter overlaps
down_prolonged_rnaseq_atac_dox <- dox_gene_promoter_ov[dox_gene_promoter_ov$gene_name %in% down_genes$gene_name]
down_no_atac_genes <- down_genes %>%
  filter(!(gene_name %in% down_prolonged_rnaseq_atac_dox$gene_name))

print(c("this is how many downregulated genes that couldn't find ATACseq peaks in dox", length(down_prolonged_rnaseq_atac_dox)))
```

    ## [1] "this is how many downregulated genes that couldn't find ATACseq peaks in dox"
    ## [2] "64"

``` r
# get TPM for all prolonged genes that are downregulated without peak
tpm_down_prolonged_genes_nopeak <- avg_and_sd_values[avg_and_sd_values$gene_name %in%
                                                     down_no_atac_genes$gene_name, ]
avg_tpm_down_prolonged_genes_nopeak <- tpm_down_prolonged_genes_nopeak[, grep("_avg",
                                                        colnames(tpm_down_prolonged_genes_nopeak))]
print("TPM summary of prolonged genes that are downregulated with no peak in dox")
```

    ## [1] "TPM summary of prolonged genes that are downregulated with no peak in dox"

``` r
print(summary(avg_tpm_down_prolonged_genes_nopeak))
```

    ##      0_avg             12_avg             24_avg             48_avg        
    ##  Min.   : 0.1564   Min.   : 0.04486   Min.   : 0.05714   Min.   : 0.03374  
    ##  1st Qu.: 0.9887   1st Qu.: 0.41813   1st Qu.: 0.36252   1st Qu.: 0.24423  
    ##  Median : 2.3017   Median : 0.95486   Median : 0.85451   Median : 0.47488  
    ##  Mean   : 7.8409   Mean   : 3.43683   Mean   : 3.09011   Mean   : 2.44647  
    ##  3rd Qu.:10.3830   3rd Qu.: 5.23241   3rd Qu.: 4.45491   3rd Qu.: 4.04695  
    ##  Max.   :39.8790   Max.   :14.55741   Max.   :14.76529   Max.   :11.19726  
    ##      96_avg       
    ##  Min.   :0.04314  
    ##  1st Qu.:0.24698  
    ##  Median :0.52429  
    ##  Mean   :2.42544  
    ##  3rd Qu.:4.24848  
    ##  Max.   :8.90358

``` r
# get TPM for all prolonged genes that are downregulated with peak
down_atac_genes <- setdiff(down_genes$gene_name, down_no_atac_genes$gene_name)
tpm_down_prolonged_genes <- avg_and_sd_values[avg_and_sd_values$gene_name %in%
                                              down_atac_genes, ]
avg_tpm_down_prolonged_genes <- tpm_down_prolonged_genes[, grep("_avg",
                                                        colnames(tpm_up_prolonged_genes))]
print("TPM summary of prolonged genes that are downregulated with ATAC peaks in dox")
```

    ## [1] "TPM summary of prolonged genes that are downregulated with ATAC peaks in dox"

``` r
print(summary(avg_tpm_down_prolonged_genes))
```

    ##      0_avg              12_avg              24_avg              48_avg        
    ##  Min.   :  0.0875   Min.   :  0.04292   Min.   :  0.01916   Min.   : 0.01915  
    ##  1st Qu.:  2.1830   1st Qu.:  1.36897   1st Qu.:  1.02910   1st Qu.: 0.80620  
    ##  Median :  5.0506   Median :  2.99619   Median :  2.41355   Median : 1.77973  
    ##  Mean   : 23.7375   Mean   : 14.66228   Mean   : 11.26169   Mean   : 7.68934  
    ##  3rd Qu.: 14.2941   3rd Qu.:  8.09607   3rd Qu.:  7.71749   3rd Qu.: 6.97074  
    ##  Max.   :308.0928   Max.   :188.11414   Max.   :132.48833   Max.   :94.47125  
    ##      96_avg       
    ##  Min.   : 0.0125  
    ##  1st Qu.: 0.7104  
    ##  Median : 1.6535  
    ##  Mean   : 6.7845  
    ##  3rd Qu.: 5.7950  
    ##  Max.   :92.8702

``` r
# plot avg TPM for downregulated-prolonged genes with peaks and no peak
plot_tpm_atac(avg_tpm_down_prolonged_genes, avg_tpm_down_prolonged_genes_nopeak,
              title="Mean TPM ± SEM of Downregulated-Prolonged Genes",
              output_file=file.path(RESULTS_PATH, "bar_tpm_prolonged_down_genes.png"))
```

### Log2FoldChange Plots

``` r
# Generate log fold change plot for upregulated genes
lfc_plot <- plot_lfc_prolonged(
  results, 
  up_atac_genes,
  output_file = file.path(RESULTS_PATH, "lfc_prolonged_up_genes.png")
)

# Generate log fold change plot for downregulated genes
lfc_plot <- plot_lfc_prolonged(
  results, 
  down_no_atac_genes,
  output_file = file.path(RESULTS_PATH, "lfc_prolonged_down_genes.png")
)
```

### TPM Plots

``` r
# Generate TPM plot for upregulated genes
tpm_plot <- plot_tpm_prolonged(
  avg_and_sd_values, 
  up_atac_genes, 
  g2s, 
  filtered_results,
  output_file = file.path(RESULTS_PATH, "tpm_prolonged_up_genes.png")
)

# Generate TPM plot for downregulated genes
tpm_plot <- plot_tpm_prolonged(
  avg_and_sd_values, 
  down_no_atac_genes, 
  g2s, 
  filtered_results,
  output_file = file.path(RESULTS_PATH, "tpm_prolonged_down_genes.png")
)
```

## SECTION 3: IGV

### MAG

Myelin associated glycoprotein
[(MAG)](https://www.uniprot.org/uniprotkb/P20916/entry) is a cell
membrane protein crucial for the formation and maintenance of myelin
sheaths in the peripheral and central nervous system. Additionally, the
absence of this protein leads to diseases such as multiple
sclerosis(MS), suggesting an important role in the integrity and
adequate function of the nervous system\[@Sato1989\].

In our analysis, we found that this gene is down-regulated in humans 96
hours after dox treatment, with no detectable ATAC-seq peak at the final
time point, suggesting reduced chromatin accessibility and gene
activity, hence protein production.This downregulation may compromise
nervous system integrity. However, more research needs to be done in the
effects of dox in MAG production.

![**IGV track with RNA-seq peaks of MAG
gene.**](results/IGV/MAG_gene_track_RNA.png) Peaks represent the
coverage of sequencing reads across this gene. Only peaks from 0 and 96
hours are shown. Note that RNA-seq signal was higher at the 0-hour
timepoint but had diminished by 96 hours, indicating a reduction in gene
expression over time.

<figure>
<img src="results/IGV/MAG_gene_track_ATAC.png"
alt="IGV track with ATAC-seq peaks of MAG gene." />
<figcaption aria-hidden="true"><strong>IGV track with ATAC-seq peaks of
MAG gene.</strong></figcaption>
</figure>

Peaks represent accessible chromatin sites. Samples at 0, 1 and 2-hours
timepoints are shown. Note that most exonic regions of this gene present
accessibility. However, there is one peak in an intronic region. Note
that the epigenetic landscape for this gene does not change overtime.

### AUTS2

The autism susceptibility gene
[(AUTS2)](https://www.uniprot.org/uniprotkb/Q8WXX7/entry) is a component
of a polycomb group multiprotein, polycomb repressive complex I (PRC1)
-like complex\[@Tamburri2020\]. This complexes are required to maintain
the transcriptionally repressive state of many genes, including Hox
genes, which are crucial throughout development\[@Gao2014\]. Mutations
in this gene have been identified in autism patients, and its
suppression in zebrafish embryos was observed to cause microcephaly
\[@Beunders2013\]. This role suggests relevance in gene regulation,
making it an interesting and possibly important gene to further explore
the effects of dox treatment. This gene was upregulated and had
accessible chromatin in our analysis, suggesting an increased activity
after dox treatment, that persisted at least 96 hours.

![**IGV track with RNA-seq peaks of AUTS2
gene.**](results/IGV/AUTS2_gene_track_RNA.png) Peaks represent the
coverage of sequencing reads across this gene. Only peaks from 0 and 96
hours are shown. Note that RNA-seq signal was diminished at 0-hour
timepoint, but higher at the 96-hour timepoint indicating an increased
gene expression over time.

![**IGV track with ATAC-seq peaks of AUTS2
gene.**](results/IGV/AUTS2_gene_track_ATAC.png) Peaks represent
accessible chromatin sites. Samples at 0, 1 and 2-hours timepoints are
shown. Note that AUTS2 gene has multiple accessible points. Note that
the epigenetic landscape for this gene does not change overtime.

### NOD1

The Nucleotide-binding oligomerization domain-containing protein 1
[(NOD1)](https://www.uniprot.org/uniprotkb/Q9Y239/entry) is a
leucine-rich molecule that can regulate both apoptosis and NF-kappaB
activation pathways\[@Inohara1999\]. Given this, it play a crucial role
in innate and adaptive immunity by recognizing Gram-negative bacteria
and viral double stranded RNA. This gene was upregulated and had open
chromatin on our analysis, suggesting an overexpression upon dox
treatment. The upregultaion of this gene might be generating an excess
immune response, perhaps against dox. More research is needed to further
explore this speculation.

![**IGV track with RNA-seq peaks of NOD1
gene.**](results/IGV/NOD1_gene_track_RNA.png) Peaks represent the
coverage of sequencing reads across this gene. Only peaks from 0 and 96
hours are shown. Note that RNA-seq signal was diminished at 0-hour
timepoint, but peaks are higher at the 96-hour timepoint indicating an
increased gene expression over time.

<figure>
<img src="results/IGV/NOD1_gene_track_ATAC.png"
alt="IGV track with ATAC-seq peaks of NOD1 gene." />
<figcaption aria-hidden="true"><strong>IGV track with ATAC-seq peaks of
NOD1 gene.</strong></figcaption>
</figure>

Peaks represent accessible chromatin sites. Samples at 0, 1 and 2-hours
timepoints are shown. Note that multiple exons across the NOD1 gene
present accessible chromatin throughout all timepoints, especially the
most upstream exon. Note that the epigenetic landscape for this gene
does not change overtime.

# References

Beunders, G., Voorhoeve, E., Golzio, C., Pardo, L. M., Rosenfeld, J. A.,
Talkowski, M. E., Simonic, I., Lionel, A. C., Vergult, S., Pyatt, R. E.,
van de Kamp, J., Nieuwint, A., Weiss, M. M., Rizzu, P., Verwer, L. E. N.
I., van Spaendonk, R. M. L., Shen, Y., Wu, B., Yu, T., … Sistermans, E.
A. (2013). Exonic deletions in AUTS2 cause a syndromic form of
intellectual disability and suggest a critical role for the C terminus.
American Journal of Human Genetics, 92(2), 210–220.
<https://doi.org/10.1016/j.ajhg.2012.12.011> Gao, Z., Lee, P., Stafford,
J. M., von Schimmelmann, M., Schaefer, A., & Reinberg, D. (2014). An
AUTS2-Polycomb complex activates gene expression in the CNS. Nature,
516(7531), 349–354. <https://doi.org/10.1038/nature13921> Inohara, N.,
Koseki, T., del Peso, L., Hu, Y., Yee, C., Chen, S., Carrio, R., Merino,
J., Liu, D., Ni, J., & Núñez, G. (1999). Nod1, an Apaf-1-like activator
of caspase-9 and nuclear factor-kappaB. The Journal of Biological
Chemistry, 274(21), 14560–14567.
<https://doi.org/10.1074/jbc.274.21.14560> Mark, K. G., & Rape, M.
(2021). Ubiquitin‐dependent regulation of transcription in development
and disease. EMBO Reports, 22(4), e51078.
<https://doi.org/10.15252/embr.202051078> Monnier, P., Martinet, C.,
Pontis, J., Stancheva, I., Ait-Si-Ali, S., & Dandolo, L. (2013). H19
lncRNA controls gene expression of the Imprinted Gene Network by
recruiting MBD1. Proceedings of the National Academy of Sciences,
110(51), 20693–20698. <https://doi.org/10.1073/pnas.1310201110> Sato,
S., Fujita, N., Kurihara, T., Kuwano, R., Sakimura, K., Takahashi, Y., &
Miyatake, T. (1989). cDNA cloning and amino acid sequence for human
myelin-associated glycoprotein. Biochemical and Biophysical Research
Communications, 163(3), 1473–1480.
<https://doi.org/10.1016/0006-291X(89)91145-5> Tamburri, S., Lavarone,
E., Fernández-Pérez, D., Conway, E., Zanotti, M., Manganaro, D., &
Pasini, D. (2020). Histone H2AK119 Mono-Ubiquitination Is Essential for
Polycomb-Mediated Transcriptional Repression. Molecular Cell, 77(4),
840-856.e5. <https://doi.org/10.1016/j.molcel.2019.11.021>
