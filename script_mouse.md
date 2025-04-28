PowerPuffGenes Report: Mouse stem cells
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
SPECIES          <- "mouse"  # Change to "human" for human data

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
    ## [2] "31787"

``` r
# Filter based on p-value and log2fold changes defined at top
filtered_results <- results %>%
  filter(padj < P_VALUE_FILTER, abs(log2FoldChange) > LOG2_FOLD_FILTER)

print(c("This are the total genes our filtered results with the se pvalue and log2 fold filter",  length(filtered_results$gene_name)))
```

    ## [1] "This are the total genes our filtered results with the se pvalue and log2 fold filter"
    ## [2] "820"

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

    ## [1] "Number of genes with expression sustained until the end: 10"

``` r
print(paste("Number of genes exhibiting transient expression:",
            length(transient_genes$gene_name)))
```

    ## [1] "Number of genes exhibiting transient expression: 464"

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
  show_gene_names = TRUE,
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

    ##                        .
    ## ESC_KO_control_0   82031
    ## ESC_WT_control_0   67410
    ## ESC_KO_control_30  71018
    ## ESC_WT_control_30  68449
    ## ESC_KO_control_60  73631
    ## ESC_WT_control_60  71773
    ## ESC_KO_control_90  94432
    ## ESC_WT_control_90  78247
    ## ESC_KO_control_120 95582
    ## ESC_WT_control_120 79053
    ## ESC_KO_control_150 85303
    ## ESC_WT_control_150 77512

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
    ## [2] "53747"

Dox samples

``` r
# common peaks in dox time points (!not time 0)
dox_samples <- names(my_consensus_peaks)[!grepl("_0$", names(my_consensus_peaks))]
dox_peaks <- my_consensus_peaks[dox_samples]
dox_common_peaks <- suppressWarnings(find_common_peaks(dox_peaks))

print(c("This is how many peaks are common in dox:", length(dox_common_peaks)))
```

    ## [1] "This is how many peaks are common in dox:"
    ## [2] "38880"

Overlap of dox and non-dox (to get unique to each condition)

``` r
# Now overlap between dox and non-dox common peaks
dox_compare_list <- list(non_dox = non_dox_common_peaks, dox = dox_common_peaks)
dox_non_dox_ov <- suppressWarnings(find_common_peaks(dox_compare_list))
print(c("This is how many peaks are common in both non- and dox treatments", length(dox_non_dox_ov)))
```

    ## [1] "This is how many peaks are common in both non- and dox treatments"
    ## [2] "37348"

``` r
# extracting peaks unique to each condition (dox non-dox)
# Peaks unique to non_dox
unique_to_non_dox <- suppressWarnings(find_my_peaks(dox_non_dox_ov, non_dox_common_peaks))
print(c("This is how many peaks are unique to non-dox condition", length(unique_to_non_dox)))
```

    ## [1] "This is how many peaks are unique to non-dox condition"
    ## [2] "16399"

``` r
# Peaks unique to dox
unique_to_dox <- suppressWarnings(find_my_peaks(dox_non_dox_ov, dox_common_peaks))

print(c("This is how many peaks are unique to dox condition", length(unique_to_dox)))
```

    ## [1] "This is how many peaks are unique to dox condition"
    ## [2] "1774"

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
    ## [2] "18753"

``` r
gr_list_gene_promoter_dox_ov <-  list(gene_promoters = gene_promoters, dox_peaks = unique_to_dox)
unique_to_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_dox_ov))

print(c("This is how many dox unique peaks overlapped gene promoters", length(unique_to_dox_gene_promoter_ov)))
```

    ## [1] "This is how many dox unique peaks overlapped gene promoters"
    ## [2] "251"

``` r
# peaks common to non-dox condition overlapped with gene promoters
gr_list_gene_promoter_dox_ov <- list(gene_promoters = gene_promoters, dox_peaks = non_dox_common_peaks)
non_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_dox_ov))

print(c("This is how many non_dox common peaks overlapped gene promoters", length(non_dox_gene_promoter_ov)))
```

    ## [1] "This is how many non_dox common peaks overlapped gene promoters"
    ## [2] "20412"

``` r
gene_promoter_lostindox <- setdiff(non_dox_gene_promoter_ov$gene_name,
                                   dox_gene_promoter_ov$gene_name)
print(c("This is how many dox lost peaks overlapped gene promoters", length(gene_promoter_lostindox)))
```

    ## [1] "This is how many dox lost peaks overlapped gene promoters"
    ## [2] "1961"

### Matching ATAC peaks to transient and prolonged genes

``` r
# filter common dox peaks genes for transient genes
gene_promoters_trans <- intersect(dox_gene_promoter_ov$gene_name, transient_genes$gene_name)
print(c("this is how many transient genes overlapped with common dox gene promoters", length(gene_promoters_trans)))
```

    ## [1] "this is how many transient genes overlapped with common dox gene promoters"
    ## [2] "153"

``` r
# get tpm for transient genes with dox gene promoter peaks
tpm_trans_genes <- avg_and_sd_values[avg_and_sd_values$gene_name %in% 
                                    gene_promoters_trans, ]
avg_tpm_trans_genes <- tpm_trans_genes[, grep("_avg", colnames(tpm_trans_genes))]
print("TPM summary of transient genes overlapped with common dox gene promoters")
```

    ## [1] "TPM summary of transient genes overlapped with common dox gene promoters"

``` r
print(summary(avg_tpm_trans_genes))
```

    ##      0_avg               12_avg             24_avg             48_avg        
    ##  Min.   :   0.0184   Min.   :   0.000   Min.   :   0.074   Min.   :   0.082  
    ##  1st Qu.:   0.7844   1st Qu.:   0.687   1st Qu.:   0.523   1st Qu.:   0.424  
    ##  Median :   1.4463   Median :   1.475   Median :   1.391   Median :   1.163  
    ##  Mean   :  20.9641   Mean   :  43.839   Mean   :  42.185   Mean   :  45.382  
    ##  3rd Qu.:   3.1572   3rd Qu.:   3.649   3rd Qu.:   3.538   3rd Qu.:   2.849  
    ##  Max.   :2566.2581   Max.   :6088.390   Max.   :5928.518   Max.   :6447.627  
    ##      96_avg        
    ##  Min.   :   0.153  
    ##  1st Qu.:   0.591  
    ##  Median :   1.237  
    ##  Mean   :  30.368  
    ##  3rd Qu.:   3.544  
    ##  Max.   :4157.756

``` r
# filter unique dox peaks genes for transient genes
unique_gene_promoters_trans <- intersect(unique_to_dox_gene_promoter_ov$gene_name, transient_genes$gene_name)
print(c("this is how many transient genes overlapped with unique dox gene promoters", length(unique_gene_promoters_trans)))
```

    ## [1] "this is how many transient genes overlapped with unique dox gene promoters"
    ## [2] "3"

``` r
if (length(unique_gene_promoters_trans) > 0) {
  print(unique_gene_promoters_trans)
}
```

    ## [1] "Neurod1"  "Prcd"     "Rps2-ps8"

``` r
# get tpm for transient genes with unique dox gene promoter peaks
unique_tpm_trans_genes <- avg_and_sd_values[avg_and_sd_values$gene_name %in% 
                                    unique_gene_promoters_trans, ]
unique_avg_tpm_trans_genes <- unique_tpm_trans_genes[, grep("_avg", colnames(unique_tpm_trans_genes))]
print("TPM summary of transient genes overlapped with unique dox gene promoters")
```

    ## [1] "TPM summary of transient genes overlapped with unique dox gene promoters"

``` r
print(summary(unique_avg_tpm_trans_genes))
```

    ##      0_avg           12_avg          24_avg           48_avg      
    ##  Min.   :1.828   Min.   :1.627   Min.   :0.7289   Min.   :0.5369  
    ##  1st Qu.:1.885   1st Qu.:1.889   1st Qu.:0.9000   1st Qu.:0.5938  
    ##  Median :1.942   Median :2.151   Median :1.0711   Median :0.6506  
    ##  Mean   :1.947   Mean   :2.045   Mean   :1.6185   Mean   :1.3397  
    ##  3rd Qu.:2.007   3rd Qu.:2.254   3rd Qu.:2.0633   3rd Qu.:1.7412  
    ##  Max.   :2.073   Max.   :2.358   Max.   :3.0556   Max.   :2.8317  
    ##      96_avg     
    ##  Min.   :1.031  
    ##  1st Qu.:1.059  
    ##  Median :1.087  
    ##  Mean   :1.158  
    ##  3rd Qu.:1.222  
    ##  Max.   :1.357

``` r
# filter upregulated-prolonged genes
up_genes <- genes_of_interest %>%
  filter(tpm_log2fc > 0)
print(c("this is how many prolonged genes that are upregulated with dox", length(up_genes$gene_name)))
```

    ## [1] "this is how many prolonged genes that are upregulated with dox"
    ## [2] "3"

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

    ##      0_avg            12_avg           24_avg           48_avg     
    ##  Min.   : 3.921   Min.   : 3.704   Min.   : 7.351   Min.   :15.80  
    ##  1st Qu.: 5.916   1st Qu.: 5.993   1st Qu.:10.061   1st Qu.:18.93  
    ##  Median : 7.911   Median : 8.282   Median :12.772   Median :22.06  
    ##  Mean   :10.543   Mean   :10.267   Mean   :17.730   Mean   :25.34  
    ##  3rd Qu.:13.853   3rd Qu.:13.548   3rd Qu.:22.919   3rd Qu.:30.11  
    ##  Max.   :19.795   Max.   :18.814   Max.   :33.066   Max.   :38.16  
    ##      96_avg     
    ##  Min.   :13.94  
    ##  1st Qu.:15.49  
    ##  Median :17.05  
    ##  Mean   :23.35  
    ##  3rd Qu.:28.06  
    ##  Max.   :39.07

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
    ## [2] "2"

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

    ##      0_avg            12_avg           24_avg          48_avg     
    ##  Min.   : 7.911   Min.   : 8.282   Min.   :12.77   Min.   :22.06  
    ##  1st Qu.:10.882   1st Qu.:10.915   1st Qu.:17.85   1st Qu.:26.08  
    ##  Median :13.853   Median :13.548   Median :22.92   Median :30.11  
    ##  Mean   :13.853   Mean   :13.548   Mean   :22.92   Mean   :30.11  
    ##  3rd Qu.:16.824   3rd Qu.:16.181   3rd Qu.:27.99   3rd Qu.:34.13  
    ##  Max.   :19.795   Max.   :18.814   Max.   :33.07   Max.   :38.16  
    ##      96_avg     
    ##  Min.   :17.05  
    ##  1st Qu.:22.55  
    ##  Median :28.06  
    ##  Mean   :28.06  
    ##  3rd Qu.:33.57  
    ##  Max.   :39.07

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

    ##      0_avg           12_avg          24_avg          48_avg         96_avg     
    ##  Min.   :3.921   Min.   :3.704   Min.   :7.351   Min.   :15.8   Min.   :13.94  
    ##  1st Qu.:3.921   1st Qu.:3.704   1st Qu.:7.351   1st Qu.:15.8   1st Qu.:13.94  
    ##  Median :3.921   Median :3.704   Median :7.351   Median :15.8   Median :13.94  
    ##  Mean   :3.921   Mean   :3.704   Mean   :7.351   Mean   :15.8   Mean   :13.94  
    ##  3rd Qu.:3.921   3rd Qu.:3.704   3rd Qu.:7.351   3rd Qu.:15.8   3rd Qu.:13.94  
    ##  Max.   :3.921   Max.   :3.704   Max.   :7.351   Max.   :15.8   Max.   :13.94

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
    ## [2] "7"

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

    ##      0_avg            12_avg           24_avg           48_avg      
    ##  Min.   :0.3459   Min.   :0.2742   Min.   :0.1946   Min.   :0.1628  
    ##  1st Qu.:1.0524   1st Qu.:0.9532   1st Qu.:0.4359   1st Qu.:0.4233  
    ##  Median :1.4309   Median :1.2785   Median :0.6628   Median :0.5662  
    ##  Mean   :2.9827   Mean   :2.4540   Mean   :1.0578   Mean   :1.0155  
    ##  3rd Qu.:3.8831   3rd Qu.:3.9738   3rd Qu.:1.2912   3rd Qu.:1.1075  
    ##  Max.   :9.2310   Max.   :5.7713   Max.   :3.0932   Max.   :3.3178  
    ##      96_avg      
    ##  Min.   :0.1434  
    ##  1st Qu.:0.3983  
    ##  Median :0.5995  
    ##  Mean   :1.0293  
    ##  3rd Qu.:1.1457  
    ##  Max.   :3.3746

``` r
lost_gene_promoters <- intersect(gene_promoter_lostindox, down_genes$gene_name)
print(c("this is how many downregulated-prolonged genes overlapped with the lost gene promoters in dox", length(lost_gene_promoters)))
```

    ## [1] "this is how many downregulated-prolonged genes overlapped with the lost gene promoters in dox"
    ## [2] "0"

``` r
if (length(lost_gene_promoters) > 0) {
  print(lost_gene_promoters)
}

# filter downregulated-prolonged genes to non dox promoter overlaps
down_prolonged_rnaseq_atac_dox <- dox_gene_promoter_ov[dox_gene_promoter_ov$gene_name %in% down_genes$gene_name]
down_no_atac_genes <- down_genes %>%
  filter(!(gene_name %in% down_prolonged_rnaseq_atac_dox$gene_name))

print(c("this is how many downregulated genes that couldn't find ATACseq peaks in dox", length(down_prolonged_rnaseq_atac_dox)))
```

    ## [1] "this is how many downregulated genes that couldn't find ATACseq peaks in dox"
    ## [2] "3"

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

    ##      0_avg            12_avg          24_avg           48_avg      
    ##  Min.   :0.7429   Min.   :0.672   Min.   :0.3274   Min.   :0.2965  
    ##  1st Qu.:1.2589   1st Qu.:1.127   1st Qu.:0.4901   1st Qu.:0.4988  
    ##  Median :2.1236   Median :2.163   Median :0.7226   Median :0.6900  
    ##  Mean   :2.4850   Mean   :2.475   Mean   :0.8636   Mean   :0.7694  
    ##  3rd Qu.:3.3497   3rd Qu.:3.511   3rd Qu.:1.0960   3rd Qu.:0.9606  
    ##  Max.   :4.9501   Max.   :4.900   Max.   :1.6816   Max.   :1.4011  
    ##      96_avg      
    ##  Min.   :0.2245  
    ##  1st Qu.:0.5058  
    ##  Median :0.8154  
    ##  Mean   :0.7788  
    ##  3rd Qu.:1.0885  
    ##  Max.   :1.2600

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

    ##      0_avg            12_avg           24_avg           48_avg      
    ##  Min.   :0.3459   Min.   :0.2742   Min.   :0.1946   Min.   :0.1628  
    ##  1st Qu.:0.8539   1st Qu.:0.7543   1st Qu.:0.4287   1st Qu.:0.3565  
    ##  Median :1.3619   Median :1.2343   Median :0.6628   Median :0.5502  
    ##  Mean   :3.6463   Mean   :2.4266   Mean   :1.3169   Mean   :1.3436  
    ##  3rd Qu.:5.2965   3rd Qu.:3.5028   3rd Qu.:1.8780   3rd Qu.:1.9340  
    ##  Max.   :9.2310   Max.   :5.7713   Max.   :3.0932   Max.   :3.3178  
    ##      96_avg      
    ##  Min.   :0.1434  
    ##  1st Qu.:0.3577  
    ##  Median :0.5720  
    ##  Mean   :1.3633  
    ##  3rd Qu.:1.9733  
    ##  Max.   :3.3746

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

## SECTION 3: IGV of interesting genes

### Usp26

The [Usp26](https://www.uniprot.org/uniprotkb/Q99MX1/entry) gene has
been involved in deubiquitination pathways, making it potentially
relevant for gene regulation. There are many transcription factors in
the cell, such as polycomb repressive complex 1 (PRC1), that are in
charge of maintaining and regulating epigenetic marks in the genome
(Tamburri et al., 2020). Many of these transcription factors, need
ubiquitination for their function, making a deubiquitination crucial for
this process as well (Mark & Rape, 2021). For example, this gene has
been found to be involved in somatic cell reprogramming through the K48
deubiquitination of two protein components of PRC1.

<figure>
<img src="results/IGV/Usp26_IGV_track2.png"
alt="IGV track with RNA-seq peaks of Usp26 gene." />
<figcaption aria-hidden="true"><strong>IGV track with RNA-seq peaks of
Usp26 gene.</strong></figcaption>
</figure>

Peaks represent the coverage of sequencing reads across this gene. Note
that peaks are higher at 0 and 12 hours, but are lower when 24 hours are
reached. They keep beung low until the last time point of 96 hours. Note
there are higher expression signals in the nost downstream exon of the
Usp26 gene.

<figure>
<img src="results/IGV/Usp26_gene_track_ATAC.png"
alt="IGV track with ATAC-seq peaks of Usp26 gene." />
<figcaption aria-hidden="true"><strong>IGV track with ATAC-seq peaks of
Usp26 gene.</strong></figcaption>
</figure>

Peaks represent accessible chromatin sites. Samples at 0, 30, 60, 90,
120 and 150 minutes after dox treatment are shown. Note that the
epigenetic landscape for this gene seems to be changing overtime, with
its higher peak presented at 90 minute after treatment timepoint. This
IGV track is showin the accessibility of Usp26 gene in its most upstream
exon. Note that this exon is not the same as the one in the RNA seq
track.

# References

Mark, K. G., & Rape, M. (2021). Ubiquitin‐dependent regulation of
transcription in development and disease. EMBO Reports, 22(4), e51078.
<https://doi.org/10.15252/embr.202051078>

Tamburri, S., Lavarone, E., Fernández-Pérez, D., Conway, E., Zanotti,
M., Manganaro, D., & Pasini, D. (2020). Histone H2AK119
Mono-Ubiquitination Is Essential for Polycomb-Mediated Transcriptional
Repression. Molecular Cell, 77(4), 840-856.e5.
<https://doi.org/10.1016/j.molcel.2019.11.021>
