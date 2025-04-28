PowerPuffGenes Report
================
2025-04-01

- [Introduction](#introduction)
- [Goals](#goals)
- [Methods](#methods)
  - [Nextflow output](#nextflow-output)
- [Data preprocessing for RNA-seq analysis with
  DESeq2](#data-preprocessing-for-rna-seq-analysis-with-deseq2)
  - [Parameters and paths definition](#parameters-and-paths-definition)
  - [Counts loading and preprocesing](#counts-loading-and-preprocesing)
  - [Metadata preparation](#metadata-preparation)
  - [DESeq execution](#deseq-execution)
  - [DESeq time-point arrangement](#deseq-time-point-arrangement)
  - [Filtering and visualization](#filtering-and-visualization)
  - [TPM analysis](#tpm-analysis)

# Introduction

- Doxycycline (or Dox) is an antibiotic drug used in biological research
  for inducing gene expression. intro in dox, experiment For this
  analysis, we used wild type mouse (Mus musculus) stem cells at
  different time points after Dox treatment was applied. We used three
  replicates for each time point: 0, 12, 24, 48 and 96 hours after
  treatment.

# Goals

- Identify specific genes that are significantly up or downregulated
  after Dox treatment at various time points.
- General analysis of gene expression. Consider repeating the analysis
  from class but normalizing based on time 0.

# Methods

\##Nextflow RNA Sequencing Pipeline - For our analysis, we used
Nextflow’s RNAseq pipeline versions
\[3.14.0.\]\[<https://nf-co.re/rnaseq/3.14.0>\] and
\[3.18.0.\]\[<https://nf-co.re/rnaseq/3.18.0>\]

- The inputs we used for running this pipeline:

- Raw fasta sample sequence files

- Mouse reference genome assembled sequence (GRCm38/mm10)

- Genome annotation

- This pipeline consists of a series of steps that process our samples
  to give an output raw read counts, and transcripts per million(tpm).
  Counts and tpm values were used to

## Nextflow output

- May describe the results from the nextflow pipeline, some QC plots and
  the type of files produced.

\##DESeq2 Analysis - Explain in general the overall idea with RNA-seq
analysis

# Data preprocessing for RNA-seq analysis with DESeq2

### Parameters and paths definition

``` r
COUNTS_PATH      <- "/scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_rnaseq/salmon.merged.gene_counts.tsv"
TPM_PATH         <- "/scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_rnaseq/salmon.merged.gene_tpm.tsv"
RESULT_PATH      <- "../result_human"
FILTER_TRESHOLD  <- 0
P_VALUE_FILTER   <- 0.01
LOG2_FOLD_FILTER <- 1
COLUMN_PREFIX    <- "gfp"
```

### Counts loading and preprocesing

``` r
# Read raw count data and make dictionary mapping gene id with gene name
counts_orig  <- read.table(COUNTS_PATH, header=TRUE, row.names=1)
g2s          <- data.frame(gene_id = rownames(counts_orig), 
                           gene_name = counts_orig[ , 1]) 
write.csv(g2s, file.path(RESULT_PATH, "g2s.csv"))

# Remove gene name column from counts and convert counts to numerical rounded matrix
counts           <- counts_orig %>% select(-gene_name)
counts_matrix    <- as.matrix(counts) 
counts_matrix    <- round(counts_matrix)

# Filter genes with no counts for any of the samples
counts_filtered  <- counts_matrix[rowSums(counts_matrix) > FILTER_TRESHOLD, ]
```

### Metadata preparation

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
```

### DESeq execution

``` r
# Prepare DESeq dataset and run DESea2
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = deseq_samples,
                              design = ~ time_point) 
dds <- DESeq(dds)
```

### DESeq time-point arrangement

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

# Create empty df to store results values (schema looks good)
results <- data.frame("gene_id" = character(),
                      "baseMean" = numeric(),
                      "log2FoldChange" = numeric(),
                      "lfcSE" = numeric(),
                      "stat" = numeric(),
                      "pvalue" = numeric(),
                      "padj" = numeric(),
                      "gene_name" = character(),
                      "comparison" = character())

# Loop through the defined contrasts
for(i in names(comparison_list)) {
  current_comparison <- comparison_list[[i]] # get comparison i
  res <- results(dds, contrast = current_comparison) # get results for time comparison i using 'contrast'
  if (is.null(res)) {
      warning(paste("Could not generate results for contrast:", contrast_name))
      next}
  # Temporary df to store the results for comparison i
  tmp_res_df <- res %>% 
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    merge(g2s, by = "gene_id", all.x = TRUE) %>% 
    mutate(comparison = i) 
  # Add the temporary df to the main results df
  results <- bind_rows(results, tmp_res_df)
}
```

### Filtering and visualization

``` r
# Filter based on p-value < 0.05 and log2 fold change > 1
filtered_results <- results %>%
  filter(padj < P_VALUE_FILTER, abs(log2FoldChange) > LOG2_FOLD_FILTER)

# Get all gene names that are significant (drop gene name repetitions)
filtered_genes <- as.data.frame(filtered_results$gene_name, collapse = "\n")
filtered_genes <- unique(filtered_genes)
colnames(filtered_genes)[1] <- "gene_name"

# Now let's write this out and do gene enrichment analysis
write.table(filtered_genes["gene_name"], row.names = FALSE, col.names = FALSE, file.path(RESULT_PATH, "filtered_genes.csv"))

save(results, filtered_results, filtered_genes, deseq_samples, g2s, file = file.path(RESULT_PATH, "DESeq2_results.RData"))
```

``` r
# Distribution of baseMean, lfcSE, p-values and fold change
hist(filtered_results$padj, xlim = c(0,P_VALUE_FILTER ), breaks = 50)
```

![](01_preprocessing_files/figure-gfm/filtered%20results%20distributions-1.png)<!-- -->

``` r
hist(filtered_results$log2FoldChange, xlim = c(-5, 5), breaks = 1000)
```

![](01_preprocessing_files/figure-gfm/filtered%20results%20distributions-2.png)<!-- -->

``` r
# T-test comparing baseMean expression of all genes vs filtered genes
basemean_all_genes <- median(results$baseMean)
basemean_fil_genes <- median(filtered_results$baseMean)
t.test(results$baseMean, filtered_results$baseMean)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  results$baseMean and filtered_results$baseMean
    ## t = 11.661, df = 3875.5, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  229.6872 322.5305
    ## sample estimates:
    ## mean of x mean of y 
    ##  610.6926  334.5837

## TPM analysis

``` r
# Read raw tpm data and delete gene name column
tpm_orig  <- read.table(TPM_PATH, header=TRUE, row.names=1)
tpm       <- tpm_orig %>% select(-gene_name)

# Filter genes with no counts for any of the samples
tpm_filtered <- tpm[rowSums(tpm) > FILTER_TRESHOLD, ]

# Remove extra timepoints and replicates (tpm)
tpm_filtered <- tpm_filtered[, colnames(tpm_filtered) %in% keep_cols, drop = FALSE]
```

``` r
# Loop to calculate avg and sd for replicates at a given time point
time_points <- c("0", "12", "24", "48", "96")
avg_and_sd_values <- list()
for (tp in time_points) {
  cols     <- grep(paste0(COLUMN_PREFIX, "_", tp, "_"), colnames(tpm_filtered))
  avg      <- rowMeans(tpm_filtered[, cols])
  sd       <- apply(tpm_filtered[, cols], 1, sd)
  sd       <- data.frame(sd)
  combined <- cbind(avg, sd)
  avg_and_sd_values <- c(avg_and_sd_values, list(combined))
}

# Convert the list to a data frame and add column names for the respective timepoint
avg_and_sd_values           <- do.call(cbind,  avg_and_sd_values)
colnames(avg_and_sd_values) <- paste0(rep(time_points, each = 2), c("_avg", "_sd"))
avg_and_sd_values           <- as.data.frame(avg_and_sd_values) %>% rownames_to_column("gene_id") %>% merge(g2s)

save(tpm_filtered, avg_and_sd_values, file = file.path(RESULT_PATH, "TPM_results.RData"))
```

``` r
# Make a distribution of TPM values in each time point - using facet-wrap
time_values <- names(avg_and_sd_values)[grep("avg", names(avg_and_sd_values))]
melt_tpm_df <- melt(avg_and_sd_values, measure.vars = time_values, value_name = "gene_id")
ggplot(melt_tpm_df, aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ variable, scales = "free_x")+
  xlim(0,500) +
  ylim(0,500) +
  theme("paperwhite")
```

![](01_preprocessing_files/figure-gfm/facet%20all%20timepoint%20TPM%20distributions-1.png)<!-- -->
