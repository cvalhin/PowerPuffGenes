##################################################
###### Process single result comparison ##########
##################################################
process_comparison <- function(dds_obj, comparison_params, comparison_name, gene_map) {
  res <- results(dds_obj, contrast = comparison_params)
  
  if (is.null(res)) {
    warning(paste("Could not generate results for contrast:", comparison_name))
    return(NULL)
  }
  
  res %>% 
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    merge(gene_map, by = "gene_id", all.x = TRUE) %>% 
    mutate(comparison = comparison_name)
}

##################################################
### Categorize transient vs prolongued genes #####
##################################################
categorize_gene_expression <- function(filtered_results) {
  # Get baseline comparisons (vs_0)
  filtered_vs_zero <- filtered_results[grep("vs_0", filtered_results$comparison), ]
  
  # Identify prolonged genes (significant at multiple time points including time 96)
  prolonged_res <- filtered_vs_zero %>% 
    group_by(gene_id) %>% 
    filter(n() > 1)
  
  prolonged_genes <- prolonged_res[grep("96_vs", prolonged_res$comparison),
                                   which(names(prolonged_res) == "gene_name")]
  
  prolonged_res_results <- merge(filtered_vs_zero, prolonged_genes)
  
  # Group prolonged genes by their triggering time points
  change_at_12 <- prolonged_res_results[grep("12_vs", prolonged_res_results$comparison),
                                        which(names(prolonged_res_results) == "gene_name")]
  change_at_24 <- prolonged_res_results[grep("24_vs", prolonged_res_results$comparison),
                                        which(names(prolonged_res_results) == "gene_name")]
  change_at_48 <- prolonged_res_results[grep("48_vs", prolonged_res_results$comparison),
                                        which(names(prolonged_res_results) == "gene_name")]
  
  # Find genes changed from specific time points to end
  change_at_12_to_end <- Reduce(intersect, list(change_at_12, change_at_24, change_at_48))
  change_at_24_to_end <- setdiff(intersect(change_at_24, change_at_48), change_at_12_to_end)
  change_at_48_to_end <- setdiff(change_at_48, 
                                 union(change_at_12_to_end, change_at_24_to_end))
  
  # Identify transient genes (only significant at one time point)
  transient_genes <- filtered_vs_zero[!(filtered_vs_zero$gene_name %in% prolonged_genes$gene_name),
                                      which(names(filtered_vs_zero) == "gene_name")]
  
  # Combine results into a list
  return(list(
    change_at_12_to_end = change_at_12_to_end,
    change_at_24_to_end = change_at_24_to_end,
    change_at_48_to_end = change_at_48_to_end,
    transient_genes = transient_genes
  ))
}

##################################################
############# Prepare Heatmap Data ###############
##################################################
prepare_heatmap_data <- function(tpm_matrix, 
                                 gene_symbol_map, 
                                 gene_filter = NULL,
                                 log_transform = TRUE,
                                 pseudo_count = 1,
                                 handle_duplicates = "keep_first") {
  
  # Apply log transformation if requested
  if (log_transform) {
    transformed_matrix <- log2(tpm_matrix + pseudo_count)
  } else {
    transformed_matrix <- tpm_matrix
  }
  
  # Add gene IDs and names
  matrix_with_genes <- transformed_matrix %>% 
    rownames_to_column("gene_id") %>% 
    merge(gene_symbol_map)
  
  # Filter genes if a filter is provided
  if (!is.null(gene_filter)) {
    if (is.data.frame(gene_filter)) {
      # If gene_filter is a data frame, extract gene names
      gene_filter <- gene_filter[[1]]
    }
    
    # Filter to genes of interest
    filtered_matrix <- matrix_with_genes[matrix_with_genes$gene_name %in% gene_filter, ]
  } else {
    filtered_matrix <- matrix_with_genes
  }
  
  # Handle duplicate gene names
  if (handle_duplicates == "keep_first") {
    # Keep only the first occurrence of each gene name
    filtered_matrix <- filtered_matrix %>%
      group_by(gene_name) %>%
      slice(1) %>%
      ungroup()
  } else if (handle_duplicates == "make_unique") {
    # Make gene names unique by appending a suffix
    filtered_matrix <- filtered_matrix %>%
      group_by(gene_name) %>%
      mutate(gene_display = if(n() > 1) paste0(gene_name, "_", row_number()) else gene_name) %>%
      ungroup()
  }
  
  # Remove rows with NA values
  filtered_matrix <- na.omit(filtered_matrix)
  
  # Create numeric matrix for heatmap
  numeric_matrix <- filtered_matrix %>% select(-gene_id, -gene_name)
  
  # If we made unique names, use those
  if (handle_duplicates == "make_unique" && "gene_display" %in% names(filtered_matrix)) {
    rownames(numeric_matrix) <- filtered_matrix$gene_display
  } else {
    # Store gene names as row names
    rownames(numeric_matrix) <- filtered_matrix$gene_name
  }
  
  return(list(
    full_data = filtered_matrix,
    numeric_matrix = numeric_matrix
  ))
}

##################################################
########### Prepare Gene Annotations #############
##################################################
prepare_gene_annotations <- function(gtf_file, 
                                     valid_chromosomes = paste0("chr", c(1:19, "X", "Y"))) {
  # Get gene annotations from gtf file
  gtf_data <- rtracklayer::import(gtf_file)
  
  # Extract gene locations
  gene_locations <- as.data.frame(gtf_data) %>%
    filter(type == "gene") %>%
    dplyr::select(chr = seqnames, start, end, gene_id, gene_name) %>%
    mutate(chr = as.character(chr))
  
  # Add 'chr' prefix if missing
  if (!all(str_detect(gene_locations$chr, "^chr"))) {
    gene_locations$chr <- paste0("chr", gene_locations$chr)
  }
  
  # Filter chromosomes
  gene_locations <- filter(gene_locations, chr %in% valid_chromosomes)
  
  # Remove versions of gene annotations
  gene_locations <- gene_locations %>%
    mutate(gene_id_base = str_remove(gene_id, "\\.\\d+$"))
  
  return(gene_locations)
}


##################################################
######## Merge results with annotations ##########
##################################################
merge_results_with_annotations <- function(de_results, gene_locations) {
  # Remove versions from results
  de_results_prepared <- de_results %>%
    mutate(gene_id_base = str_remove(gene_id, "\\.\\d+$"))
  
  # Merge DE results with gene locations
  annotated_results <- inner_join(de_results_prepared, gene_locations, by = "gene_id_base")
  
  # Select and rename relevant columns
  annotated_results <- annotated_results %>%
    dplyr::select(
      gene_id = gene_id.x, 
      log2FoldChange, 
      padj, 
      gene_name = gene_name.x,
      comparison, 
      chr, 
      start, 
      end
    )
  
  return(annotated_results)
}

##################################################
############## t-test for one gene ###############
##################################################
calculate_stats <- function(gene_row) {
  tpm_0 <- as.numeric(gene_row[t0_cols])
  mean_0 = mean(tpm_0)
  tpm_96 <- as.numeric(gene_row[t96_cols])
  mean_96 = mean(tpm_96)
  log2FC = log2(mean(tpm_96) + 1e-3) - log2(mean(tpm_0) + 1e-3)
  t_result <- t.test(tpm_0, tpm_96, paired = FALSE)
  return(data.frame( pvalue = t_result$p.value,
                     tpm_log2fc = log2FC)
  )
}