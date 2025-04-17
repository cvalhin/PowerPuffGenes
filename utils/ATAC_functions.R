###################################
# ATACSEQ IMPORT PEAKS
###################################

#' FUNCTION to import peak files from ATACseq into List of GRanges objects
#'
#' @description 
#' Any peak file can be loaded in and made into a GRanges from one file 
#' to thouusands of peak files. This function will load and get into GRange format.
#' 
#' @param peak_file, file path to the peak files - generic name "consensus_file_path".  
#' @param  regex, # Warning custom to application # there is a regex to get "sample name" for input into function 
#' 
#' 
import_peaks <- function(consensus_file_path) {
  
  # List all MACS2 broadPeak files in the directory
  peak_files <- list.files(consensus_file_path, full.names = TRUE, pattern = ".broadPeak")
  
  # Extract sample names directly from file names
  names <- str_extract(basename(peak_files), ".*(?=.mLb.clN_peaks.broadPeak)")

  # Apply transformation
  sample_names <- sapply(names, function(x) {
    if (grepl("[0-9.]+h_", x)) {
      # Extract number before 'h_'
      hour_str <- str_extract(x, "[0-9.]+(?=h_)")
      mins <- as.integer(as.numeric(hour_str) * 60)
      # Replace "<number>h" with "<minutes>"
      x <- sub("([0-9.]+)h", mins, x)
    }
    return(x)
  })
  
  # Initialize an empty list to store filtered GRanges objects
  peak_list <- list()
  
  # Import, filter, and store GRanges objects
  for (i in seq_along(peak_files)) {
    # Import peaks as GRanges
    peaks <- rtracklayer::import(peak_files[i])
    
    # Filter to keep only seqnames starting with "chr"
    filtered_peaks <- peaks[grepl("^chr", as.character(seqnames(peaks)))]
    
    # Add filtered peaks to the list with a sample-specific name
    peak_list[[sample_names[i]]] <- filtered_peaks
  }
  
  return(peak_list)
}

###################################
# FIND CONSENSUS PEAKS
###################################

#' Find Consensus Peaks Across Samples
#'
#' @description
#' This function identifies consensus peaks shared across multiple samples at each time point.
#' The consensus peaks are defined as those that are common across samples from the same time 
#' point. If a sample group has only one replicate, it returns the peaks directly for that group.
#' 
#' @param peak_list A list of GRanges objects, each containing peaks for a different sample. The list
#'   is named by sample names, which should contain time points in the format "_<time>_".
#' 
#' @return A list of consensus peaks for each unique sample group across time points. 
#'   Each entry corresponds to a sample name, and if multiple replicates exist, consensus peaks
#'   are calculated using the `find_common_peaks` function.
#' 
#' @examples
#' consensus_peaks <- find_consensus_peaks(peak_list)
#'
find_consensus_peaks <- function(peak_list) {
  
  # Extract unique time points based on sample names in the peak_list
  time_points <- sort(unique(as.numeric(str_extract(names(peak_list), "(?<=_)[0-9]+(?=_REP)"))))
  
  # Initialize an empty list to store consensus peaks
  consensus_peaks <- list()
  
  # Loop over each time point to find consensus peaks
  for (i in seq_along(time_points)) {
    t <- time_points[i]
    
    # Find sample names corresponding to the current time point
    cols <- names(peak_list)[grep(paste0("_", t ,"_"), names(peak_list))]
    col_names <- unique(str_remove(cols, "_REP[0-9]+$"))  # Remove replicate numbers
    
    # Loop over each sample group and compute consensus peaks
    for (j in seq_along(col_names)) {
      sample_group <- col_names[j]
      
      # If there are multiple replicates for the sample group, find common peaks
      if (length(grep(sample_group, names(peak_list))) > 1) {
        consensus_peaks[[sample_group]] <- suppressWarnings(
          find_common_peaks(peak_list[c(grep(sample_group, names(peak_list)))]))
      } else {
        # If there's only one replicate, keep the peaks for that sample group
        consensus_peaks[sample_group] <- peak_list[c(grep(sample_group, names(peak_list)))]
      }
    }
  }
  
  return(consensus_peaks)
}


###################################
# FIND COMMON PEAKS
###################################

#' FUNCTION to identify common peaks across multiple ATAC-seq samples
#'
#' @description 
#' This function processes a list of GRanges objects, where each GRanges represents 
#' ATAC-seq peaks from a single sample. It identifies genomic regions (peaks) 
#' that are shared across all samples by iteratively calculating overlaps. 
#' The resulting common peaks are labeled with unique identifiers for downstream 
#' analysis or visualization in tools like IGV.
#'
#' @param gr_list A list of GRanges objects. Each GRanges represents the peaks from 
#'                a single ATAC-seq sample. The list should have names corresponding 
#'                to the sample identifiers (e.g., "sample1", "sample2", etc.).
#' 
#' @return A GRanges object containing the genomic intervals (peaks) that are common 
#'         across all input GRanges objects. Each interval is labeled with a unique 
#'         identifier in the metadata column `name`, using the format "common_peak_<number>".
#'
#' @details 
#' 1. The function validates that the input is a list of GRanges objects.
#' 2. It iteratively identifies overlaps across all GRanges objects in the list.
#' 3. Only intervals shared across all samples are retained at each step.
#' 4. The resulting peaks are assigned unique names for easy tracking.
#'
#' @examples
#' # Example list of GRanges objects
#' gr_list <- list(
#'   sample1 = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 50), end = c(20, 70))),
#'   sample2 = GRanges(seqnames = "chr1", ranges = IRanges(start = c(10, 60), end = c(30, 80))),
#'   sample3 = GRanges(seqnames = "chr1", ranges = IRanges(start = c(15, 55), end = c(25, 75)))
#' )
#'
#' # Find common peaks
#' common_peaks <- find_common_peaks(gr_list)
#'
#' # View results
#' print(common_peaks)
#'
#' @export
find_common_peaks <- function(gr_list) {
  # Validate input
  if (!is.list(gr_list) || !all(sapply(gr_list, inherits, "GRanges"))) {
    stop("Input must be a list of GRanges objects.")
  }
  
  # Start with the first GRanges object
  common_peaks <- gr_list[[1]]
  
  # Iteratively find overlaps across all GRanges objects
  for (i in 2:length(gr_list)) {
    current_gr <- gr_list[[i]]
    
    # Find overlaps
    overlaps <- findOverlaps(common_peaks, current_gr)
    
    # Subset to overlapping regions
    common_peaks <- subsetByOverlaps(common_peaks, current_gr)
  }
  
  # Assign custom names to the common peaks
  mcols(common_peaks)$name <- paste0("common_peak_", seq_along(common_peaks))
  
  return(common_peaks)
}

###################################
# FIND UNION PEAKS
###################################

#' Find Union of Peaks Across Multiple ATAC-seq Samples
#'
#' @description 
#' This function identifies the union of all genomic regions (peaks) present in 
#' any of the input ATAC-seq samples. Overlapping peaks across samples are merged 
#' into a single region using `reduce()`.
#'
#' @param gr_list A list of GRanges objects. Each GRanges represents ATAC-seq peaks
#'                from a single sample. The list should be named for clarity.
#' 
#' @return A GRanges object representing the merged union of all peak regions across
#'         all input samples. Each region is labeled with a unique name in the metadata
#'         column `name`, formatted as "union_peak_<number>".
#'
#' @details 
#' - The function validates input as a list of GRanges objects.
#' - All peaks from all samples are concatenated and reduced (i.e., merged).
#' - Useful for creating a consensus peak set that represents any signal across samples.
#'
#' @examples
#' gr_list <- list(
#'   sample1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 50), c(20, 70))),
#'   sample2 = GRanges(seqnames = "chr1", ranges = IRanges(c(15, 60), c(35, 80)))
#' )
#' union_peaks <- find_union_peaks(gr_list)
#' print(union_peaks)
#'
#' @export
# Define function to find union peaks using a for loop
find_union_peaks <- function(gr_list) {
  # Initialize combined GRanges object
  combined_gr <- GRanges()
  
  # Loop through the gr_list and combine each GRanges element
  for (gr in gr_list) {
    # Ensure each element is a GRanges object
    if (inherits(gr, "GRanges")) {
      combined_gr <- c(combined_gr, gr)
    } else {
      warning("One of the elements in the list is not a GRanges object.")
    }
  }
  
  # Reduce to get union of overlapping or adjacent peaks
  union_peaks <- reduce(combined_gr)
  
  # Assign unique names for each peak in the result
  mcols(union_peaks)$name <- paste0("union_peak_", seq_along(union_peaks))
  
  return(union_peaks)
}

###################################
# FIND My PEAKS
###################################


#' Find Unique Peaks in a GRanges Object
#'
#' This function identifies peaks that are unique to one GRanges object 
#' (`original_peaks`) by excluding peaks that overlap with another 
#' GRanges object (`common_peaks`). It is designed for genomic peak 
#' analyses, such as identifying condition-specific peaks in comparative 
#' datasets.
#'
#' @param common_peaks A GRanges object representing the common peaks 
#' across multiple conditions or datasets.
#' @param original_peaks A GRanges object representing the peaks for a 
#' specific condition or dataset.
#' 
#' @return A GRanges object containing peaks that are unique to 
#' `original_peaks`, i.e., those that do not overlap with any peaks 
#' in `common_peaks`.
#'
#' @examples
#' # Load required package
#' library(GenomicRanges)
#'
#' # Example GRanges objects
#' common_peaks <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 50, 100), end = c(10, 60, 110)))
#' original_peaks <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(5, 70, 120), end = c(15, 80, 130)))
#'
#' # Identify unique peaks
#' unique_peaks <- find_my_peaks(common_peaks, original_peaks)
#'
#' # Inspect results
#' unique_peaks
#'
#' @export
find_my_peaks <- function(common_peaks, original_peaks) {
  # Find overlaps
  overlaps <- findOverlaps(original_peaks, common_peaks)
  
  # Identify peaks in original_peaks that are not in common_peaks
  unique_peaks <- original_peaks[-queryHits(overlaps)]
  
  # Return unique peaks
  return(unique_peaks)
}

