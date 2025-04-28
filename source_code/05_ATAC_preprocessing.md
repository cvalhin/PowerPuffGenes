PowerPuffGenes Report
================
2025-04-16

#### Objective: To determine if chromatin accessibility changes due to dox exposure in human stem cells.

We have performed a time course series of experiements measuring
chromating accessibility (ATACseq peaks)upon exposure to dox. These are
0, 30, 60, 90, 120, 150 minutes. Each time point has a replicate. The
fastq sequencing files were processed by the NF_CORE ATACseq pipeline v-
2.12 (<https://nf-co.re/atacseq/2.1.2/>).

We will be using the MACS2 output files from the NF_CORE Pipeline For
example: -broad.peak -consensus.broad.peak -consensus.feature.counts

#### Approach:

#### Parameter and path definition

``` r
PEAK_PATH         <- "/scratch/Shares/rinnclass/MASTER_CLASS/DATA/human_atacseq"
ANNOTATION_FILE   <- "/scratch/Shares/rinnclass/MASTER_CLASS/DATA/genomes/Homo_sapiens/Gencode/v38/gencode.v38.annotation.gtf"
RESULT_PATH       <- "../result_human"
```

##### (a) Loading in ATACseq peak files with custom function import_peaks (list of GRanges output)

``` r
# creating a file list also needed for import_peaks function to get sample name associated with file
fl <- list.files(PEAK_PATH, full.names = TRUE, pattern = ".broadPeak")

# running import_peaks
my_peaks <- import_peaks(PEAK_PATH)

# finding consensus peaks among replicates
my_consensus_peaks <- find_consensus_peaks(my_peaks)

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

##### (b) finding number of peaks common in all samples using find_common_peaks custom funciton

``` r
# run find_common_peaks function
common_peaks <-  suppressWarnings(find_common_peaks(my_consensus_peaks))
```

Here are the number of peaks that overlap in all samples: 58068

##### (c) Peaks that are unique to dox and non-dox conditions

Using find_common_overlaps peaks that are specific to dox or non-dox
will be identified.

Non-dox samples

``` r
# common peaks in non-dox (0 time point)
# TODO: What does KO mean?
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

print(c("This is how many peaks are common in dox",length(dox_common_peaks)))
```

    ## [1] "This is how many peaks are common in dox"
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

##### (d) Creating mouse gene, lincrna, mRNA annotation GRange objects

``` r
# Loading gencode genome annotation as GRanges 
gencode_gr <- rtracklayer::import(ANNOTATION_FILE)

# all genes
gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
gene_promoters <- promoters(gencode_genes, upstream = 2000, downstream = 2000)
```

##### (e) Compare overlaps of dox and non-dox peaks with gene annotaitons

Now we will overlap our dox and non-dox unique peaks with genome
annotations (gene promoters) First we will find number of overlaps with
gene promoters and then genes that had changed in RNAseq

``` r
# gr_list of promoters and peaks common to non_dox condition
gr_list_gene_promoter_non_dox_ov <- list( gene_promoters = gene_promoters, non_dox_peaks = non_dox_common_peaks)
non_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_non_dox_ov))

print("This is how many non-dox_unique peaks overlapped gene promoters")
```

    ## [1] "This is how many non-dox_unique peaks overlapped gene promoters"

``` r
length(non_dox_gene_promoter_ov)
```

    ## [1] 26274

``` r
# peaks common to dox condition overlapped with gene promoters
gr_list_gene_promoter_dox_ov <- list( gene_promoters = gene_promoters, dox_peaks = dox_common_peaks)
dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_dox_ov))

print(c("This is how many dox_unique peaks overlapped gene promoters", length(dox_gene_promoter_ov)))
```

    ## [1] "This is how many dox_unique peaks overlapped gene promoters"
    ## [2] "22896"

``` r
# peaks common in dox and non-dox conditions overlapped with gene promoters
gr_list_gene_promoter_dox_non_dox_ov <- list( gene_promoters = gene_promoters, dox_peaks = dox_non_dox_ov)
dox_non_dox_gene_promoter_ov <- suppressWarnings(find_common_peaks(gr_list_gene_promoter_dox_non_dox_ov))

print(c("This is how many shared (common) dox, non-dox peaks overlapped gene promoters", length(dox_non_dox_gene_promoter_ov)))
```

    ## [1] "This is how many shared (common) dox, non-dox peaks overlapped gene promoters"
    ## [2] "22878"

``` r
save(dox_non_dox_gene_promoter_ov, dox_gene_promoter_ov, non_dox_gene_promoter_ov,
     file=file.path(RESULT_PATH, "ATAC_gene_promoter.RData"))
```
