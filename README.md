# PowerPuffGenes

**RNA-seq and ATAC-seq Time Course Analysis of Doxycycline-Induced Gene Expression in Mouse and Human Stem Cells**

## Overview

This project explores the temporal dynamics of gene expression and chromatin accessibility in response to doxycycline (dox) exposure using RNA-seq and ATAC-seq in mouse and human stem cells. It leverages data from Dr. Rinn’s lab, focusing on both transcriptional and epigenetic changes across multiple time points.

## Goals

- Identify genes significantly regulated by dox administration.
- Classify genes into transient and prolonged responders.
- Analyze chromosomal distribution of expression changes.
- Integrate RNA-seq and ATAC-seq data to understand regulatory mechanisms.
- Highlight key genes of biological interest.

## Key Results

- PCA analysis revealed progressive transcriptomic shifts over 96 hours of dox exposure.
- Volcano and circos plots showed widespread, time-dependent differential expression.
- A subset of genes displayed sustained expression changes across the full time course.
- ATAC-seq analysis indicated chromatin accessibility changes related to prolonged gene regulation.
- Important genes such as **Usp26**, **Cbr3**, **MAG**, **AUTS2**, and **NOD1** were identified with significant changes in expression and/or chromatin accessibility.
- Full reference list available in [Results.md](Results.md).

---

## Tools and Pipelines Used

- **RNA-seq pipeline:** [nf-core/rnaseq v3.18](https://nf-co.re/rnaseq)
- **ATAC-seq pipeline:** [nf-core/atacseq v2.1.2](https://nf-co.re/atacseq/2.1.2)
- **Statistical analysis:** DESeq2
- **Visualization:** PCA plots, volcano plots, circos plots, heatmaps, IGV tracks

## How to Reproduce

1. Process RNA-seq FASTQ files with `nf-core/rnaseq`.
2. Process ATAC-seq FASTQ files with `nf-core/atacseq`.
3. Follow RNA-seq and ATAC-seq analysis workflow in [script_mouse.Rmd](script_mouse.Rmd) or [script_human.Rmd](script_human.Rmd)
4. Visualize specific genes of interest using IGV.

---

## License

This project is released under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

Special thanks to Dr. Rinn’s lab for providing the dataset and to the [Computational Genomic Lab](https://www.lncrna.io/teaching) for guidance on RNA-seq analysis and interpretation.


