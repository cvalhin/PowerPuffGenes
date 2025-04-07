# PowerPuffGenes: RNA-seq analysis of dox induction

## Introduction
Understanding the temporal dynamics of gene expression is crucial for deciphering complex biological processes like stem cell differentiation or response to stimuli. This project utilizes a doxycycline (dox)-inducible system in mouse stem cells, where dox triggers gene expression via the rtTA-TRE mechanism. While powerful, concerns exist that dox itself might exert off-target effects on the transcriptome. Leveraging Dr. Rinn's lab RNA-seq dataset spanning five time points (0, 12, 24, 48, 96 hours) with three replicates. We employ the nf-core/rnaseq pipeline for data processing and DESeq2 for differential expression analysis, providing a practical application of bioinformatics skills to time-course data exploration.

## Goals
1) Identify genes significantly regulated following dox administration.
2) Characterize the temporal profiles of these regulated genes, distinguishing between transient and sustained/prolonged responses.
3) Investigate the role, significance, and chromosomal location of these genes.

## Results
After performing DESeq2, we found that from the 55401 total genes, 794 were significantly expressed genes. From our significantly expressed genes, we generated a list of prolonged genes that represent candidates that were up or down regulated after Dox treatment at an early time point (12 or 24 hours) and maintained their trend until the 96 hour time point (last collected). A total of 27 prolonged genes were found. We looked at prolonged up or down regulated genes, because if their levels are high or low throughout all time points, they might be causing a significant change in important genes.


We explored prolonged genes roles in the cell by looking at their gene ontologies (GO), and we found three genes with strong evidence of important regulatory roles in mice, H19, Usp26, and Ier3. From those three genes H19 had prolonged up regulation, while Ier3 and Usp26 had prolonged down regulation. Their prolonged change in expression suggests that these genes might be changing permanently, compromising cell integrity by changing regulatory pathways. However, additional time points should be sampled to test this hypothesis.


### H19 

This gene is a long non-coding RNA (lncRNA) gene. Strong evidence suggests that this gene is involved in the control of embryonic growth. In addition, it has been found to regulate nine genes of an imprinted gene network [@monnier_h19_2013]. Its prolonged up regulation after Dox treatment suggests a potential role in maintaining cellular processes related to growth and development, possibly by modulating the activity of imprinted genes that are critical for cell differentiation and tissue patterning. 

![H19 gene track](results/IGV/h19_gene_track.png)
**IGV track with peaks of H19 gene.** Peaks represent the coverage of sequencing reads across this gene. Note that peaks are low at 0 and 12 hours, but start getting higher at 24. By the time they have reached 48 hours, they are high and they keep this high trend until the last sampled time point of 96 hours. 

### Ier3

The immediate early response [(Ier3)](https://www.uniprot.org/uniprotkb/P46694/entry) gene may play a role in the ERK pathway for multiple cellular processes, such as growth,  proliferation, differentiation, and survival. This role suggests relevance in gene regulation, making it an interesting and possibly important gene to further explore the effects of Dox treatment. Dox might be affecting the normal function of this gene by down regulating its expression at least 96 hours after treatment. 

![Ier3 gene track](results/IGV/Ier3_IGV_track.png)
**IGV track with peaks of Ier3 gene.** Peaks represent the coverage of sequencing reads across this gene. Note that peaks are higher in the beginning at 0 hours, and they decrease as time increases. They remain low until the last time point sampled of 96 hours. 


### Usp26

The [Usp26](https://www.uniprot.org/uniprotkb/Q99MX1/entry) gene has been involved in deubiquitination pathways, making it potentially relevant for gene regulation. There are many transcription factors in the cell, such as  polycomb repressive complex 1 (PRC1), that are in charge of maintaining and regulating epigenetic marks in the genome [@tamburri_histone_2020]. Many of these transcription factors, need ubiquitination for their function, making a deubiquitination crucial for this process as well [@mark_ubiquitindependent_2021]. For example, this gene has been found to be involved in somatic cell reprogramming through the K48 deubiquitination of two protein components of PRC1.


![Usp26 gene track](results/IGV/Usp26_IGV_track2.png)
**IGV track with peaks of Usp26 gene.** Peaks represent the coverage of sequencing reads across this gene. Note that peaks are higher at 0 and 12 hours, but are lower when 24 hours are reached. They keep beung low until the last time point of 96 hours.

- You can look at the first draft for the project by :
  - Opening the main_project.pdf
  - Knitting the main_project.Rmd 
  - Opening the main_project.html in web browser