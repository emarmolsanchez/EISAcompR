# EISAcompR

EISAcompR is a comprehensive and user-friendly pipeline for performing Exon/Intron Split Analysis (EISA) using RNA-Seq data, as reported by Gaidatzis *et al.* (2015) [[1]]. This document is intended to give a technical supplementary description about how to run the EISAcompR pipeline through a detailed explanation of all the modules that form part of this tool.

&nbsp;
&nbsp;

-------------------------------------------------------------------------------------------------------------------------------

# Index

- [Prerequisites](https://github.com/emarmolsanchez/EISACompR/#prerequisites)
- [Functions](https://github.com/emarmolsanchez/EISACompR/#functions)
    - [makeEISAgtfs](https://github.com/emarmolsanchez/EISACompR/#makeeisagtfs)
    - [writeEISAgtfs](https://github.com/emarmolsanchez/EISACompR/#writeeisagtfs)
    - [getEISAcounts](https://github.com/emarmolsanchez/EISACompR/#geteisacounts)
    - [getEISAcompR](https://github.com/emarmolsanchez/EISACompR/#geteisacomp)

&nbsp;
&nbsp;

# Prerequisites

The following R libraries are required for running the EISACompR pipeline:

+ [GenomicFeatures] [[2]]
+ [IRanges] [[2]]
+ [Rsubread] [[3]]
+ [edgeR] [[4]]


&nbsp;
&nbsp;

# Install 

devtools::install_github(repo="emarmolsanchez/EISACompR") 

Note: if EISACompR R package is not in the default branch, use instead: 

devtools::install_github(repo="lauzingaretti/EISACompR@R_package") 

devtools::install_github(repo="emarmolsanchez/EISACompR@R_package") 

where @ indicates the branch allocating the R package structure. 


# Functions

The next R functions are implemented for running EISA on a set of genes from RNA-Seq data or any other quantification method able to discriminate between exonic and intronic genomic ranges. All required steps are thoroughly explained and exemplified as follows:

&nbsp;

## makeEISAgtfs

The first EISACompR function aims to build custom annotation files in GTF format for late reads quantification spanning either exoninc or intronic genomic ranges. Reference GTF annotation for your species of interest is needed in order to generate the corresponding exonic/intronic regions for running EISA analysis.

This function requires 2 arguments:

+ PATH to reference GTF annotation file.
+ Exon junction boundary correction threshold (10 bp by default).

Example of usage:

```r

GTFs <- makeEISAgtfs(annotFile="PATH_to_GTF", boundaryFix=10)

```
&nbsp;

Once the function has run, an EISAcompR object will be stored at predefined object `GTFs`, storing both exonic and intronic GTFs generated after splitting the reference GTF provided. Intronic ranges overlapping any other exonic regions from alternative splicing isoforms or other annotated genes spanning the same region will be removed. Additionally, exonic ranges will be enlarged by 10 bp (or any other bp threshold established by the user), in order to avoid counting reads mapping to exon/intron junctions as intronic.

&nbsp;
&nbsp;

## writeEISAgtfs

Users can export the exon/intron annotation files in GTF format by using the ad-hoc writeEISAgtfs function.

This function requires two arguments:

+ EISACompR object previously stored with exon/intron split GTFs.
+ PATH to desired output.

Example of usage:

```r

writeEISAgtfs(GTFs, "~/")

```
&nbsp;

This command will create two files called `exons.gtf` and `introns.gtf` at `HOME` PATH, which can be used as canonical GTF annotation for gene quantification using any preferred quantifier tool able to process GTF files as input annotation format.

&nbsp;
&nbsp;

## getEISAcounts

This function allows users to quantify the level of expression of their aligned reads by making use of custom-created exonic/intronic GTF annotation files. The featureCounts built-in function from Rsubread package is used for such purpose. Users may optionally decide to use any other preferred quantification pipeline provided the use of exonic/intronic GTF annotation as a reference.

This function requires five arguments:

+ A character vector giving PATH to input BAM/SAM files containing read mapping results.
+ PATH to exonic/intronic GTF annotation file.
+ Strandness (either 0 = unstranded, 1 = stranded forward or 2 = stranded reversed).
+ Number of available threads for parallel computation.
+ Boolean specifying if input sequencing data is of type Paired-end (TRUE/FALSE).

Example of usage:

```r

exon_counts <- getEISAcounts(files=vector_of_input_BAM/SAM, annotFile="PATH_to_exon_GTF", strandness=0/1/2, nthreads=1, PairedEnd=TRUE)
intron_counts <- getEISAcounts(files=vector_of_input_BAM/SAM, annotFile="PATH_to_intron_GTF", strandness=0/1/2, nthreads=1, PairedEnd=TRUE)

```
&nbsp;

Once the function has run, it will create a raw count matrix stored at `counts` object with quantification estimates of each mapped read to the corresponding exonic/intronic regions and gene assignment.

&nbsp;
&nbsp;

## getEISAcompR

This is the main function used to calculate exon/intron split estimates and compute transcriptional and post-transcriptional components of each gene, as well as to infer the significance of each regulatory component independently. The pipeline is implemented for two-group contrast (tipically control vs treated) and includes the possibility of user-defined batch correction. Users must prepare their exon/intron counts, design and batch (optional) matrices including data for same number of selected samples. In the event that design matrix includes an additional batch effect apart from sample name (1st) and group assignment (2nd), it will be considered as independent batch effect. If no batch correction needs to be implemented, please only use a two-column design matrix.

This function requires seven arguments:

+ Exonic raw counts (genes in rows and samples in columns).
+ Intronic raw counts (genes in rows and samples in columns).
+ Design matrix (1st = sample names; 2nd = group assignment + optionally one additional column with batch effect).
+ Boolean to perform filtering based on expression criteria to remove lowly expressed genes (TRUE/FALSE).
+ Percentage of samples showing minimum expression threshold for filtering (50% by default).
+ counts-per-million (CPM) expression threshold for filtering lowly expressed genes (1 CPM by default).

Example of usage:

```r

eisa <- getEISAcomp(exons=exon_counts, introns=intron_counts, design=design_matrix, filterExpr=TRUE, percent=0.5, cpm=1)

```
&nbsp;

Once the function has run, an object ot type `EISAcompR` will be created. This object will contain four tables:

+ resDE = Differential Expression analysis for exonic counts using glmQLFTest from [edgeR] package.
+ resTc = EISA for the transcriptional component effect on each analyzed gene.
+ resPTc = EISA for the post-transcriptional component effect on each analyzed gene.
+ Expr_Int = Normalized log2 expression matrix for Intronic counts.
+ Expr_Ex = Normalized log2 expression matrix for Exonic counts.

&nbsp;
&nbsp;

For transcriptional (Tc) and post-transcriptional (PTc) components, generally, the higher their absolute values (either showing negative or positive regulatory influence), the more relevant regulatory effects could be inferred. Please be aware that the abscence of significance in PTc component might indicate the presence of mixed transcriptional and post-transcriptional componentes affecting the same gene. Any PTc component showing significant interaction with any other transcriptional (Tc) influence detected at intronic levels will be shown as not significant. Users should compare canonical differential expression (DE) results and significance for their genes of interest, as well as their Tc and PTc components, in order to extract meaningfull information about the putative regulatory influence affecting their genes of interest.

&nbsp;
&nbsp;

Output interpretation:

+ Significant PTc + Significant DE = Gene showing Post-transcriptional regulatory signal
+ Significant Tc + Significant DE = Gene showing Transcriptional regulatory signal
+ Non-significant PTc + Significant Tc + Significant DE = Gene showing mixed Transcriptional and Post-transcriptional regulatory signal

&nbsp;
&nbsp;

-------------------------------------------------------------------------------------------------------------------------------

# References

1. [Gaidatzis D et al. (2015) Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. *Nature Biotechnology*, 33, 722–729.]

2. [Lawrence M et al. (2013) Software for computing and annotating genomic ranges. *PLoS Computational Biology*, 9, e1003118.]

3. [Liao Y et al. (2019) The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleid Acids Research*, 47, e47.]

4. [Robinson MD et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26, 139–140.]

[GenomicFeatures]:https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
[IRanges]:https://bioconductor.org/packages/release/bioc/html/IRanges.html
[Rsubread]:https://bioconductor.org/packages/release/bioc/html/Rsubread.html
[edgeR]:https://bioconductor.org/packages/release/bioc/html/edgeR.html


[1]:https://www.nature.com/articles/nbt.3269
[2]:https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118
[3]:https://academic.oup.com/nar/article/47/8/e47/5345150
[4]:https://academic.oup.com/bioinformatics/article/26/1/139/182458

[Gaidatzis D et al. (2015) Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. *Nature Biotechnology*, 33, 722–729.]:https://www.nature.com/articles/nbt.3269
[Lawrence M et al. (2013) Software for Computing and Annotating Genomic Ranges. *PLoS Computational Biology*, 9, e1003118.]:https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118
[Liao Y et al. (2019) The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleid Acids Research*, 47, e47.]:https://academic.oup.com/nar/article/47/8/e47/5345150
[Robinson MD et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26, 139–140.]:https://academic.oup.com/bioinformatics/article/26/1/139/182458
