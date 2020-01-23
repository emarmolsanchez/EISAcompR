# EISACompR

EISACompR is a comprehensive and user-friendly pipeline for performing Exon/Intron Split Analysis (EISA) using RNA-Seq data, as reported by Gaidatzis *et al.* (2015) [[1]]. This document is intended to give a technical supplementary description about how to run the EISACompR pipeline through a detailed explanation of all the modules that form part of this tool.

&nbsp;
&nbsp;

-------------------------------------------------------------------------------------------------------------------------------

# Index

- [Prerequisites](https://github.com/emarmolsanchez/EISACompR/#prerequisites)
- [Functions](https://github.com/emarmolsanchez/EISACompR/#functions)
    - [makeEISAgtfs](https://github.com/emarmolsanchez/EISACompR/#makeeisagtfs)
    - [writeEISAgtfs](https://github.com/emarmolsanchez/EISACompR/#writeeisagtfs)
    - [getEISAcounts](https://github.com/emarmolsanchez/EISACompR/#geteisacounts)
    - [getEISAcomp](https://github.com/emarmolsanchez/EISACompR/#geteisacomp)

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

Once the function has run, an EISACompR object will be stored at predefined object `GTFs`, storing both exonic and intronic GTFs generated after splitting the reference GTF provided. Intronic ranges overlapping any other exonic regions from alternative splicing isoforms or other annotated genes spanning the same region will be removed. Additionally, exonic ranges will be enlarged by 10 bp (or any other bp threshold established by the user), in order to avoid counting reads mapping to exon/intron junctions as intronic.

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

counts <- getEISAcounts(files=vector_of_input_BAM/SAM, annotFile="PATH_to_exon/intron_GTFs", strandness=0/1/2, nthreads=1, PairedEnd=TRUE)

```

Once the function has run, it will create a raw count matrix stored at `counts` object with quantification estimates of each mapped read to the corresponding exonic/intronic regions and gene assignment.

&nbsp;
&nbsp;

## getEISAcomp













-------------------------------------------------------------------------------------------------------------------------------

[GenomicFeatures]:https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
[IRanges]:https://bioconductor.org/packages/release/bioc/html/IRanges.html
[Rsubread]:https://bioconductor.org/packages/release/bioc/html/Rsubread.html
[edgeR]:https://bioconductor.org/packages/release/bioc/html/edgeR.html


[1]:https://www.nature.com/articles/nbt.3269
[2]:https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118
[3]:https://academic.oup.com/nar/article/47/8/e47/5345150
[4]:https://academic.oup.com/bioinformatics/article/26/1/139/182458
