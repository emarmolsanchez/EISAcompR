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

makeEISAgtfs(annotFile="PATH_to_GTF", boundaryFix=10)
















-------------------------------------------------------------------------------------------------------------------------------

[GenomicFeatures]:https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
[IRanges]:https://bioconductor.org/packages/release/bioc/html/IRanges.html
[Rsubread]:https://bioconductor.org/packages/release/bioc/html/Rsubread.html
[edgeR]:https://bioconductor.org/packages/release/bioc/html/edgeR.html


[1]:https://www.nature.com/articles/nbt.3269
[2]:https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118
[3]:https://academic.oup.com/nar/article/47/8/e47/5345150
[4]:https://academic.oup.com/bioinformatics/article/26/1/139/182458
