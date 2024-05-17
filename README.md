# EISAcompR

EISAcompR is a comprehensive and user-friendly pipeline for performing Exon/Intron Split Analysis (EISA) using RNA-Seq data, as reported by Gaidatzis *et al.* (2015) [[1]]. This document is intended to give a technical supplementary description about how to run the EISAcompR pipeline through a detailed explanation of all the modules that form part of this tool.

We have implemented a framework to prioritize miRNA-driven post-transcriptional regulatory signals jointly with enrichment and covariation analyses to infer which miRNAs might be putatively affecting the post-transcriptional downregulation of targeted mRNAs highlighted by EISA. For further details and citation purposes, please refer to https://onlinelibrary.wiley.com/doi/10.1111/age.13238

&nbsp;

## Contact

emilio.marmol.sanchez@gmail.com

&nbsp;
&nbsp;

-------------------------------------------------------------------------------------------------------------------------------

# Index

- [Prerequisites](https://github.com/emarmolsanchez/EISACompR/#prerequisites)
- [Functions](https://github.com/emarmolsanchez/EISACompR/#functions)
    - [makeEISAgtfs](https://github.com/emarmolsanchez/EISACompR/#makeeisagtfs)
    - [writeEISAgtfs](https://github.com/emarmolsanchez/EISACompR/#writeeisagtfs)
    - [getEISAcounts](https://github.com/emarmolsanchez/EISACompR/#geteisacounts)
    - [getEISAcompR](https://github.com/emarmolsanchez/EISACompR/#geteisacompr)
    - [getCES](https://github.com/emarmolsanchez/EISACompR/#getces)

&nbsp;
&nbsp;

# Prerequisites

The following R libraries are required for running the EISACompR pipeline:

+ [GenomicFeatures] [[2]]
+ [IRanges] [[2]]
+ [Rsubread] [[3]]
+ [edgeR] [[4]]
+ [PCIT] [[5]]
+ [psych]
+ [reshape2]


&nbsp;
&nbsp;

# Install 

```r
devtools::install_github(repo="emarmolsanchez/EISAcompR") 

```

&nbsp;
&nbsp;


# Functions

The next R functions are implemented for running EISA on a set of genes from RNA-Seq data or any other quantification method able to discriminate between exonic and intronic genomic ranges. All required steps are thoroughly explained and exemplified as follows:

&nbsp;

## makeEISAgtfs

The first EISACompR function aims to build custom annotation files in GTF format for late reads quantification spanning either exoninc or intronic genomic ranges. Reference GTF annotation for your species of interest is needed in order to generate the corresponding exonic/intronic regions for running EISA analysis.

&nbsp;

**This function requires 2 arguments:**

+ PATH to reference GTF annotation file.
+ Exon junction boundary correction threshold (10 bp by default).

&nbsp;

**Example of usage:**

```r

GTFs <- makeEISAgtfs(annotFile="PATH_to_GTF", boundaryFix=10)

```
&nbsp;

Once the function has run, it will generate both exonic and intronic ad hoc GTFs generated after splitting the reference GTF provided. Intronic ranges overlapping any other exonic regions from alternative splicing isoforms or other annotated genes spanning the same region will be removed. Additionally, exonic ranges will be enlarged by 10 bp (or any other bp threshold established by the user), in order to avoid counting reads mapping to exon/intron junctions as intronic.

Most commonly, the majority of quantifier tools take only the exonic fraction of each feature (genes/transcripts). In this way, intronic regions are ignored when obtaining gene abundance estimations using canonical GTF annotation files. Users might also prefer to estimate intronic ranges based on considering the whole gene-body of each feature and the exonic fraction separately. Then, the exonic fraction can be subtracted to the whole gene-body quantification for each feature. If users would like to test such alternative approach, we recommend to use the gene-body function from FeatureCounts tool [[3]] and its exonic mode (default). The intronic estimates will be the result of subtracting the exonic counts matrix from the gene-body counts matrix.

&nbsp;
&nbsp;
&nbsp;

## writeEISAgtfs

Users can export the exon/intron annotation files in GTF format by using the writeEISAgtfs function.

&nbsp;

**This function requires two arguments:**

+ Object previously stored with exon/intron split GTFs generated by `makeEISAgtfs`.
+ PATH to desired output.

&nbsp;

**Example of usage:**

```r

writeEISAgtfs(GTFs, "~/")

```
&nbsp;

This command will create two files called `exons.gtf` and `introns.gtf` at the predefined `PATH`, which can be used as canonical GTF annotation for gene quantification using any preferred quantifier tool able to process GTF files as input annotation format (tested on FeatureCounts tool [[3]])

&nbsp;
&nbsp;
&nbsp;

## getEISAcounts

This function allows users to quantify the level of expression of their aligned reads by making use of custom-created exonic/intronic GTF annotation files. The featureCounts built-in function from Rsubread package [[3]] is used for such purpose. Users may optionally decide to use any other preferred quantification pipeline provided they use the exonic/intronic GTF annotation files previously generated as a reference. Then, users can input their count matrices in further steps.

&nbsp;

**This function requires five arguments:**

+ A character vector giving PATH to input BAM/SAM files containing read mapping results.
+ PATH to exonic/intronic GTF annotation file.
+ Strandness (either 0 = unstranded, 1 = stranded forward or 2 = stranded reversed).
+ Number of available threads for parallel computation.
+ Boolean specifying if input sequencing data is of type Paired-end (TRUE/FALSE).

&nbsp;

**Example of usage:**

```r

exon_counts <- getEISAcounts(files=vector_of_input_BAM/SAM, annotFile="PATH_to_exon_GTF", strandness=0/1/2, nthreads=1, PairedEnd=TRUE)
intron_counts <- getEISAcounts(files=vector_of_input_BAM/SAM, annotFile="PATH_to_intron_GTF", strandness=0/1/2, nthreads=1, PairedEnd=TRUE)

```
&nbsp;

Once the function has run, it will create a raw count matrix stored at `counts` object with quantification estimates of each mapped read to the corresponding exonic/intronic regions and gene assignment.

&nbsp;
&nbsp;
&nbsp;

## getEISAcompR

This is the main function used to calculate exon/intron split estimates and to compute the transcriptional and post-transcriptional components of each gene, as well as to infer the significance of each regulatory component independently. The pipeline is implemented for two-group contrast (tipically control vs treated) and includes the possibility of user-defined batch correction. We, however, recommend to remove any possible sources of batch bias before inputing the expression matrices to the function. Users must prepare their exon/intron counts, design and batch (optional) matrices including data for same number of selected samples. In the event that design matrix includes an additional batch effect apart from sample name (1st) and group assignment (2nd), it will be considered as independent batch effect. If no batch correction needs to be implemented, please only use a two-column design matrix. 

While the likelihood ratio test (LRT) is a more obvious choice for obtaining EISA and DE significance estimates, the quasi-likelihood F-test (QLF) is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It provides more robust and reliable error rate control when the number of replicates is small. If the experimental design includes a good replication number (typically more than 10 replicates per group), users can opt for the more flexible LRT approach. If the design is limited in replication, we recommend to implement the QLF method.

&nbsp;

**This function requires seven arguments:**

+ Exonic raw counts (genes in rows and samples in columns).
+ Intronic raw counts (genes in rows and samples in columns).
+ Design matrix (1st column = sample names; 2nd column = group assignment + optionally one additional column with batch effect).
+ Boolean to perform filtering based on expression criteria to remove lowly expressed genes (TRUE/FALSE).
+ Type of test to perform for inferring significance in EISA and DE results. Quasi-likelihood F-test ("QLF") or likelihood ratio test ("LRT"). Default is QLF.
+ Percentage of samples showing minimum expression threshold for filtering (50% by default).
+ counts-per-million (CPM) expression threshold for filtering lowly expressed genes (1 CPM by default).

&nbsp;

**Example of usage:**

```r

eisa <- getEISAcompR(Exons=exon_counts$counts, Introns=intron_counts$counts, design=design_matrix, filterExpr=TRUE, model="QLF", percent=0.5, cpm=1)

```
&nbsp;

Once the function has run, a list of tabulated results will be created containing the following four tables:

+ resDE = Differential Expression analysis for exonic counts using glmQLFTest from edgeR [[4]] package.
+ resTc = EISA for the transcriptional component effect on each analyzed gene.
+ resPTc = EISA for the post-transcriptional component effect on each analyzed gene.
+ Expr_Int = Normalized log2 expression matrix for Intronic counts.
+ Expr_Ex = Normalized log2 expression matrix for Exonic counts.

&nbsp;

For transcriptional (Tc) and post-transcriptional (PTc) components, generally, the higher their absolute values (either showing negative or positive regulatory influence), the more relevant regulatory effects could be inferred. Please be aware that the abscence of significance in PTc scores might indicate the presence of mixed transcriptional and post-transcriptional componentes affecting the same gene. Any PTc component showing significant interaction with any other transcriptional (Tc) influence detected at intronic levels will be shown as not significant. Only PTc components with opposing direction to Tc components or when no Tc interacting component is detected will appear as significant.

Users should compare canonical differential expression (DE) results and significance for their genes of interest, as well as their Tc and PTc components, in order to extract meaningfull information about the putative regulatory influence affecting their genes of interest.

&nbsp;

Output interpretation:

+ Significant PTc + Significant DE = Gene showing Post-transcriptional regulatory signal
+ Significant Tc + Significant DE = Gene showing Transcriptional regulatory signal
+ Non-significant PTc + Significant Tc + Significant DE = Gene showing mixed Transcriptional and Post-transcriptional regulatory signal

&nbsp;
&nbsp;
&nbsp;

## getCES

This function is intented to calculate the Covariation Enrichment Score (CES) value for each gene in a given set of genes compared to a defined background. The CES value reflects the frequency with which the mRNA expression correlation between members of the defined gene set of interest (in this case, genes showing high post-transcriptional signals and downregulation driven by miRNAs) is different compared with other genes defined as background (in this case, differentially expressed genes with no evidence of miRNA-related downregulatory effects). Other set of genes might be used as background, such as, for instance, the whole set of expressed genes excluding those included in the gene set of interest. However, including a high number of genes as background might overestimate CES values and result in high memory usage and prolonged running time. We recommend to limit the analysis to DE genes and a reduced number of genes of interest with high post-transcriptional signals.

In this way, we can represent the variability in covariation within a set of genes as a fold chang comparing the increase or reduction of significant covariation events relative to the overall gene expression background. Genes sets showing a coordinated high post-transcriptional downregulation are expected to show an increased covariation among each other when compared with their covariation with other genes.

To achive this purpose, we have implemented a network-oriented filtering criteria based on Partial Correlations and Information Theory ([PCIT]) approach as proposed by Reverter *et al.* (2008) [[5]]. By using first-order partial correlation coefficients estimated for each trio of genes along with an information theory approach, this tool identifies meaningful gene-to-gene nteractions. This approach aims to determine truly informative correlations between node pairs (genes in our context), once the influence of other nodes in the network has been considered. Alternative methods based on naive *P*-value and multiple testing corrected *P*-value with the False Discovery Rate (FDR) method [[6]] are also provided as alternative to the [PCIT] algorithm.

&nbsp;

**This function requires five arguments**

+ Normalized counts in log2 scale belonging to the gene set of interest (genes in rows and samples in columns).
+ Normalized counts in log2 scale belonging to differentially expresssed genes excluding genes present in the set of interest (genes in rows and samples in columns).
+ Method to calculate pairwise correlation among genes (available methods are "pearson", "spearman" and "kendall". Spearman by default). 
+ Correlation threshold for prioritizing significant pairwise covariation events (0.6 by default).
+ Network inference algorithm to prioritize significant covariation events (available methods are "pcit", "pvalue", "fdr". PCIT by default).

&nbsp;

**Example of usage:**

```r

CES <- getCES(ExprSet=Gene_set_counts, ExprBack=DE_genes_counts, method="spearman", cor=0.6, Filter="pcit")

```
&nbsp;

Once the function has run, it will calculate a CES value for each gene included in the gene set of interest `ExprSet`, as a measure of the fold change of significant covariation events when comparing them with each other within the gene set of interest and with the whole set of genes used as background.

&nbsp;
&nbsp;
&nbsp;

-------------------------------------------------------------------------------------------------------------------------------

# References

1. [Gaidatzis D et al. (2015) Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. *Nature Biotechnology*, 33, 722–729.]

2. [Lawrence M et al. (2013) Software for computing and annotating genomic ranges. *PLoS Computational Biology*, 9, e1003118.]

3. [Liao Y et al. (2019) The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleid Acids Research*, 47, e47.]

4. [Robinson MD et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26, 139–140.]

5. [Reverter A et al. (2008) Combining partial correlation and an information theory approach to the reversed engineering of gene co-expression networks. *Bioinformatics*, 24, 2491-97.]

6. [Benjamini Y & Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B (Methodological)*, 57, 289-300.]

[GenomicFeatures]:https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
[IRanges]:https://bioconductor.org/packages/release/bioc/html/IRanges.html
[Rsubread]:https://bioconductor.org/packages/release/bioc/html/Rsubread.html
[edgeR]:https://bioconductor.org/packages/release/bioc/html/edgeR.html
[PCIT]:https://github.com/nathanhaigh/pcit
[psych]:https://cran.r-project.org/web/packages/psych
[reshape2]:https://cran.r-project.org/web/packages/reshape2


[1]:https://www.nature.com/articles/nbt.3269
[2]:https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118
[3]:https://academic.oup.com/nar/article/47/8/e47/5345150
[4]:https://academic.oup.com/bioinformatics/article/26/1/139/182458
[5]:https://academic.oup.com/bioinformatics/article/24/21/2491/192682
[6]:https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1995.tb02031.x

[Gaidatzis D et al. (2015) Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. *Nature Biotechnology*, 33, 722–729.]:https://www.nature.com/articles/nbt.3269
[Lawrence M et al. (2013) Software for Computing and Annotating Genomic Ranges. *PLoS Computational Biology*, 9, e1003118.]:https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118
[Liao Y et al. (2019) The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleid Acids Research*, 47, e47.]:https://academic.oup.com/nar/article/47/8/e47/5345150
[Robinson MD et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26, 139–140.]:https://academic.oup.com/bioinformatics/article/26/1/139/182458
[Reverter A et al. (2008) Combining partial correlation and an information theory approach to the reversed engineering of gene co-expression networks. *Bioinformatics*, 24, 2491-97.]:https://academic.oup.com/bioinformatics/article/24/21/2491/192682
[Benjamini Y & Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B (Methodological)*, 57, 289-300.]:https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1995.tb02031.x

