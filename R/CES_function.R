#'@title getEISAcompR
#'@description
#'This function is intented to calculate the CES value for each gene
#'in a given set of genes compared to a defined background.
#'The CES value reflects the frequency with which the mRNA expression
#'correlation between members of the defined gene set of interest
#'-in this case, genes showing high post-transcriptional signals
#'and downregulation driven by miRNAs- is different compared with
#'other genes defined as background -in this case, differentially
#'expressed genes with no evidence of miRNA-related downregulatory
#'effects-. Other set of genes might be used as background, such as,
#'for instance, the whole set of expressed genes excluding those
#'included in the gene set of interest. However, including a high
#'number of genes as background might overestimate CES values and
#'result in high memory usage and prolonged running time.
#'We recommend to limit the analysis to DE genes and a reduced
#'number of genes of interest with high post-transcriptional signals.

#'In this way, we can represent the variability in
#'covariation within a set of genes as a fold chang comparing the
#'increase or reduction of significant covariation events relative
#' to the overall gene expression background. Genes sets showing a
#' coordinated high post-transcriptional downregulation are expected
#' to show an increased covariation among each other when compared
#' with their covariation with other genes.

#'To achive this purpose, we have implemented a network-oriented
#'filtering criteria based on Partial Correlations and Information
#'Theory approach as proposed by Reverter et al. (2008).
#'By using first-order partial correlation coefficients estimated
#'for each trio of genes along with an information theory approach,
#'this tool identifies meaningful gene-to-gene nteractions.
#'This approach aims to determine truly informative correlations
#'between node pairs -genes in our context-, once the influence
#'of other nodes in the network has been considered.
#'Alternative methods based on naive P-value and multiple testing
#'corrected P-value with the False Discovery Rate method
#'are also provided as alternative to the PCIT algorithm.
#'
#'

#'@param ExprSet Normalized counts in log2 scale belonging to the gene set of interest (genes in rows and samples in columns).
#'@param ExprBack Normalized counts in log2 scale belonging to differentially expresssed genes excluding genes present in the set of interest (genes in rows and samples in columns).
#'@param method Method to calculate pairwise correlation among genes (available methods are "pearson", "spearman" and "kendall". Spearman by default).
#'@param cor  Correlation threshold for prioritizing significant pairwise covariation events (0.6 by default).
#'@param Filter Network inference algorithm to prioritize significant covariation events (available methods are "pcit", "pvalue", "fdr". PCIT by default).
#'@import PCIT
#'@import reshape2
#'@import psych
#'
#'@return {Once the function has run, it will calculate a CES value
#' for each gene included in the gene set of interest ExprSet, as a
#' measure of the fold change of significant covariation events when
#'  comparing them with each other within the gene set of interest and
#'  with the whole set of genes used as background.}
#'

#'@export getCES
#'@references
#' \enumerate{
#'  \item Gaidatzis D et al. (2015) Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. Nature Biotechnology, 33, 722–729.
#'  \item Lawrence M et al. (2013) Software for computing and annotating genomic ranges. PLoS Computational Biology, 9, e1003118.
#'  \item Liao Y et al. (2019) The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleid Acids Research, 47, e47.
#'  \item Robinson MD et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139–140.
#'  \item Reverter A et al. (2008) Combining partial correlation and an information theory approach to the reversed engineering of gene co-expression networks. Bioinformatics, 24, 2491-97.
#'  \item Benjamini Y & Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B (Methodological), 57, 289-300.

#'}
#'
#'@author Emilio Mármol Sánchez
#'@examples
#' {
#' \dontrun{
#' CES <- getCES(ExprSet=Gene_set_counts, ExprBack=DE_genes_counts, method="spearman", cor=0.6, Filter="pcit")
#' }
#'
#' }
#'@name getCES
#'@rdname getCES-getCES



getCES <- function(ExprSet, ExprBack, method="spearman", cor=0.6, Filter="pcit"){

  Expr <- rbind(ExprSet, ExprBack)

  if(Filter == "pcit"){

    corr_expr <- cor(t(Expr), method=method)
    pcit_expr <- pcit(corr_expr)

    #Get occurrence tables
    signif <- idx(pcit_expr)
    nonsignif <- idxInvert(nrow(corr_expr), signif)
    corr_expr[nonsignif] <- 0
    corr_expr[corr_expr < -cor] <- 1
    corr_expr[corr_expr == 1] <- 0
    corr_expr[corr_expr > cor] <- 1
    corr_expr[corr_expr < 1] <- 0
    corr_expr <- as.data.frame(corr_expr)

   } else if(Filter == "pvalue"){

     corr_expr <- corr.test(t(Expr), method=method, adjust="none", ci=FALSE)
     corr_pvals <- corr_expr$p
     corr_pvals <- melt(corr_pvals)
     corr_expr <- corr_expr$r

     #Get occurrence tables
     signif = which(corr_pvals$value < 0.05)
     nonsignif = idxInvert(nrow(corr_expr), signif)
     corr_expr[nonsignif] <- 0
     corr_expr[corr_expr < -cor] <- 1
     corr_expr[corr_expr == 1] <- 0
     corr_expr[corr_expr > cor] <- 1
     corr_expr[corr_expr < 1] <- 0
     corr_expr <- as.data.frame(corr_expr)

   } else if(Filter == "fdr"){

     corr_expr <- corr.test(t(Expr), method=method, adjust="fdr", ci=FALSE)
     corr_pvals <- corr_expr$p.adj
     corr_pvals <- melt(corr_pvals)
     corr_expr <- corr_expr$r

     #Get occurrence tables
     signif = which(corr_pvals$value < 0.05)
     nonsignif = idxInvert(nrow(corr_expr), signif)
     corr_expr[nonsignif] <- 0
     corr_expr[corr_expr < -cor] <- 1
     corr_expr[corr_expr == 1] <- 0
     corr_expr[corr_expr > cor] <- 1
     corr_expr[corr_expr < 1] <- 0
     corr_expr <- as.data.frame(corr_expr)

   }


  #Get sets and CES
  corr_expr_set <- corr_expr[1:nrow(ExprSet),1:nrow(ExprSet)]
  Ratios_expr_set <- rowSums(corr_expr_set)/nrow(ExprSet)

  corr_expr_overall <- corr_expr[1:nrow(ExprSet),]
  Ratios_expr_overall <- rowSums(corr_expr_overall)/nrow(Expr)

  CES <- Ratios_expr_set / Ratios_expr_overall


  return(CES)

}

