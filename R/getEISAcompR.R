#'@title getEISAcompR
#'@description  This is the main function used to calculate exon/intron split estimates
#'and compute transcriptional and post-transcriptional components of each gene,
#'as well as to infer the significance of each regulatory component independently.
#'The pipeline is implemented for two-group contrast (tipically control vs treated)
#'and includes the possibility of user-defined batch correction.
#'Users must prepare their exon/intron counts, design and batch (optional)
#'matrices including data for same number of selected samples.
#'In the event that design matrix includes an additional batch
#'effect apart from sample name (1st) and group assignment (2nd),
#'it will be considered as independent batch effect.
#'Further details:
#'If no batch correction needs to be implemented, please only use a two-column design matrix.
#'For transcriptional (Tc) and post-transcriptional (PTc) components, generally, the higher their absolute values (either showing negative or positive regulatory influence), the more relevant regulatory effects could be inferred. Please be aware that the abscence of significance in PTc component might indicate the presence of mixed transcriptional and post-transcriptional componentes affecting the same gene. Any PTc component showing significant interaction with any other transcriptional (Tc) influence detected at intronic levels will be shown as not significant. Users should compare canonical differential expression (DE) results and significance for their genes of interest, as well as their Tc and PTc components, in order to extract meaningfull information about the putative regulatory influence affecting their genes of interest.
#'
#'Output interpretation:
#'Significant PTc + Significant DE = Gene showing Post-transcriptional regulatory signal
#'Significant Tc + Significant DE = Gene showing Transcriptional regulatory signal
#'Non-significant PTc + Significant Tc + Significant DE = Gene showing mixed Transcriptional and Post-transcriptional regulatory signal
#'

#'@param exons Exonic raw counts (genes in rows and samples in columns).
#'@param introns Intronic raw counts (genes in rows and samples in columns).
#'@param design Design matrix (1st = sample names; 2nd = group assignment + optionally one additional column with batch effect).
#'@param filterExpr Boolean to perform filtering based on expression criteria to remove lowly expressed genes (TRUE/FALSE).
#'@param percent Percentage of samples showing minimum expression threshold for filtering (50\% by default).
#'@param cpm counts-per-million (CPM) expression threshold for filtering lowly expressed genes (1 CPM by default).


#'@import edgeR
#'@importFrom stats model.matrix
#'@importFrom edgeR calcNormFactors
#'@importFrom edgeR DGEList
#'@importFrom edgeR estimateDisp
#'@importFrom edgeR glmQLFit
#'@importFrom edgeR glmQLFTest
#'@importFrom edgeR topTags
#'@importFrom edgeR glmFit
#'@importFrom edgeR glmLRT
#'@importFrom methods setClass
#'@importFrom methods new

#'@return {object with four tables.
#'resDE = Differential Expression analysis for exonic counts using glmQLFTest from edgeR package.
#'resTc = EISA for the transcriptional component effect on each analyzed gene.
#'resPTc = EISA for the post-transcriptional component effect on each analyzed gene.
#'Expr_Int = Normalized log2 expression matrix for Intronic counts.
#'Expr_Ex = Normalized log2 expression matrix for Exonic counts}
#'
#'@export getEISAcompR
#'@references
#' \enumerate{
#'  \item Gaidatzis D et al. (2015) Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. Nature Biotechnology, 33, 722–729.
#'  \item Lawrence M et al. (2013) Software for computing and annotating genomic ranges. PLoS Computational Biology, 9, e1003118.
#'  \item Liao Y et al. (2019) The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleid Acids Research, 47, e47.
#'  \item Robinson MD et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139–140.
#'}
#'
#'@author Emilio Mármol Sánchez
#'@examples
#' {
#' \dontrun{
#'eisa <- getEISAcomp(exons=exon_counts, introns=intron_counts, design=design_matrix,
#' filterExpr=TRUE, percent=0.5, cpm=1)
#' }
#'
#' }
#'@name getEISAcompR
#'@rdname getEISAcompR-getEISAcompR

getEISAcompR <- function(exons, introns, design,
                         filterExpr=TRUE, percent=0.5, cpm=1){


  ## Define design
  design0 <- design[,1:2]
  design.m <- stats::model.matrix(~0+design0[,2])
  ncol_design <- ncol(design)
  i1 <- as.vector(which(design.m[,1]==1))
  i2 <- as.vector(which(design.m[,2]==1))
  exons.1 <- exons[,i1]
  exons.2 <- exons[,i2]
  introns.1 <- introns[,i1]
  introns.2 <- introns[,i2]

  exonsf <- exons
  intronsf <- introns

  ## Normalize
  exonsf0 <- mean(colSums(exonsf))*(exonsf)/colSums(exonsf)
  intronsf0 <- mean(colSums(intronsf))*(intronsf)/colSums(intronsf)


  if(filterExpr==TRUE){

    ## Filter by Expression
    message("Filtering by Expression thresholds ... ", appendLF=FALSE)

    n <- ncol(exonsf0)
    exons_dge <- edgeR::DGEList(exonsf0)
    exons_keep <- rowSums(cpm(exons_dge)>cpm) >= round(n*percent)
    exonsf0 <- exonsf0[exons_keep,]

    introns_dge <- edgeR::DGEList(intronsf0)
    introns_keep <- rowSums(cpm(introns_dge)>cpm) >= round(n*percent)
    intronsf0 <- intronsf0[introns_keep,]

    message("OK")

    ## Log2 Transform
    message("Computing DiffEx/DiffInt Estimates ... ", appendLF=FALSE)

    exonsf.1 <- log2((exonsf0[,i1])+1)
    exonsf.2 <- log2((exonsf0[,i2])+1)
    intronsf.1 <- log2((intronsf0[,i1])+1)
    intronsf.2 <- log2((intronsf0[,i2])+1)
    exonsflog <- cbind(exonsf.1, exonsf.2)
    intronsflog <- cbind(intronsf.1, intronsf.2)
    mergedflog <- merge(exonsflog, intronsflog, by="row.names")
    rownames(mergedflog) <- mergedflog$Row.names
    mergedflog <- mergedflog[,-1]
    exonsflog <- mergedflog[,1:n]
    colnames(exonsflog) <- colnames(exonsf0)
    intronsflog <- mergedflog[,(n+1):(n*2)]
    colnames(intronsflog) <- colnames(intronsf0)

    ## Diff Exons/Introns
    MeanEx.1 <- rowMeans(exonsf.1)
    MeanEx.2 <- rowMeans(exonsf.2)
    DiffEx <- as.data.frame(MeanEx.2 - MeanEx.1)
    rownames(DiffEx) <- rownames(exonsf0)
    colnames(DiffEx) <- c("DiffEx")

    MeanInt.1 <- rowMeans(intronsf.1)
    MeanInt.2 <- rowMeans(intronsf.2)
    DiffInt <- as.data.frame(MeanInt.2 - MeanInt.1)
    rownames(DiffInt) <- rownames(intronsf0)
    colnames(DiffInt) <- c("DiffInt")


    ## Merge DiffEx y DiffInt
    merged_Diff <- merge(DiffEx, DiffInt, by="row.names")

    ## Estimate Tc and PTc components
    merged_Diff$DiffExDiffInt <- merged_Diff$DiffEx - merged_Diff$DiffInt
    rownames(merged_Diff) <- merged_Diff$Row.names
    merged_Diff <- merged_Diff[,-1]

    message("OK")

    ## Filter by Expression and Correct
    iex <- which(rownames(exonsf) %in% rownames(merged_Diff))
    iint <- which(rownames(intronsf) %in% rownames(merged_Diff))

    exonsf <- exonsf[iex,]
    intronsf <- intronsf[iint,]

    merged_ei <- data.frame(exonsf, intronsf)

    factor.group <- rep(design.m[,2],2)
    factor.group2 <- design.m[,2]
    factor.type <- c(rep(1,n), rep(0,n))


    ## Perform E/I DE Analysis
    message("Performing E/I DE Analysis ... ", appendLF=FALSE)

    if(ncol(design) == 3){

      design.batch <- stats::model.matrix(~0+design[,2]+design[,3])
      design.batch2 <- rbind(design.batch[,], design.batch[,])
      design.batch2 <- design.batch2[,3]
      design.batch <- design.batch[,3]

      designf <- stats::model.matrix(~factor.group*factor.type+design.batch2)
      exons_dge <- edgeR::DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- edgeR::calcNormFactors(exons_dge, method="TMM")
      exons_dge <- edgeR::estimateDisp(exons_dge, designf, robust=T)
      exons_fit <-edgeR::glmQLFit(exons_dge, designf)
      exons_qlf <- edgeR::glmQLFTest(exons_fit)
      exons_results <- edgeR::topTags(exons_qlf, n=nrow(exons_dge))

      designf2 <- stats::model.matrix(~factor.group2+design.batch)
      introns_dge <- edgeR::DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- edgeR::calcNormFactors(introns_dge, method="TMM")
      introns_dge <- edgeR::estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <-edgeR::glmQLFit(introns_dge, designf2)
      introns_qlf <- edgeR::glmQLFTest(introns_fit)
      introns_results <- edgeR::topTags(introns_qlf, n=nrow(introns_dge))

      designf3 <- stats::model.matrix(~factor.group2+design.batch)
      DE_dge <- edgeR::DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- edgeR::calcNormFactors(DE_dge, method="TMM")
      DE_dge <- edgeR::estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <-edgeR::glmQLFit(DE_dge, designf3)
      DE_qlf <- edgeR::glmQLFTest(DE_fit)
      DE_results <- edgeR::topTags(DE_qlf, n=nrow(DE_dge))

      message("OK")

    } else if(ncol(design) == 2) {

      designf <- stats::model.matrix(~factor.group*factor.type)
      exons_dge <- edgeR::DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- edgeR::calcNormFactors(exons_dge, method="TMM")
      exons_dge <- edgeR::estimateDisp(exons_dge, designf, robust=T)
      exons_fit <- edgeR::glmFit(exons_dge, designf)
      exons_qlf <- edgeR::glmLRT(exons_fit)
      exons_results <- edgeR::topTags(exons_qlf, n=nrow(exons_dge))

      designf2 <- stats::model.matrix(~factor.group2)
      introns_dge <- edgeR::DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- edgeR::calcNormFactors(introns_dge, method="TMM")
      introns_dge <- edgeR::estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <- edgeR::glmFit(introns_dge, designf2)
      introns_qlf <- edgeR::glmLRT(introns_fit)
      introns_results <- edgeR::topTags(introns_qlf, n=nrow(introns_dge))

      designf3 <- stats::model.matrix(~factor.group2)
      DE_dge <- edgeR::DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- edgeR::calcNormFactors(DE_dge, method="TMM")
      DE_dge <- edgeR::estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <- edgeR::glmFit(DE_dge, designf3)
      DE_qlf <- edgeR::glmLRT(DE_fit)
      DE_results <- edgeR::topTags(DE_qlf, n=nrow(DE_dge))

      message("OK")

    }


    ## Merge Results
    merged_PTc <- merge(merged_Diff, exons_results$table, by="row.names")
    merged_PTc <- data.frame(cbind(merged_PTc$logFC, merged_PTc$DiffEx,
                                   scale(merged_PTc$DiffExDiffInt),
                                   merged_PTc$PValue, merged_PTc$FDR))
    rownames(merged_PTc) <- rownames(merged_Diff)
    colnames(merged_PTc) <- c("log2FC", "DiffEx", "PTc", "Pvalue", "FDR")

    merged_Tc <- merge(merged_Diff, introns_results$table, by="row.names")
    merged_Tc <- data.frame(cbind(merged_Tc$logFC, merged_Tc$DiffInt,
                                  scale(merged_Tc$DiffInt),
                                  merged_Tc$PValue, merged_Tc$FDR))
    rownames(merged_Tc) <- rownames(merged_Diff)
    colnames(merged_Tc) <- c("log2FC", "DiffInt", "Tc", "Pvalue", "FDR")

    log2CPM_exonsf.1 <- log2(cpm(exons_dge[,i1])+1)
    log2CPM_exonsf.2 <- log2(cpm(exons_dge[,i2])+1)
    Meanlog2CPM.1 <- data.frame(rowMeans(log2CPM_exonsf.1))
    Meanlog2CPM.2 <- data.frame(rowMeans(log2CPM_exonsf.2))

    Table_DE <- data.frame(cbind(DE_results$table$logFC, DE_results$table$PValue,
                                 DE_results$table$FDR))
    rownames(Table_DE) <- DE_results$table$genes

    merged_DE <- merge(Table_DE, Meanlog2CPM.1, by="row.names")
    merged_DE <- merge(merged_DE, Meanlog2CPM.2, by.x="Row.names",
                       by.y="row.names")
    rownames(merged_DE) <- merged_DE$Row.names
    merged_DE <- merged_DE[,-1]
    colnames(merged_DE) <- c("log2FC", "Pvalue", "FDR", "log2CPM.1", "log2CPM.2")


    ## Store results in object
    methods::setClass("EISAcompR",
             slots = list(resPTc = "data.frame", resTc = "data.frame", resDE = "data.frame",
                          Expr_Int = "data.frame", Expr_Ex = "data.frame"))
    results <- methods::new("EISAcompR", resPTc = data.frame(merged_PTc),
                   resTc = data.frame(merged_Tc), resDE = data.frame(merged_DE),
                   Expr_Int = data.frame(intronsflog), Expr_Ex = data.frame(exonsflog))


  } else if(filterExpr==FALSE){

    ## Log2 Transform
    message("Computing DiffEx/DiffInt Estimates ... ", appendLF=FALSE)

    n <- ncol(exonsf0)

    exonsf.1 <- log2((exonsf0[,i1])+1)
    exonsf.2 <- log2((exonsf0[,i2])+1)
    intronsf.1 <- log2((intronsf0[,i1])+1)
    intronsf.2 <- log2((intronsf0[,i2])+1)
    exonsflog <- cbind(exonsf.1, exonsf.2)
    intronsflog <- cbind(intronsf.1, intronsf.2)
    mergedflog <- merge(exonsflog, intronsflog, by="row.names")
    rownames(mergedflog) <- mergedflog$Row.names
    mergedflog <- mergedflog[,-1]
    exonsflog <- mergedflog[,1:n]
    colnames(exonsflog) <- colnames(exonsf0)
    intronsflog <- mergedflog[,(n+1):(n*2)]
    colnames(intronsflog) <- colnames(intronsf0)

    ## Diff Exons/Introns
    MeanEx.1 <- rowMeans(exonsf.1)
    MeanEx.2 <- rowMeans(exonsf.2)
    DiffEx <- as.data.frame(MeanEx.2 - MeanEx.1)
    rownames(DiffEx) <- rownames(exonsf0)
    colnames(DiffEx) <- c("DiffEx")

    MeanInt.1 <- rowMeans(intronsf.1)
    MeanInt.2 <- rowMeans(intronsf.2)
    DiffInt <- as.data.frame(MeanInt.2 - MeanInt.1)
    rownames(DiffInt) <- rownames(intronsf0)
    colnames(DiffInt) <- c("DiffInt")


    ## Merge DiffEx y DiffInt
    merged_Diff <- merge(DiffEx, DiffInt, by="row.names")

    ## Estimate Tc and PTc components
    merged_Diff$DiffExDiffInt <- merged_Diff$DiffEx - merged_Diff$DiffInt
    rownames(merged_Diff) <- merged_Diff$Row.names
    merged_Diff <- merged_Diff[,-1]

    message("OK")

    ## Filter by Expression and Correct
    iex <- which(rownames(exonsf) %in% rownames(merged_Diff))
    iint <- which(rownames(intronsf) %in% rownames(merged_Diff))

    exonsf <- exonsf[iex,]
    intronsf <- intronsf[iint,]

    merged_ei <- data.frame(exonsf, intronsf)

    factor.group <- rep(design.m[,2],2)
    factor.group2 <- design.m[,2]
    factor.type <- c(rep(1,n), rep(0,n))


    ## Perform E/I DE Analysis
    message("Performing E/I DE Analysis ... ", appendLF=FALSE)

    if(ncol(design) == 3) {

      design.batch <- stats::model.matrix(~0+design[,2]+design[,3])
      design.batch2 <- rbind(design.batch[,], design.batch[,])
      design.batch2 <- design.batch2[,3]
      design.batch <- design.batch[,3]

      designf <- stats::model.matrix(~factor.group*factor.type+design.batch2)
      exons_dge <- edgeR::DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- edgeR::calcNormFactors(exons_dge, method="TMM")
      exons_dge <- edgeR::estimateDisp(exons_dge, designf, robust=T)
      exons_fit <-edgeR::glmQLFit(exons_dge, designf)
      exons_qlf <- edgeR::glmQLFTest(exons_fit)
      exons_results <- edgeR::topTags(exons_qlf, n=nrow(exons_dge))

      designf2 <- stats::model.matrix(~factor.group2+design.batch)
      introns_dge <- edgeR::DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- edgeR::calcNormFactors(introns_dge, method="TMM")
      introns_dge <- edgeR::estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <-edgeR::glmQLFit(introns_dge, designf2)
      introns_qlf <- edgeR::glmQLFTest(introns_fit)
      introns_results <- edgeR::topTags(introns_qlf, n=nrow(introns_dge))

      designf3 <- stats::model.matrix(~factor.group2+design.batch)
      DE_dge <- edgeR::DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- edgeR::calcNormFactors(DE_dge, method="TMM")
      DE_dge <- edgeR::estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <-edgeR::glmQLFit(DE_dge, designf3)
      DE_qlf <- edgeR::glmQLFTest(DE_fit)
      DE_results <- edgeR::topTags(DE_qlf, n=nrow(DE_dge))

      message("OK")

    } else if(ncol(design) == 2) {

      designf <- stats::model.matrix(~factor.group*factor.type)
      exons_dge <- edgeR::DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- edgeR::calcNormFactors(exons_dge, method="TMM")
      exons_dge <- edgeR::estimateDisp(exons_dge, designf, robust=T)
      exons_fit <-edgeR::glmQLFit(exons_dge, designf)
      exons_qlf <- edgeR::glmQLFTest(exons_fit)
      exons_results <- edgeR::topTags(exons_qlf, n=nrow(exons_dge))

      designf2 <- stats::model.matrix(~factor.group2)
      introns_dge <- edgeR::DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- edgeR::calcNormFactors(introns_dge, method="TMM")
      introns_dge <- edgeR::estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <-edgeR::glmQLFit(introns_dge, designf2)
      introns_qlf <- edgeR::glmQLFTest(introns_fit)
      introns_results <- edgeR::topTags(introns_qlf, n=nrow(introns_dge))

      designf3 <- stats::model.matrix(~factor.group2)
      DE_dge <- edgeR::DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- edgeR::calcNormFactors(DE_dge, method="TMM")
      DE_dge <- edgeR::estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <-edgeR::glmQLFit(DE_dge, designf3)
      DE_qlf <- edgeR::glmQLFTest(DE_fit)
      DE_results <- edgeR::topTags(DE_qlf, n=nrow(DE_dge))

      message("OK")

    }


    ## Merge Results
    merged_PTc <- merge(merged_Diff, exons_results$table, by="row.names")
    merged_PTc <- data.frame(cbind(merged_PTc$logFC, merged_PTc$DiffEx,
                                   scale(merged_PTc$DiffExDiffInt),
                                   merged_PTc$PValue, merged_PTc$FDR))
    rownames(merged_PTc) <- rownames(merged_Diff)
    colnames(merged_PTc) <- c("log2FC", "DiffEx", "PTc", "Pvalue", "FDR")

    merged_Tc <- merge(merged_Diff, introns_results$table, by="row.names")
    merged_Tc <- data.frame(cbind(merged_Tc$logFC, merged_Tc$DiffInt,
                                  scale(merged_Tc$DiffInt),
                                  merged_Tc$PValue, merged_Tc$FDR))
    rownames(merged_Tc) <- rownames(merged_Diff)
    colnames(merged_Tc) <- c("log2FC", "DiffInt", "Tc", "Pvalue", "FDR")

    log2CPM_exonsf.1 <- log2(cpm(exons_dge[,i1])+1)
    log2CPM_exonsf.2 <- log2(cpm(exons_dge[,i2])+1)
    Meanlog2CPM.1 <- data.frame(rowMeans(log2CPM_exonsf.1))
    Meanlog2CPM.2 <- data.frame(rowMeans(log2CPM_exonsf.2))

    Table_DE <- data.frame(cbind(DE_results$table$logFC, DE_results$table$PValue,
                                 DE_results$table$FDR))
    rownames(Table_DE) <- DE_results$table$genes

    merged_DE <- merge(Table_DE, Meanlog2CPM.1, by="row.names")
    merged_DE <- merge(merged_DE, Meanlog2CPM.2, by.x="Row.names",
                       by.y="row.names")
    rownames(merged_DE) <- merged_DE$Row.names
    merged_DE <- merged_DE[,-1]
    colnames(merged_DE) <- c("log2FC", "Pvalue", "FDR", "log2CPM.1", "log2CPM.2")


    ## Store results in object
    methods::setClass("EISAcompR",
             slots = list(resPTc = "data.frame", resTc = "data.frame", resDE = "data.frame",
                          Expr_Int = "data.frame", Expr_Ex = "data.frame"))
    results <- methods::new("EISAcompR", resPTc = data.frame(merged_PTc),
                   resTc = data.frame(merged_Tc), resDE = data.frame(merged_DE),
                   Expr_Int = data.frame(intronsflog), Expr_Ex = data.frame(exonsflog))

  }

  message("")
  return(results)

}
