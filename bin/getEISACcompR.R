########################################
#
#  getEISAcomp function from EISAcompR
#              EMS 2020
#
########################################

getEISAcompR <- function(exons, introns, design,
                        filterExpr=TRUE, percent=0.5, cpm=1){
  suppressMessages(require(edgeR))
  
  ## Define design
  design0 <- design[,1:2]
  design.m <- model.matrix(~0+design0[,2])
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
    exons_dge <- DGEList(exonsf0)
    exons_keep <- rowSums(cpm(exons_dge)>cpm) >= round(n*percent)
    exonsf0 <- exonsf0[exons_keep,]
    
    introns_dge <- DGEList(intronsf0)
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
      
      design.batch <- model.matrix(~0+design[,2]+design[,3])
      design.batch2 <- rbind(design.batch[,], design.batch[,])
      design.batch2 <- design.batch2[,3]
      design.batch <- design.batch[,3]
      
      designf <- model.matrix(~factor.group*factor.type+design.batch2)
      exons_dge <- DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- calcNormFactors(exons_dge, method="TMM")
      exons_dge <- estimateDisp(exons_dge, designf, robust=T)
      exons_fit <- glmQLFit(exons_dge, designf)
      exons_qlf <- glmQLFTest(exons_fit)
      exons_results <- topTags(exons_qlf, n=nrow(exons_dge))
      
      designf2 <- model.matrix(~factor.group2+design.batch)
      introns_dge <- DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- calcNormFactors(introns_dge, method="TMM")
      introns_dge <- estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <- glmQLFit(introns_dge, designf2)
      introns_qlf <- glmQLFTest(introns_fit)
      introns_results <- topTags(introns_qlf, n=nrow(introns_dge))
      
      designf3 <- model.matrix(~factor.group2+design.batch)
      DE_dge <- DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- calcNormFactors(DE_dge, method="TMM")
      DE_dge <- estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <- glmQLFit(DE_dge, designf3)
      DE_qlf <- glmQLFTest(DE_fit)
      DE_results <- topTags(DE_qlf, n=nrow(DE_dge))
      
      message("OK")
      
    } else if(ncol(design) == 2) {
      
      designf <- model.matrix(~factor.group*factor.type)
      exons_dge <- DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- calcNormFactors(exons_dge, method="TMM")
      exons_dge <- estimateDisp(exons_dge, designf, robust=T)
      exons_fit <- glmFit(exons_dge, designf)
      exons_qlf <- glmLRT(exons_fit)
      exons_results <- topTags(exons_qlf, n=nrow(exons_dge))
      
      designf2 <- model.matrix(~factor.group2)
      introns_dge <- DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- calcNormFactors(introns_dge, method="TMM")
      introns_dge <- estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <- glmFit(introns_dge, designf2)
      introns_qlf <- glmLRT(introns_fit)
      introns_results <- topTags(introns_qlf, n=nrow(introns_dge))
      
      designf3 <- model.matrix(~factor.group2)
      DE_dge <- DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- calcNormFactors(DE_dge, method="TMM")
      DE_dge <- estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <- glmFit(DE_dge, designf3)
      DE_qlf <- glmLRT(DE_fit)
      DE_results <- topTags(DE_qlf, n=nrow(DE_dge))
      
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
    merged_Tc <- data.frame(cbind(merged_Tc$logFC, merged_Tc$DiffEx, 
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
    setClass("EISACompR",
             slots = list(resPTc = "data.frame", resTc = "data.frame", resDE = "data.frame",
                          Expr_Int = "data.frame", Expr_Ex = "data.frame"))
    results <- new("EISACompR", resPTc = data.frame(merged_PTc), 
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
      
      design.batch <- model.matrix(~0+design[,2]+design[,3])
      design.batch2 <- rbind(design.batch[,], design.batch[,])
      design.batch2 <- design.batch2[,3]
      design.batch <- design.batch[,3]
      
      designf <- model.matrix(~factor.group*factor.type+design.batch2)
      exons_dge <- DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- calcNormFactors(exons_dge, method="TMM")
      exons_dge <- estimateDisp(exons_dge, designf, robust=T)
      exons_fit <- glmQLFit(exons_dge, designf)
      exons_qlf <- glmQLFTest(exons_fit)
      exons_results <- topTags(exons_qlf, n=nrow(exons_dge))
      
      designf2 <- model.matrix(~factor.group2+design.batch)
      introns_dge <- DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- calcNormFactors(introns_dge, method="TMM")
      introns_dge <- estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <- glmQLFit(introns_dge, designf2)
      introns_qlf <- glmQLFTest(introns_fit)
      introns_results <- topTags(introns_qlf, n=nrow(introns_dge))
      
      designf3 <- model.matrix(~factor.group2+design.batch)
      DE_dge <- DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- calcNormFactors(DE_dge, method="TMM")
      DE_dge <- estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <- glmQLFit(DE_dge, designf3)
      DE_qlf <- glmQLFTest(DE_fit)
      DE_results <- topTags(DE_qlf, n=nrow(DE_dge))
      
      message("OK")
      
    } else if(ncol(design) == 2) {
      
      designf <- model.matrix(~factor.group*factor.type)
      exons_dge <- DGEList(counts=merged_ei, genes=rownames(merged_ei))
      exons_dge <- calcNormFactors(exons_dge, method="TMM")
      exons_dge <- estimateDisp(exons_dge, designf, robust=T)
      exons_fit <- glmQLFit(exons_dge, designf)
      exons_qlf <- glmQLFTest(exons_fit)
      exons_results <- topTags(exons_qlf, n=nrow(exons_dge))
      
      designf2 <- model.matrix(~factor.group2)
      introns_dge <- DGEList(counts=intronsf, genes=rownames(intronsf))
      introns_dge <- calcNormFactors(introns_dge, method="TMM")
      introns_dge <- estimateDisp(introns_dge, designf2, robust=T)
      introns_fit <- glmQLFit(introns_dge, designf2)
      introns_qlf <- glmQLFTest(introns_fit)
      introns_results <- topTags(introns_qlf, n=nrow(introns_dge))
      
      designf3 <- model.matrix(~factor.group2)
      DE_dge <- DGEList(counts=exonsf, genes=rownames(exonsf))
      DE_dge <- calcNormFactors(DE_dge, method="TMM")
      DE_dge <- estimateDisp(DE_dge, designf3, robust=T)
      DE_fit <- glmQLFit(DE_dge, designf3)
      DE_qlf <- glmQLFTest(DE_fit)
      DE_results <- topTags(DE_qlf, n=nrow(DE_dge))
      
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
    merged_Tc <- data.frame(cbind(merged_Tc$logFC, merged_Tc$DiffEx, 
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
    setClass("EISACompR",
             slots = list(resPTc = "data.frame", resTc = "data.frame", resDE = "data.frame",
                          Expr_Int = "data.frame", Expr_Ex = "data.frame"))
    results <- new("EISACompR", resPTc = data.frame(merged_PTc), 
                   resTc = data.frame(merged_Tc), resDE = data.frame(merged_DE),
                   Expr_Int = data.frame(intronsflog), Expr_Ex = data.frame(exonsflog))
    
  }
  
  message("")
  return(results)
  
}
