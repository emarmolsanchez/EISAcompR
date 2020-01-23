########################################
#
#  getEISAcomp function from EISACompR
#              EMS 2020
#
########################################

getEISAcomp <- function(exons, introns, design, capOut=TRUE,
                  filterExpr=TRUE, percent=0.5, cpm=1){
  suppressMessages(require(edgeR))
  
  ## Define design
  design.m <- model.matrix(~0+design[,2])
  i1 <- as.vector(which(design.m[,1]==1))
  i2 <- as.vector(which(design.m[,2]==1))
  exons.1 <- exons[,i1]
  exons.2 <- exons[,i2]
  introns.1 <- introns[,i1]
  introns.2 <- introns[,i2]
  
  if(capOut==TRUE){
    
    
    capOutlier <- function(x){
    qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
    caps <- quantile(x, probs=c(.1, .90), na.rm = T)
    H <- 1.5 * IQR(x, na.rm = T)
    x[x < (qnt[1] - H)] <- caps[1]
    x[x > (qnt[2] + H)] <- caps[2]
    return(x)
    
    }
  
  ## Cap Outliers
  message("")
  message("Capping Outliers ... ", appendLF=FALSE)
  
  exons.1_c <- round(t(apply(exons.1, 1, capOutlier)))
  exons.2_c <- round(t(apply(exons.2, 1, capOutlier)))
  introns.1_c <- round(t(apply(introns.1, 1, capOutlier)))
  introns.2_c <- round(t(apply(introns.2, 1, capOutlier)))
  
  message("OK")
  message("")
  
  ## merge
  exonsf <- cbind(exons.1_c, exons.2_c)
  intronsf <- cbind(introns.1_c, introns.2_c)
  
  } else if(capOut==FALSE){
    
    exonsf <- exons
    intronsf <- introns
  }
  
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
    merged_Diff$PTc <- scale(merged_Diff$DiffEx - merged_Diff$DiffInt)
    merged_Diff$Tc <- scale(merged_Diff$DiffInt)
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
    message("Performing DE Analysis ... ", appendLF=FALSE)
    
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
    
    message("OK")
    
    ## Merge Results
    
    merged_1 <- merge(merged_Diff, exons_results, by="row.names")
    merged_2 <- merge(merged_Diff, introns_results, by="row.names")
    Table1 <- cbind(merged_1$logFC, merged_1$DiffEx, merged_1$PTc,
                    merged_1$PValue, merged_1$FDR)
    rownames(Table1) <- merged_1$genes
    colnames(Table1) <- c("log2FC", "DiffEx", "PTc", "Pvalue","FDR")
    
    Table2 <- cbind(merged_2$logFC, merged_2$DiffInt, merged_2$Tc,
                    merged_2$PValue, merged_2$FDR)
    rownames(Table2) <- merged_2$genes
    colnames(Table2) <- c("log2FC", "DiffInt", "Tc", "Pvalue","FDR")
    
    ## Store results in object
    setClass("EISACompR",
             slots = list(resTc = "data.frame", resPTc = "data.frame",
                          Expr_Int = "data.frame", Expr_Ex = "data.frame"))
    results <- new("EISACompR", resTc = data.frame(Table2), resPTc = data.frame(Table1),
                   Expr_Int = data.frame(intronsflog), Expr_Ex = data.frame(exonsflog))
    
    
    } else if(filterExpr==FALSE){
    
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
      merged_Diff$PTc <- scale(merged_Diff$DiffEx - merged_Diff$DiffInt)
      merged_Diff$Tc <- scale(merged_Diff$DiffInt)
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
      message("Performing DE Analysis ... ", appendLF=FALSE)
      
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
      
      message("OK")
      
      ## Merge Results
      
      merged_1 <- merge(merged_Diff, exons_results, by="row.names")
      merged_2 <- merge(merged_Diff, introns_results, by="row.names")
      Table1 <- cbind(merged_1$logFC, merged_1$DiffEx, merged_1$PTc,
                        merged_1$PValue, merged_1$FDR)
      rownames(Table1) <- merged_1$genes
      colnames(Table1) <- c("log2FC", "DiffEx", "PTc", "Pvalue","FDR")
      
      Table2 <- cbind(merged_2$logFC, merged_2$DiffInt, merged_2$Tc,
                      merged_2$PValue, merged_2$FDR)
      rownames(Table2) <- merged_2$genes
      colnames(Table2) <- c("log2FC", "DiffInt", "Tc", "Pvalue","FDR")
      
      ## Store results in object
      setClass("EISACompR",
               slots = list(resTc = "data.frame", resPTc = "data.frame",
                            Expr_Int = "data.frame", Expr_Ex = "data.frame"))
      results <- new("EISACompR", resTc = data.frame(Table2), resPTc = data.frame(Table1),
                     Expr_Int = data.frame(intronsflog), Expr_Ex = data.frame(exonsflog))
      
}

  message("")
  return(results)
  
}