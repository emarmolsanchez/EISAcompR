########################################
#
#  getEISAcounts function from EISACompR
#              EMS 2020
#
########################################


getEISAcounts <- function(files, annotFile, strandness, nthreads, PairedEnd=TRUE){
  suppressMessages(require(Rsubread))
  
  message("Loading Annotation File ... ", appendLF=FALSE)
  
  gtf <- read.table(annotFile, sep="\t")
  
  message("OK")
  
  message("Building Count Matrix ...", appendLF=FALSE)
  
  if(gtf$V3[1]=="exon"){
  
      counts <- featureCounts(files, annot.ext=annotFile,
                      strandSpecific = strandness,
                      isGTFAnnotationFile = T,
                      GTF.featureType = "exon",
                      nthreads=nthreads, isPairedEnd=PairedEnd,
                      countMultiMappingReads = F)
      
  
  } else if(gtf$V3[1]=="intron" && PairedEnd==TRUE){
    
      counts <- featureCounts(files, annot.ext=annotFile,
                            strandSpecific = strandness,
                            isGTFAnnotationFile = T,
                            GTF.featureType = "intron",
                            nthreads=nthreads, isPairedEnd=T,
                            countMultiMappingReads = F)
      
    
  } else if(gtf$V3[1]=="exon" && PairedEnd==FALSE){
    
    counts <- featureCounts(files, annot.ext=annotFile,
                            strandSpecific = strandness,
                            isGTFAnnotationFile = T,
                            GTF.featureType = "exon",
                            nthreads=nthreads, isPairedEnd=F,
                            countMultiMappingReads = F)
    
  } else if(gtf$V3[1]=="intron" && PairedEnd==FALSE){
    
    counts <- featureCounts(files, annot.ext=annotFile,
                            strandSpecific = strandness,
                            isGTFAnnotationFile = T,
                            GTF.featureType = "intron",
                            nthreads=nthreads, isPairedEnd=F,
                            countMultiMappingReads = F)
  }
  
  return(counts)
  
}
