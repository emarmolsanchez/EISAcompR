#'@title getEISAcounts
#'@description  This function allows users to quantify the level of
#'expression of their aligned reads by making use of custom-created exonic/intronic
#'GTF annotation files. The featureCounts built-in function from Rsubread package is
#'used for such purpose. Users may optionally decide to use any other preferred
#'quantification pipeline provided the use of exonic/intronic GTF annotation
#'as a reference.
#'@param files PATH to input BAM/SAM files containing read mapping results.
#'@param annotFile PATH to exonic/intronic GTF annotation file.
#'@param strandness Strandness. Character. You have to select one out "unstranded","stranded_forward","stranded_reversed".
#'@param nthreads Number of available threads for parallel computation.
#'@param PairedEnd Boolean specifying if input sequencing data is of type Paired-end (TRUE/FALSE).

#'@import Rsubread
#'@importFrom Rsubread featureCounts
#'@importFrom utils read.table

#'
#'@return results. It will create a raw count matrix stored at counts object with quantification estimates of each mapped read to the corresponding exonic/intronic regions and gene assignment.
#'@export getEISAcounts
#'@name getEISAcounts
#'@rdname getEISAcounts-getEISAcounts
#'@author Emilio Mármol Sánchez
#'@examples
#' {
#' \dontrun{
#'exon_counts <- getEISAcounts(files=vector_of_input_BAM/SAM, annotFile="PATH_to_exon_GTF",
#'strandness="unstranded", nthreads=4, PairedEnd=TRUE)
#'intron_counts <- getEISAcounts(files=vector_of_input_BAM/SAM, annotFile="PATH_to_intron_GTF",
#'strandness="unstranded", nthreads=4, PairedEnd=TRUE)
#' }
#'
#' }
#'@name getEISAcounts
#'@rdname getEISAcounts

getEISAcounts <- function(files, annotFile, strandness, nthreads, PairedEnd=TRUE){

if (!strandness%in%c("unstranded","stranded_forward","stranded_reversed")){
  stop("strandness has to be one of unstranded,stranded_forward,stranded_reversed")
}

if(strandness=="unstranded"){
 strandness=0
}
if(strandness=="stranded_forward"){
strandness=1
}
if(strandness=="stranded_reversed"){
    strandness=2
}

  message("Loading Annotation File ... ", appendLF=FALSE)

  gtf <- utils::read.table(annotFile, sep="\t")

  message("OK")

  message("Building Count Matrix ...", appendLF=FALSE)

  if(gtf$V3[1]=="exon"){

      counts <- Rsubread::featureCounts(files, annot.ext=annotFile,
                      strandSpecific = strandness,
                      isGTFAnnotationFile = T,
                      GTF.featureType = "exon",
                      nthreads=nthreads, isPairedEnd=PairedEnd,
                      countMultiMappingReads = F)


  } else if(gtf$V3[1]=="intron" && PairedEnd==TRUE){

      counts <- Rsubread::featureCounts(files, annot.ext=annotFile,
                            strandSpecific = strandness,
                            isGTFAnnotationFile = T,
                            GTF.featureType = "intron",
                            nthreads=nthreads, isPairedEnd=T,
                            countMultiMappingReads = F)


  } else if(gtf$V3[1]=="exon" && PairedEnd==FALSE){

    counts <- Rsubread::featureCounts(files, annot.ext=annotFile,
                            strandSpecific = strandness,
                            isGTFAnnotationFile = T,
                            GTF.featureType = "exon",
                            nthreads=nthreads, isPairedEnd=F,
                            countMultiMappingReads = F)

  } else if(gtf$V3[1]=="intron" && PairedEnd==FALSE){

    counts <- Rsubread::featureCounts(files, annot.ext=annotFile,
                            strandSpecific = strandness,
                            isGTFAnnotationFile = T,
                            GTF.featureType = "intron",
                            nthreads=nthreads, isPairedEnd=F,
                            countMultiMappingReads = F)
  }

  return(counts)

}

