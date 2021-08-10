#'@title makeEISAgtfs
#'@description  This is the first function of EISACompR. makeEISAgtfs aims to
#'build custom annotation files in GTF format for late reads
#'quantification spanning either exonic or intronic genomic ranges.
#'Reference GTF annotation for your species of interest is needed in
#' order to generate the corresponding exonic/intronic regions for
#'  running EISA analysis.

#'@param annotFile PATH to reference GTF annotation file.
#'@param path_temp_files Home is the default "~/". This argument indicates the path were the temp_dir with exons, introns info should be created.
#'@param boundaryFix Exon junction boundary correction threshold (10 bp by default).
#'@param show_message (logical) if TRUE the proccess history is printed. Default TRUE.

#'@import GenomicFeatures
#'@import IRanges
#'@importFrom GenomicFeatures makeTxDbFromGFF
#'@importFrom GenomicFeatures exonicParts
#'@importFrom GenomicFeatures intronicParts
#'@importFrom IRanges subsetByOverlaps
#'@importFrom utils read.table
#'@importFrom methods setClass
#'@importFrom methods new
#'
#'@return results. It returns a EIsacompR object with the exons and introns information
#'@export makeEISAgtfs
#'@author Emilio Mármol Sánchez
#'@examples
#' {
#' \dontrun{
#'GTFS<-makeEISAgtfs(annotFile="path_to_ann_file")
#' }
#'
#' }
#'
#'@name makeEISAgtfs
#'@rdname makeEISAgtfs-makeEISAgtfs
makeEISAgtfs <- function(annotFile, path_temp_files="~/",boundaryFix=10, show_message=TRUE){
  if(show_message){
    message("")
    message("Loading GTF file ... ", appendLF=FALSE)
  }

  setwd(path_temp_files)
  Dir0 <- "tmp_EISA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  workdir <- paste0(path_temp_files,Dir0)
  setwd(workdir)
  options(warn=-1)
  if(show_message){
    message("Done")

  }

  GR <- GenomicFeatures::makeTxDbFromGFF(annotFile)

  ## GR split GTF
  if(show_message){
    message("Splitting E/I features ... ", appendLF=FALSE)
    }

  exons <- GenomicFeatures::exonicParts(GR, linked.to.single.gene.only = F)
  introns <- GenomicFeatures::intronicParts(GR, linked.to.single.gene.only = F)
  if(show_message){
    message("Done")
  }


  ## Filter intronic ranges overlapping exons
  if(show_message){
    message("Removing overlapping exonic loci ... ", appendLF=FALSE)
  }


  introns <- IRanges::subsetByOverlaps(introns, exons, invert=T, ignore.strand=T)
  if(show_message){
    message("Done")
    ## Write exons/introns GTF
    message("Generating E/I gtfs ... ", appendLF=FALSE)

  }



  GR2gtf(exons, "exons.gtf", feature.type="exon")
  GR2gtf(introns, "introns.gtf", feature.type="intron")

  ## Read exons/introns GTF
  #add fill
  exons <- utils::read.table("exons.gtf", sep="\t",fill=TRUE)
  introns <- utils::read.table("introns.gtf", sep="\t",fill=TRUE)

  if(show_message){
    message("Done")
    ## Remove single nt positions
    message("Removing singleton positions ... ", appendLF=FALSE)

  }

  i_ex <- which(exons$V4==exons$V5)

  if (length(i_ex)!= 0){
    exons <- exons[-i_ex,]
  } else if (length(i_ex)==0) {
    exons <- exons
  }

  i_int <- which(introns$V4==introns$V5)

  if (length(i_int)!= 0){
    introns <- introns[-i_int,]
  } else if (length(i_int)==0) {
    introns <- introns}


  if(show_message){
    message("Done")
    ## Modify boundaries
    message("Adjusting boundaries ... ", appendLF=FALSE)

  }

  exons$V4 <- exons$V4-boundaryFix
  exons$V5 <- exons$V5+boundaryFix
  exons$V4 <- ifelse(exons$V4<0, 1, exons$V4)

  i_ex <- which(exons$V4==exons$V5)

  if (length(i_ex)!=0){
  exons <- exons[-i_ex,]
  } else if (length(i_ex)==0){
    exons <- exons}

  i_ex <- which(exons$V4>exons$V5)

  if (length(i_ex)!=0){
  exons <- exons[-i_ex,]
  } else if (length(i_ex)==0){
    exons <- exons}

  introns$V4 <- introns$V4+boundaryFix
  introns$V5 <- introns$V5-boundaryFix
  introns$V4 <- ifelse(introns$V4<0, 1, introns$V4)

  i_int <- which(introns$V4==introns$V5)

  if (length(i_int)!=0){
    introns <- introns[-i_int,]
  } else if (length(i_int)==0){
    introns <- introns}

  i_int <- which(introns$V4>introns$V5)

  if (length(i_int)!=0){
    introns <- introns[-i_int,]
  } else if (length(i_int)==0){
    introns <- introns}

  exons <- exons[-grep("gene_id c", exons$V9), ]
  introns <- introns[-grep("gene_id c", introns$V9), ]


  ## Store results in object
  methods::setClass("EISAcompR",
           slots = list(exonsGTF = "data.frame", intronsGTF = "data.frame"))
  results <- methods::new("EISAcompR", exonsGTF = exons, intronsGTF = introns)

  if(show_message){
    message("Done")
    }


  setwd(path_temp_files)
  unlink("tmp_EISA", recursive=T)

  return(results)


}
