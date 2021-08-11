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
#'@importFrom data.table fread
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

  Exons_f <- GenomicFeatures::exonicParts(GR, linked.to.single.gene.only = F)
  Introns_f <- GenomicFeatures::intronicParts(GR, linked.to.single.gene.only = F)
  if(show_message){
    message("Done")
  }


  ## Filter intronic ranges overlapping Exons_f
  if(show_message){
    message("Removing overlapping exonic loci ... ", appendLF=FALSE)
  }


  Introns_f <- IRanges::subsetByOverlaps(Introns_f, Exons_f, invert=T, ignore.strand=T)
  if(show_message){
    message("Done")
    ## Write Exons/Introns GTF
    message("Generating E/I gtfs ... ", appendLF=FALSE)

  }



  GR2gtf(Exons_f, "Exons.gtf", feature.type="exon")
  GR2gtf(Introns_f, "Introns.gtf", feature.type="intron")

  ## Read Exons/Introns GTF
  #add fill
  #considering use fread instead
  Exons_fR <- data.table::fread("Exons.gtf", sep="\t",fill=TRUE)
  Introns_fR <- data.table::fread("Introns.gtf", sep="\t",fill=TRUE)

  if(show_message){
    message("Done")
    ## Remove single nt positions
    message("Removing singleton positions ... ", appendLF=FALSE)

  }

  i_ex <- which(Exons_fR$V4==Exons_fR$V5)

  if (length(i_ex)!= 0){
    Exons_fR <- Exons_fR[-i_ex,]
  } else if (length(i_ex)==0) {
    Exons_fR <- Exons_fR
  }

  i_int <- which(Introns_fR$V4==Introns_fR$V5)

  if (length(i_int)!= 0){
    Introns_fR <- Introns_fR[-i_int,]
  } else if (length(i_int)==0) {
    Introns_fR <- Introns_fR}


  if(show_message){
    message("Done")
    ## Modify boundaries
    message("Adjusting boundaries ... ", appendLF=FALSE)

  }

  Exons_fR$V4 <- Exons_fR$V4-boundaryFix
  Exons_fR$V5 <- Exons_fR$V5+boundaryFix
  Exons_fR$V4 <- ifelse(Exons_fR$V4<0, 1, Exons_fR$V4)

  i_ex <- which(Exons_fR$V4==Exons_fR$V5)

  if (length(i_ex)!=0){
  Exons_fR <- Exons_fR[-i_ex,]
  } else if (length(i_ex)==0){
    Exons_fR <- Exons_fR}

  i_ex <- which(Exons_fR$V4>Exons_fR$V5)

  if (length(i_ex)!=0){
  Exons_fR <- Exons_fR[-i_ex,]
  } else if (length(i_ex)==0){
    Exons_fR <- Exons_fR}

  Introns_fR$V4 <- Introns_fR$V4+boundaryFix
  Introns_fR$V5 <- Introns_fR$V5-boundaryFix
  Introns_fR$V4 <- ifelse(Introns_fR$V4<0, 1, Introns_fR$V4)

  i_int <- which(Introns_fR$V4==Introns_fR$V5)

  if (length(i_int)!=0){
    Introns_fR <- Introns_fR[-i_int,]
  } else if (length(i_int)==0){
    Introns_fR <- Introns_fR}

  i_int <- which(Introns_fR$V4>Introns_fR$V5)

  if (length(i_int)!=0){
    Introns_fR <- Introns_fR[-i_int,]
  } else if (length(i_int)==0){
    Introns_fR <- Introns_fR}

  Exons_fR <- Exons_fR[-grep("gene_id c", Exons_fR$V9), ]
  Introns_fR <- Introns_fR[-grep("gene_id c", Introns_fR$V9), ]


  ## Store results in object
  methods::setClass("EISAcompR",
           slots = list(exonsGTF = "data.frame", intronsGTF = "data.frame"))
  results <- methods::new("EISAcompR", exonsGTF = Exons_fR, intronsGTF = Introns_fR)

  if(show_message){
    message("Done")
    }


  setwd(path_temp_files)
  unlink("tmp_EISA", recursive=T)

  return(results)


}
