########################################
#
#  makeEISAgtfs function from EISACompR
#              EMS 2020
#
########################################


makeEISAgtfs <- function(annotFile, boundaryFix=10){
  
  message("")
  message("Loading GTF file ... ", appendLF=FALSE)
  suppressMessages(require(GenomicFeatures))
  suppressMessages(require(IRanges))
  setwd("~/")
  Dir0 <- "tmp_EISA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  workdir <- "~/tmp_EISA/"
  setwd(workdir)
  options(warn=-1)
  message("OK")
  
  GR <- makeTxDbFromGFF(annotFile)

  ## GR split GTF
  message("Splitting E/I features ... ", appendLF=FALSE)
  
  exons <- exonicParts(GR, linked.to.single.gene.only = F)
  introns <- intronicParts(GR, linked.to.single.gene.only = F)
  
  message("OK")

  ## Filter intronic ranges overlapping exons
  message("Removing overlapping exonic loci ... ", appendLF=FALSE)
  
  introns <- subsetByOverlaps(introns, exons, invert=T, ignore.strand=T)
  
  message("OK")
  
  ## Write exons/introns GTF
  message("Generating E/I gtfs ... ", appendLF=FALSE)
  
  ## Auxiliar functions
  makeGtfAttributes <- function(df, cols=NULL) {
    if (is.null(cols))
      cols = colnames(df)
    # make sure that gene_id and transcript_id are the first two columns
    mandatory = c("gene_id", "transcript_id")
    o = match(c(mandatory, setdiff(cols, mandatory)), cols)
    if (any(is.na(o[1:length(mandatory)]))) {
      o = o[!is.na(o)]
    }
    cols = cols[o]
    return(paste(apply(sapply(cols, function(s) {
      content = df[,s]
      if (is.character(content) | is.factor(content)) {
        content = paste('"', content, '"', sep="")
      }
      paste(gsub(".", "_", s, fixed=T), content, sep=" ")
    }), 1, paste, collapse="; "), ";", sep=""))
  }
  
  GR2gtf <- function(regions, filename, feature.type, src="EISACompR",
                     score=".", phase=".", attributes=NULL, ...) {
    
    strnd = as.character(strand(regions))
    strnd[strnd == "*"] = "."
    
    tab = data.frame(as.character(seqnames(regions)), src, feature.type,
                     as.numeric(start(regions)), as.numeric(end(regions)),
                     score, strnd, phase,
                     makeGtfAttributes(as.data.frame(elementMetadata(regions)),
                                       cols=attributes), stringsAsFactors=F)
    
    write.table(tab, file=filename, sep="\t", quote=F, row.names=F,
                col.names=F, ...)
  }
  
  
  GR2gtf(exons, "exons.gtf", feature.type="exon")
  GR2gtf(introns, "introns.gtf", feature.type="intron")
  
  ## Read exons/introns GTF
  exons <- read.table("exons.gtf", sep="\t")
  introns <- read.table("introns.gtf", sep="\t")
  
  message("OK")
  
  ## Remove single nt positions
  message("Removing singleton positions ... ", appendLF=FALSE)
  
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
  
  message("OK")
  
  ## Modify boundaries
  message("Adjusting boundaries ... ", appendLF=FALSE)
  
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
  setClass("EISACompR",
           slots = list(exonsGTF = "data.frame", intronsGTF = "data.frame"))
  results <- new("EISACompR", exonsGTF = exons, intronsGTF = introns)
  
  message("OK")
  
  setwd("~/")
  unlink("tmp_EISA", recursive=T)
  
  return(results)
  
  
}
