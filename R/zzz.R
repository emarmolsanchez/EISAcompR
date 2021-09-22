
## Auxiliar functions
#'@importFrom BiocGenerics strand
#'@importFrom BiocGenerics start
#'@importFrom BiocGenerics end
#'@importForm S4Vectors elementMetadata
#'@importFrom GenomeInfoDb seqnames
#'@importFrom utils write.table
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


GR2gtf <- function(regions, feature.type, src="EISAcompR",
                   score=".", phase=".", attributes=NULL, ...) {

  strnd = as.character(BiocGenerics::strand(regions))
  strnd[strnd == "*"] = "."

  tab = data.frame(as.character(GenomeInfoDb::seqnames(regions)), src, feature.type,
                   as.numeric(BiocGenerics::start(regions)), as.numeric(BiocGenerics::end(regions)),
                   score, strnd, phase,
                   makeGtfAttributes(as.data.frame(S4Vectors::elementMetadata(regions)),
                                     cols=attributes), stringsAsFactors=F)

  colnames(tab) = c("chr", "src", "feature", "start", "end", "score",
                      "strnd", "phase", "attrb")
  
    return(tab)
}
