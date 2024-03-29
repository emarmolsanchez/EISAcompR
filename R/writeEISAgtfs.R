#'@title writeEISAgtfs
#'@description  This function allows to export the exon/intron annnotation files in GTF format.

#'@param eisaR.obj EISACompR object previously stored with exon/intron split GTFs.
#'@param out.path PATH to desired output. Default is "~/".
#'@export writeEISAgtfs
#'@author Emilio Mármol Sánchez
#'@examples
#' {
#' \dontrun{
#'writeEISAgtfs(GTFs, "~/")
#' }
#'
#' }
#'@name writeEISAgtfs
#'@rdname writeEISAgtfs-writeEISAgtfs
#'@importFrom utils write.table

removeNewlines <- function(subGTF){
  corrctedGTF <- subGTF %>%
                 dplyr::mutate(attrb = gsub("\n", "", attrb),
                            start = as.integer(start),
                            end = as.integer(end)
                          )
   return(corrctedGTF)
  }

writeEISAgtfs <- function(eisaR.obj, out.path="~/"){

  exonGTF <- removeNewlines(eisaR.obj$exonsGTF)
             
  intronGTF <- removeNewlines(eisaR.obj$intronsGTF)

  ## Write GTFs
  utils::write.table(exonGTF, paste0(out.path,"Exons.gtf"), quote=F,
              sep="\t", col.names=F, row.names=F)
  utils::write.table(intronGTF, paste0(out.path,"Introns.gtf"), quote=F,
              sep="\t", col.names=F, row.names=F)


}
