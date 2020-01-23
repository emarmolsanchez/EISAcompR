########################################
#
#  writeEISAgtfs function from EISACompR
#              EMS 2020
#
########################################


writeEISAgtfs <- function(eisaR.obj, out.path){

  exonGTF <- eisaR.obj@exonsGTF
  intronGTF <- eisaR.obj@intronsGTF

  ## Set output path
  setwd(out.path)
  
  ## Write GTFs
  write.table(exonGTF, "exons.gtf", quote=F, 
              sep="\t", col.names=F, row.names=F)
  write.table(intronGTF, "introns.gtf", quote=F, 
              sep="\t", col.names=F, row.names=F)
  
  
}
