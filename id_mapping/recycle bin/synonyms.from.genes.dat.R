accessions=c()
synonyms=c()
for(i in 1:length(genes.dat)){
  accessions=c(accessions,names(genes.dat)[[i]])
  if(identical(grepl("SYNONYMS - ",genes.dat[[1]]),integer(0))==FALSE){ ##if synonyms exist for the gene
    synonyms[i]=paste(genes.dat[[i]][["SYNONYMS - "]],collapse = ",")
  }else(
    synonyms[i]="No synonyms"
  )
}
genes.dat.synonyms=cbind(accessions,synonyms)

