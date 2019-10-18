#Goal: parse the excel file from copying-pasting from CGSC: http://cgsc2.biology.yale.edu/StrainQuery.php
##(I don't know how to use the url to parse it directly)

path <- rstudioapi::getActiveDocumentContext()$path #This sets the directory to where this script locates
Encoding(path) <- "UTF-8"
setwd(dirname(path)) 


library(xlsx)
library(stringr)
strainInfo=read.xlsx2("keioStrainsFromCGSC.xlsx",sheetIndex=1)
strainInfo$allele=sapply(strainInfo$Genotype,FUN=function(genotype){ #This guy gets like: "xerC757(del)::kan". 757 is the allele No.
  geneWithKan=str_extract(genotype,"Δ[a-zA-Z]{2,5}[A-Z]{0,1}-{0,1}[0-9]{1,4}::kan")  
  #This should only search and replace the 1st pattern found (because I put ::kan at the end)
  geneWithKan=sub("Δ","",geneWithKan) #This removes "Δ"
  geneWithKan=sub("::kan","(del)::kan",geneWithKan) #This adds (del) infront of "::kan"
  })

strainInfo$Name=str_extract(strainInfo$Name,"JW[0-9]{4}")

strainInfo=strainInfo[,c("Name","allele")]

#write.table(strainInfo,"JWstrainInfoFromCGSC.txt",row.names = F,col.names = F,quote=F,sep="\t")
