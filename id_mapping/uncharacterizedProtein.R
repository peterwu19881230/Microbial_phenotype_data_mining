#Uncharacterized proteins in Nichols'


#Get EcoCyc identifiers and protein info of uncharacterized proteins
gene_association.ecocyc=read.table("Data/EcoCyc_21.1//data/gene_association.ecocyc",sep="\t",quote = "",stringsAsFactors = F) 
## =>quote = "" is to disable quoting altogether => Fix the ' symbol (or other quoting symbols if any) being misread
accession_proteinDescrip=unique(gene_association.ecocyc[,c("V10","V11")])
dim(accession_proteinDescrip)
accession_uncharacterizedProtein=accession_proteinDescrip[grepl("uncharacterized protein",accession_proteinDescrip$V10),]
library(stringr)
accession_uncharacterizedProtein$V11=str_extract(accession_uncharacterizedProtein$V11,"ECK[0-9]{4}")

dim(accession_uncharacterizedProtein) #Only 22 uncharacterized proteins

#Get Nichols' id & EcoCycID 
id_ECK=unique(ECK_1st_table[,c("ids","EcoCycID")])


#Do CorrVSAnnot after removing uncharacterized proteins?


