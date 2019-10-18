#Goal: prepare data for creating relational database for Nichols'-GO-Uniprot 
#==> This is the 1st table

source("id_mapping/clean_names.R")
library(stringr)

#Add Unique IDs and all synonyms based on genes.dat (from EcoCyc)
##Extract info from  genes.dat object created by parse.genes.dat.R
##(!)Note that some of the ECKs are pseudogenes. They are not in genes.dat

uniqueID=list()
synonyms=list()
accessions=list()
accessions_bNumber=list()
for(i in 1:length(genes.dat)){
  
  ###retrieve Unique ID
  uniqueID[i]=names(genes.dat)[[i]]
  
  ###retrieve b number
  accessions_bNumber[i]=paste(genes.dat[[i]][["ACCESSION-1 - "]],collapse = ",")
  
  ###retrieve ECK 
  accessions[i]=paste(genes.dat[[i]][["ACCESSION-2 - "]],collapse = ",") 
  ###subsetting a list element that doesn't exist won't return "subscript out of bounds". Instead, it returns NULL
  
  ###retrieve synonyms if any
  synonyms[i]=paste(genes.dat[[i]][["SYNONYMS - "]],collapse = ",") 
  

  }


indices=list()
for(i in 1:3979){
  ECK=str_extract(id_ECKs_CorrectedECKs_AssociatedGeneNames$`originalECKs`[i],"^ECK[0-9]{4}") 
  
  index=which(accessions %in% ECK)
  
  if( index %>% length >= 1 &is.numeric(index)==T){
    indices[[i]]=index
    }else{
    indices[[i]]=""
  } 
  
}
rm(ECK,index) 

###This loop checks if there are more than 2 ECKs found
two.or.more=list()
j=1
for(i in 1:length(indices)){
  if(length(indices[[i]])>=2){
    two.or.more[j]=indices[i]
    names(two.or.more)[j]=paste("index from indices= ",i)
    j=j+1
  }
}
two.or.more ###indices[687]=1859 2015  (ECK1366 was found twice in genes.dat) 

###This block prepares data to cbind(). Mainly because ECK1366 found 2 data in genes.dat
EcoCycID=c()
other.synonyms=c()
bNumber=c()
for(i in 1:3979){
  EcoCycID[i]=paste(uniqueID[indices[[i]]],collapse=",") ###Allowing concatnating a vector into 1 string
  other.synonyms[i]=paste(synonyms[indices[[i]]],collapse=",") 
  if(grepl("^,",other.synonyms[i])==T | other.synonyms[i]=="NULL" | other.synonyms[i]==""){
    other.synonyms[i]=NA ###If there are no synomym, the result will be NA
  }  ###remove "," without any synonym attached
  bNumber[i]=paste(accessions_bNumber[indices[[i]]],collapse=",")
}


###drop all levels in the column: sorted_ECK_missing_gene_names_added
id_ECKs_CorrectedECKs_AssociatedGeneNames$sorted_ECK_missing_gene_names_added=as.character(id_ECKs_CorrectedECKs_AssociatedGeneNames$sorted_ECK_missing_gene_names_added)


##The result obj. Note that ECK1366-LOMR' has 2 EcoCyc identifiers => we should use G6692 because LOMR' means the disruption on the N terminal (if it's 'LOMR it's C terminal)
ECK_1st_table=cbind(id_ECKs_CorrectedECKs_AssociatedGeneNames,
      EcoCycID=EcoCycID,
      other.synonyms=other.synonyms,bNumber,stringsAsFactors=F)

##Correct the row: ECK1366-LOMR'
ECK_1st_table[ECK_1st_table$ids=="687",]$EcoCycID="G6692"
ECK_1st_table[ECK_1st_table$ids=="687",]$bNumber="b1369"

##Add Ecocyc ID where strains don't have ECK number: tp2: G0-8894, tpkE70: G0-8906 , istR-1: G0-10202
ECK_1st_table[ECK_1st_table$associated_gene_names=="tp2","EcoCycID"]="G0-8894"
ECK_1st_table[ECK_1st_table$associated_gene_names=="tpkE70","EcoCycID"]="G0-8906"
ECK_1st_table[ECK_1st_table$associated_gene_names=="istR-1","EcoCycID"]="G0-10202"


ECK_1st_table=as.data.frame(ECK_1st_table)
names(ECK_1st_table)[2]="Original Name"




#The following code: 
##1. add a column for JW numbers by an original file provided by Nichols' paper
##2. add columns for the link to CGSC synonym, strain availability, position of deletion
##(problems are explained in each block)


library(stringr)
##1
strainInfo=read_xls("Data/inline-supplementary-material-4(Baba et al., 2006).xls",sheet=1,skip=2,col_names = T)[,1:3]
colnames(strainInfo)=c("ECK","gene","JW")
ECK_1st_table$ECK=str_extract(ECK_1st_table$sorted_ECK_missing_gene_names_added,"^ECK[0-9]{4}")

temp=merge(ECK_1st_table,strainInfo,by="ECK",all.x=T,all.y=F)
temp$ids=as.character(temp$ids)
temp=temp[order(as.numeric(temp$ids)),]

### A problem:
sum(duplicated(temp$ids)) #7 ECKs have more than 1 JW names (according to the file: inline-supplementary-material-4(Baba et al., 2006).xls). 

##This table gives the info of all those ECKs that have 2 JWs => The genes deleted are the same but the positions are different
duplicated.ECK=temp[temp$ECK %in% c("ECK0614","ECK1007","ECK1157","ECK1409","ECK1718","ECK1975","ECK2025"),]
##For example ECK0614: 
## http://cgsc2.biology.yale.edu/Strain.php?ID=115890 
## http://cgsc2.biology.yale.edu/Strain.php?ID=115891

##2
##Add 3 columns (CGSC synonym, strain availability, position of deletion) from: JWstrainInfoFromCGSC.txt created by parse_CGSC_Strains.R
JWstrainInfoFromCGSC=read.table("Data/JWstrainInfoFromCGSC.txt")
colnames(JWstrainInfoFromCGSC)=c("JW","JW2","position","strain_availibility")
temp2=merge(temp,JWstrainInfoFromCGSC,by="JW",all.x=T,all.y=F)
  

temp2=temp2[,c("ids","ECK","JW","JW2","associated_gene_names","position","Original Name","sorted_ECK_missing_gene_names_added","EcoCycID","bNumber","other.synonyms","strain_availibility")]

### A problem:
##JW1878 has 2 forms from CGSC: JW1878-4, JW1878-4/pWT2594 => I should remove the one with plasmid. It's unlikely Nichols used the strain with a plasmid encoding gfp and other things
JWstrainInfoFromCGSC$JW[which(duplicated(JWstrainInfoFromCGSC$JW))]

temp2=temp2[-which(temp2$JW2=="JW1878-4/pWTZ594"),]
temp2$bNumber[temp2$bNumber=="NULL"]=NA #clean the bNumber column

ECK_1st_table=temp2

##Check to see if the table is correct using Supplemental Information: Bacterial Strains from Nichols' paper
## 3737 Keio, 117 SPA, 5 DAS, 9 alleles, 9 linked controls, 2 truncations, and 100 sRNAs
sum(grepl(" - SPA",ECK_1st_table$sorted_ECK_missing_gene_names_added)) 
sum(grepl(" - DAS",ECK_1st_table$sorted_ECK_missing_gene_names_added)) 
sum(grepl(" - Linked",ECK_1st_table$sorted_ECK_missing_gene_names_added))
sum(grepl(" - Truncation",ECK_1st_table$sorted_ECK_missing_gene_names_added))

## I have to correct things that are not in the 3737 Keio (They do have ECK at the front but they weren't from Keio. However, the code above considered them to be):
hobbs=read_xls("Data/hobbs2009_tables1_smallrnasandproteins.xls",sheet=1,skip=9,col_names=F)[,1:2]
colnames(hobbs)=c("EcoCycID","genotype")
hobbsInNichols=merge(ECK_1st_table,hobbs,by="EcoCycID")

## add strain background column


save(ECK_1st_table,file="Data/ECK_1st_table.RData")


















##(!)Some of the ECKs don't match an EcoCyc gene ID. Here they are.

###The table
NULLgenes=ECK_1st_table[ECK_1st_table$EcoCycID=="NULL",]

###A list of ECKs for those genes
Nullgenes.list=sapply(NULLgenes$`Original Name`,str_extract,pattern="^ECK[0-9][0-9][0-9][0-9]") %>% as.character

###Output the list to a text file
###write.table(Nullgenes.list,"Nullgenes.list.txt",quote=F,row.names = F,col.names = F)

###Create a smart table and map them to EcoCyc gene ID (Object ID) and pathways. Luckily, none of them are in a pathway
eightyNulls=read.table("Data/80NULLs.txt")

##Also, I Got all the ECKs and made a smart table on EcoCyc 


###A list of all ECKs: ECKs_rows_fixed (defined in clean_names.R)


ECKs.list=sapply(ECKs_rows_fixed,str_extract,pattern="^ECK[0-9][0-9][0-9][0-9]") 
###write.table(ECKs.list,"Data/ECKs.list.txt",quote=F,row.names = F,col.names = F)





