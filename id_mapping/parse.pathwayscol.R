#Goal: parse pathwayscol.xlsx
#1. read the file
#2. use ECK_1st_table to map EcoCyc IDs to pathways
#3. objs at the end: 
'
EcoCycID.Pwys: This is a table representing pathwayscol in another way
ECK.Pathway_table
ECK.Pathway: ids-Data (pathway ID)-Data Type 
'


#1
library(xlsx)
pathwayscol=read.xlsx("Data/pathwayscol.xlsx",colIndex = c(1:2,112:221),rowIndex = 32:456,sheetName = "pathwayscol",stringsAsFactors=F)
#check if the reading is correct: write.csv(pathwayscol,"pathwayscol.csv")

#2

##EcoCyc ID - Pathway ID - Pathway Name
EcoCycID.Pwys=data.frame()
for(i in 1:dim(pathwayscol)[1]){
  
  j=3 #Gene IDs start from the 3rd column
  while(is.na(pathwayscol[i,j])==F){ ##empty cells were read as NAs from read.xlsx()
    EcoCycID=pathwayscol[i,j]
    
    
    EcoCycID.Pwys=rbind(EcoCycID.Pwys,
                        c(EcoCycID,pathwayscol$UNIQUE.ID[i],pathwayscol$NAME[i]),
                        stringsAsFactors=F)
    
    j=j+1
    if(j>dim(pathwayscol)[2]) break
  }
  
}
colnames(EcoCycID.Pwys)=c("EcoCycID","PathwayID","PathwayName")


###genes to pathways is a many to many relationship:
###length(unique(EcoCycID.Pwys$`EcoCyc ID`))  ###1074 unique genes
###length(EcoCycID.Pwys$`EcoCyc ID`) ###2681
###length(unique(EcoCycID.Pwys$`Pathway ID`)) ###424 unique pathways
###length(EcoCycID.Pwys$`Pathway ID`) ###2681

##Combine EcoCycID.Pwys with Nichols id -- This is a comprehensive obj that has Nichols' ID, all accession No., Pathway ID and Pathway name
name_ACC1_ACC2=sapply(seq_along(genes.dat),FUN=function(i){
  c(names(genes.dat[i]),genes.dat[[i]]$'ACCESSION-1 - ',genes.dat[[i]]$'ACCESSION-2 - ')
})
name_ACC1_ACC2=Reduce(rbind,name_ACC1_ACC2)
name_ACC1_ACC2=apply(name_ACC1_ACC2,2,as.character) #I want to get rid of the factor attribute
colnames(name_ACC1_ACC2)=c("EcoCycID","ACCESSION-1 - ","ACCESSION-2 - ")


ECK_1st_table$ids=as.numeric(ECK_1st_table$ids)#I want to get rid of the factor attribute



ECK.Pathway_table=merge(ECK_1st_table,EcoCycID.Pwys,by="EcoCycID",all.x=T)
ECK.Pathway_table$ids=as.character(ECK.Pathway_table$ids)
ECK.Pathway_table=ECK.Pathway_table[order(as.numeric(ECK.Pathway_table$ids)),
                                    c("ids","EcoCycID","Original Name",
                                       "sorted_ECK_missing_gene_names_added",
                                       "associated_gene_names",
                                       "other.synonyms",
                                       "PathwayID",
                                       "PathwayName"
                                       )] #reorder the rows (by id) and columns


#Add a column for no. of genes in pathways (defined by pathwayscol.xlsx)
no_gene_inPathway=apply(pathwayscol,1,FUN=function(row){
  n_gene_inPathway=sum(!is.na(row))-2 # -2 is to not count the first 2 columns
})

pathwayID_noGeneInPathway=cbind(pathwayscol$UNIQUE.ID,no_gene_inPathway) %>% as.data.frame
pathwayID_noGeneInPathway$V1=as.character(pathwayID_noGeneInPathway$V1)
pathwayID_noGeneInPathway$no_gene_inPathway=as.character(pathwayID_noGeneInPathway$no_gene_inPathway)

colnames(pathwayID_noGeneInPathway)=c("PathwayID","NoGeneInPathway_byPathwaysCol")


ECK.Pathway_table=left_join(ECK.Pathway_table,pathwayID_noGeneInPathway,by="PathwayID")

##Create the SQL table: id (From ECK)- Pathway - Data Type (Pathway)
ECK.Pathway=cbind(ECK.Pathway_table[,c("ids","PathwayID")],"Pathway",stringsAsFactors=F)
colnames(ECK.Pathway)=c("ids","Data","Data Type")




#save(EcoCycID.Pwys,ECK.Pathway_table,ECK.Pathway,file="Data/sourced/parse.pathwayscol.RData")


