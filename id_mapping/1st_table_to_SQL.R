##Create a table (ids - associated_gene_names - ECK - all_gene_synonyms) on MySQL

library(dplyr)
library(RMySQL)

load("ECK_1st_table.RData") #the dependency file


gene_name=ECK_1st_table$associated_gene_names 
other_synonym=strsplit(ECK_1st_table$other.synonyms,split=",",fixed=T) 

all_gene_synonyms=c()
for(i in 1:length(gene_name)){
  if(identical(character(0),other_synonym[[i]])){
    all_gene_synonyms[i]=gene_name[i]
  }else{
    all_gene_synonyms[i]=union(gene_name[i],other_synonym[[i]]) %>% paste(collapse=",") 
  }
  
}

table_=cbind(ECK_1st_table[,c("ids","associated_gene_names","ECK")],
             all_gene_synonyms) %>% unique #unique is to get rid of the propagated rows resulted from other colnames in ECK_1st_table



#Connect and create the table (if exists, rewrite)
mydb = dbConnect(MySQL(), user='peterwu1230', password='*********', dbname='chemgen', host='127.0.0.1',port=3308)
dbWriteTable(mydb, name='Genotypes', value=table_,row.names=FALSE,overwrite=T) ###I already have ids, so there's no need for row names







