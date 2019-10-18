#Goal: 
#1.Make a .csv that contains the unique associated gene names in Nichols'

UniqueGeneNames<-unique(ECK_name_columns[,"associated_gene_names"]) #Total number: 3932

library(stringr)


formatted_UniqueGeneNames<-sapply(UniqueGeneNames,function(x)
  {
  sub("[A-Za-z]{3}",str_extract(x,"[A-Za-z]{3}") %>% tolower,x)
  }
)

#clean the exception of the above regular expression
formatted_UniqueGeneNames["TP2"]<-"tp2"

write.csv(formatted_UniqueGeneNames,file="formatted_unique_gene_names.csv")


# Make a .csv that contains the all the ECK designation (The missing genes for some ECK have been added,and the duplicates will be labeled with -1, -2) and their associated gene names
ECK_formattedGeneNames<-ECK_name_columns[,c("sorted_ECK_missing_gene_names_added","associated_gene_names")]

library(stringr)

ECK_formattedGeneNames[,"associated_gene_names"]<-sapply(ECK_formattedGeneNames[,"associated_gene_names"],function(x)
{
  sub("[A-Za-z]{3}",str_extract(x,"[A-Za-z]{3}") %>% tolower,x)
}
)

#clean the exception of the above regular expression
ECK_formattedGeneNames[,"associated_gene_names"][ECK_formattedGeneNames[,"associated_gene_names"]=="TP2"]<-"tp2"

write.csv(ECK_formattedGeneNames,file="ECK_formattedGeneNames.csv")


