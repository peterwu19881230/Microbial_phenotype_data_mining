#Goal: Map synonyms for Nichols' gene names


##Read genes.col file
###genes.col<-read.table("ECK_name_mapping/genes.col.txt",skip=31,fill=TRUE,sep="\t",header=TRUE,stringsAsFactors = FALSE) => This will only get me 18XX rows. I don't know why
genes.col<-read.csv("ECK_name_mapping/genes.col.csv",skip=31,header=TRUE,stringsAsFactors = FALSE) 

## the sep argument says that "tab" was used to separte the original data

##Extract all the names needed
geneNames_genes.col<-genes.col[,c("NAME","SYNONYMS","SYNONYMS.1","SYNONYMS.2","SYNONYMS.3")]


##Concatenate them
all_gene_synonyms<-c()
for(i in 1:dim(geneNames_genes.col)[1]){
  all_gene_synonyms[i]<-geneNames_genes.col[i,1]
          for(j in 2:dim(geneNames_genes.col)[2]){
            if(geneNames_genes.col[i,j]!=""){
              all_gene_synonyms[i]<-paste(all_gene_synonyms[i],geneNames_genes.col[i,j],sep=",")
            }
            
          }
}
##There is a hole in geneNames_genes.col: There are 2 NULLs in row . Don't know why. I removed it.
##Result: A vector called: all_gene_synonyms




