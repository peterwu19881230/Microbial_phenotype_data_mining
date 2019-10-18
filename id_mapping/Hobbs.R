#Hobbs table

library(xlsx)
library(stringr)
library(dplyr)

GSOdata=read.xlsx2("/Users/peterwu/TAMU Google Drive/OMP/Nichols - phenotype landscape paper/strains/hobbs2009_tables1_smallrnasandproteins.xls",
                   sheetIndex=1,startRow = 10,header = F,colIndex =1:2)

names(GSOdata)=c("GSO","strain_allele")


geneNames=GSOdata$strain_allele %>% sub(pattern="^MG1655 D",replacement="") %>% sub(pattern="::kan$",replacement="") 
# => From the above line the names aren't all clean

length(geneNames[geneNames %in% ECK_1st_table$associated_gene_names] %>% unique) #107


