#Goal: create a subset of Nichols's strains such that all of them are deletions -> make the corresponding subset of strainCC file


#1st part
library(stringr)

#ECK
all_deletions_ECK_only<-str_extract(most_common_genes,"^ECK[0-9][0-9][0-9][0-9]")

#gene name
ECKXXXXremoved<-str_extract(most_common_genes,"-(.*)") ###This is based on the assumption that all the full names start with ECKXXXX-
all_deletions_gene_names<-str_extract(ECKXXXXremoved,"[A-Z]{3,4}") ###This is a little dangerous since whether the names are actual gene names cannot be verified. Will verify after mapping to gene names in the .GAF file
rm(ECKXXXXremoved)


all_deletions_nameTable<-cbind(most_common_genes, all_deletions_ECK_only, all_deletions_gene_names)


#2nd part
all_deletions_cor_strains<-cor_strains[most_common_genes,most_common_genes]


#ECK
ECK_all_deletions_cor_strains<-all_deletions_cor_strains
rownames(ECK_all_deletions_cor_strains)<-all_deletions_ECK_only
colnames(ECK_all_deletions_cor_strains)<-rownames(ECK_all_deletions_cor_strains)

#gene name
gene_names_all_deletions_cor_strains<-all_deletions_cor_strains #This is what Curtis wants
rownames(gene_names_all_deletions_cor_strains)<-all_deletions_gene_names
colnames(gene_names_all_deletions_cor_strains)<-rownames(gene_names_all_deletions_cor_strains)





