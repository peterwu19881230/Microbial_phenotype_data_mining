#Goal: rename the colnames of the cor_strains obj I produced from Nichols' original data

#ECK No. only
cp_cor_strains<-cor_strains
rownames(cp_cor_strains)<-str_extract(rownames(cp_cor_strains),"^ECK[0-9][0-9][0-9][0-9]")
colnames(cp_cor_strains)<-rownames(cp_cor_strains)
strainCC_ECK_No_only<-cp_cor_strains
rm(cp_cor_strains)

#Gene names only
cp_cor_strains<-cor_strains

ECKXXXXremoved<-str_extract(rownames(cp_cor_strains),"-(.*)") ###This is based on the assumption that all the full names start with ECKXXXX-

rownames(cp_cor_strains)<-str_extract(ECKXXXXremoved,"[A-Z]{3,4}") ###This is a little dangerous since whether the names are actual gene names cannot be verified. Will verify after mapping to gene names in the .GAF file

rm(ECKXXXXremoved)

colnames(cp_cor_strains)<-rownames(cp_cor_strains)
strainCC_gene_names_only<-cp_cor_strains
rm(cp_cor_strains)