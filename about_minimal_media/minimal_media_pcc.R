#All the unique nutrients in minimal media (https://docs.google.com/spreadsheets/d/16ZYwtAqf0Yl8oN4erbLnehFpXGOPbmGmNFG8vhI7F2k/edit#gid=1339220469)
TF=( names(uniqueChemIndex) %in% c("NH4Cl (MOPS)","Iron excess-FeSO4","Iron starvation-FeSO4","Acetate (M9)",
                                   "Glucosamine (M9)","Glucose (M9)","Glycerol (M9)","Maltose (M9)","N-acetyl Glucosamine","Succinate (M9)") )
minimal_uniq=names(uniqueChemIndex)[TF]

cond_indices_minimal=uniqueChemIndex[minimal_uniq] %>% as.numeric

#experiments:

##do corr VS annot excluding the minial media condition

##(I suspect this will give us a better result. There are only 10 minimal media conditions out of 324 conditions)

##correlate genes only based on the minimal media condition
cor_matrix=cor(t(All_Data_NAimputed[,cond_indices_minimal])) ##warning: sd is 0 => might be scores for some strains are all 0
##anyIncomplete(cor_matrix) #verify that NAs have been generated
cor_matrix[is.na(cor_matrix)]=-9 #change where cor=NA to -9 
dist_=as.dist(cor_matrix)
minimal_cor_table=meltANDsort_dist(dist_)
minimal_cor_table=minimal_cor_table[!(minimal_cor_table[,3]==-9),] 

#distribution of the PCCs
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
png(filename=paste(currentDir,"minimal_cor.png",sep="/"))
hist(minimal_cor_table[,3],xlab="PCC",main="")
dev.off()



ids_genes=id_allAttributes[,c("ids","associated_gene_names")] %>% unique
#sum(duplicated(ids_genes)) #does any id match more than 1 gene names => no => I can use match() without worries

minimal_cor_table$gene1=ids_genes$associated_gene_names[match(minimal_cor_table$`Object 1`,ids_genes$ids)]
minimal_cor_table$gene2=ids_genes$associated_gene_names[match(minimal_cor_table$`Object 2`,ids_genes$ids)]

minimal_cor_table=minimal_cor_table[,c(4,5,3)]
names(minimal_cor_table)=c("gene1","gene2","pcc")
minimal_cor_table$pcc=round(minimal_cor_table$pcc,2)

currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
#write.table(minimal_cor_table,file="Data/minimal_cor_table.txt",sep="\t",quote=F,row.names = F)



#distribution of the fitness scores given the minimal media
all_fitness=All_Data_NAimputed %>% as.matrix %>% as.numeric
fitness_minimal_media=All_Data_NAimputed[,cond_indices_minimal] %>% as.matrix %>% as.numeric


plot(density(all_fitness),col="blue")
lines(density(fitness_minimal_media),col="red")








