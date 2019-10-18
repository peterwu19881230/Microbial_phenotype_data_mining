#See how overrepresentation works on phenotype clusters of strains


#I already have this obj sourced, but the original code generating this was lost. 
#Just want to make sure this obj is correct
hclust_pcc_complete=(1-abs(cor_strains)) %>% as.dist %>% hclust(method="complete") 

#Did I define this obj correctly?
k=500
trees=hclust_pcc_complete %>% cutree(k=k)
##no. of clusters =10 is an arbitrary no. I picked

#Find the UniProt ID for each cluster 
UniProtID=list()
for(i in 1:k){
  id=which(trees==i) #Get id for each small tree (Ref: https://www.biostars.org/p/308265/)
  UniProtID[[i]]=id.UniProt$Data[id.UniProt$ids %in% id]
}


#file.create("Data/k_clusters_uniprot.txt") #initiate an empty file
#for(i in 1:k){
#  write(paste(UniProtID[[i]], collapse=" "),"Data/k_clusters_uniprot.txt",append=T)
#}







