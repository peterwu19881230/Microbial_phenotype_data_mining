#Dimension reduction for annotations

##========================================================================================================================================

#Get dummy variables for all annotations (rows: strains, columns: whether the strain has that particular annotation or not)

annots=c("Pwy","pcomplex","regulator","operon","kegg_modules")

annot_dummies=matrix(, nrow = 3979, ncol = 0) #I missed the 1st argument on purpose
col_names=c()
for(annot in annots){
  id_annot=id_allAttributes[,c("ids",annot)] %>% unique()
  attribute_list=attr_list(id_annot[[1]],id_annot[[2]])
  
  uniqueAttr_vec=unlist(attribute_list) %>% unique
  uniqueAttr_vec=uniqueAttr_vec[!is.na(uniqueAttr_vec)]
  
  idRow_attrCol=sapply(attribute_list,FUN = function(attribute){ #idRow_attrCol will be created as a matrix
    uniqueAttr_vec %in% attribute  #c("A","B","C") %in% NA will give FALSE, so no worries here
  }) %>% t
  
  annot_dummies=cbind(annot_dummies,idRow_attrCol)
  col_names=c(col_names,uniqueAttr_vec)
}

colnames(annot_dummies)=col_names

annot_dummies=annot_dummies[order(as.numeric(rownames(annot_dummies))),]

save(annot_dummies,file="Data/annot_dummies.RData")

annot_dummies[]=as.numeric(annot_dummies)
write.csv(annot_dummies,file="Data/annot_dummies.csv",row.names = F)
##========================================================================================================================================


##dim(annot_dummies) 5031 distinct annotations
##hist(rowSums(annot_dummies))




## MCA (similar to PCA but variables are categorical)
###Ref: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/114-mca-multiple-correspondence-analysis-in-r-essentials/
###install.packages("FactoMineR")


## I used the first 200 annotations to test. Doesn't seem to work because the 1st variable only explains 5.4% variance
###mca=FactoMineR::MCA(annot_dummies[,1:200],ncp=1,graph=T)
###factoextra::fviz_screeplot(mca,addlabels = TRUE)


## LDA


## Hierarchical clustering

### based on hamming distance
annot_hamming=hamming_distance(annot_dummies) %>% as.dist #this takes about 1 min

annot_hclust=stats::hclust(annot_hamming)

annot_20_clusters=stats::cutree(annot_hclust,k=20) #I am trying 20 clusters because in Campos' they said by t-SNE there's about 20 clusters
annot_20_clusters=annot_20_clusters[order(as.numeric(names(annot_20_clusters)))]

### This shows how many genes are grouped into each cluster
table(annot_20_clusters)

annot_20_clusters_hclust_hamming=annot_20_clusters
save(annot_20_clusters_hclust_hamming,file="Data/annot_20_clusters_hclust_hamming.RData")
write.table(as.matrix(annot_20_clusters_hclust_hamming),sep=",",file="Data/annot_20_clusters_hclust_hamming.csv",col.names = F,row.names = F)

#======Prepare a list of genes that seem to be more well labeled======
load("Data/annot_20_clusters_hclust_hamming.RData")

annot_20_clusters_subset=annot_20_clusters[annot_20_clusters!=2]
table(annot_20_clusters_subset)
#save(annot_20_clusters_subset,file="Data/annot_20_clusters_subset_hclust_hamming.RData")


#=====================================================================


### based on mutual information
### (!) Have to rerun because I corrected the obj: annot_dummies
start.time=Sys.time()
annot_mi= my_mutual_info_dist(annot_dummies)  
end.time=Sys.time()
end.time-start.time #Time difference of 19.22066 hours on my PC
#save(annot_mi,file="Data/pairwise_mi_for_Nichols_annot.RData") #the above takes extremely long so I saved it first


annot_hclust=stats::hclust(annot_mi)

annot_20_clusters=stats::cutree(annot_hclust,k=20) #I am trying 20 clusters because in Campos' they said by t-SNE there's about 20 clusters
annot_20_clusters=annot_20_clusters[order(as.numeric(names(annot_20_clusters)))]

### This shows how many genes are grouped into each cluster
table(annot_20_clusters) #This result shows mi is probably not the right method to do hclust here






# K-mean clustering





##other methods

###CorEx
###https://stats.stackexchange.com/questions/159705/would-pca-work-for-boolean-binary-data-types
###https://arxiv.org/pdf/1410.7404.pdf