#Generate sup table that has: strain1 - strain2 - pcc - ternary MI - co-annotatinos (separated by ";")

tab1=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","pcc","mi_ternary")] %>% as.matrix
tab2=id_allAttributes[,c("ids","Pwy","pcomplex","regulator","operon","kegg_modules")] %>% as.matrix

#from tab2 find co-annotations for all pairs. Separate co-annotations by ";" when a pair of genes have more than 1 annotation
all_annot=tab2[,-1] %>% as.vector %>% unique



id_annot_list=list()
for(id in unique(tab2[,1])){
  annot=tab2[tab2[,1]==id,-1] %>% as.vector %>% unique
  annot=annot[!is.na(annot)]
  id_annot_list[[id]]=annot
}

temp=sapply(id_annot_list,FUN=function(annots){
  all_annot %in% annots
})

id_allAnnots=t(temp)
colnames(id_allAnnots)=all_annot


start.time=Sys.time()

table=id_allAnnots[,1:100] #remove [,1:X] after test is complete
pairwise_table_list=list()
for(i in 1:dim(table)[2]){
  
  X=t(table[,i])
  out=t(X)%*%X %>% melt
  
  if(i==1){
    pairwise_table_list[[i]]=out
  }else{
    pairwise_table_list[[i]]=out[,3]
  }
  
}

end.time=Sys.time()
end.time-start.time #Time difference of 54.77922 secs

##(!)Error: vector memory exhausted (limit reached?)
##add code to turn pairwise_table_list into a table: pairwise_table

colnames(pairwise_table)=c("strain1","strain2",colnames(table))

temp=apply(pairwise_table,1,FUN=function(row_){
  id_allAnnots[row_]
})

rm(pairwise_table)



#merge tab1, tab2 by left_join

#save the result to Data/sup_pairwise_table.txt