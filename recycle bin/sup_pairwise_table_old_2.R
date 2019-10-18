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






#Parallel computation: divide the id-annotation table -> compute pairwise co-annotation tables -> join
library(foreach)
library(doParallel)
library(parallel)
library(doSNOW) 

table_list=list(id_allAnnots[,1:1000],
                id_allAnnots[,1001:2000],
                id_allAnnots[,2001:3000],
                id_allAnnots[,3001:4000],
                id_allAnnots[,4001:5000],
                id_allAnnots[,5001:6000],
                id_allAnnots[,6001:7000],
                id_allAnnots[,7001:7439]
)

start.time=Sys.time()

numCores=detectCores()-1
cl=makeCluster(numCores)
registerDoSNOW(cl)
clusterExport(cl,"table_list") 

pairwise_table_list=foreach(i=1:length(table_list)) %dopar%{
  require(dplyr)
  require(reshape2)
  require(stringr)
  
  table=table_list[[i]]
  
  result=NULL
  for(j in 1:dim(table)[2]){
    X=t(table[,j])
    out=t(X)%*%X %>% melt

    #============================This if-else causes the runtime to change from 24 min on my PC to forever==================================        
    if(j==1){
      result=out[,1:2]
      name=ifelse(out[,3],colnames(table)[j],"")
      result=cbind(result,name)
    }else{
      name=ifelse(out[,3],colnames(table)[j],"")
      result[,3]=paste(result[,3],name,sep=";")
      result[,3]=str_remove(result[,3], ";$") #clean ";" if co-annotation field is "" (empty)
    }
  }
    #==============================================================
  
  result[,3]=str_remove(result[,3], "^;") #clean ";" at the beginning (It was formed when this happends: "" +  annot )
  
  return(result)
}


stopCluster(cl)

end.time=Sys.time()
end.time-start.time #Time difference of ? on my PC


##turn pairwise_table_list into a table: pairwise_table

pairwise_table=Reduce("cbind",pairwise_table_list)

colnames(pairwise_table)=c("strain1","strain2",colnames(table))

temp=apply(pairwise_table,1,FUN=function(row_){
  co_annotations_row=colnames(id_allAnnots)[row_[-c(1,2)]]
  paste(co_annotations_row,sep=";")
})

#rm(pairwise_table)



#merge tab1, tab2 by left_join




#save the result to Data/sup_pairwise_table.txt