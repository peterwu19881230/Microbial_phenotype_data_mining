#Generate sup table that has: strain1 - strain2 - pcc - ternary MI - co-annotatinos (separated by ";")


tab1=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","pcc","mi_ternary")] %>% as.matrix

##for each kind of annotation get a coannotation column
#====================================================
id_allCoannotationColumns=matrix(NA,nrow=7914231,ncol=5); count=1
for(annot in c("Pwy","pcomplex","regulator_name","operon_name","kegg_modules")){
  tab2=id_allAttributes[,c("ids",annot)] %>% as.matrix
  
  
  #from tab2 find co-annotations for all pairs. Separate co-annotations by ";" when a pair of genes have more than 1 annotation
  all_annot=tab2[,-1] %>% as.vector %>% unique %>% removeNA
  
  
  
  id_annot_list=list()
  for(id in unique(tab2[,1])){
    annot=tab2[tab2[,1]==id,-1] %>% as.vector %>% unique %>% removeNA
    id_annot_list[[id]]=annot
  }
  
  temp=sapply(id_annot_list,FUN=function(annots){
    all_annot %in% annots
  })
  
  id_allAnnots=t(temp)
  colnames(id_allAnnots)=all_annot
  
  
  
  id_allAnnotsInOneCol=apply(id_allAnnots,1,FUN=function(row_TF){
    paste(all_annot[row_TF],collapse=";")
  })
  
  id_allAnnotsInOneCol=id_allAnnotsInOneCol[order(names(id_allAnnotsInOneCol) %>% as.numeric)] #sort by the strain ids
  
  
  id_allCoannotationColumns[,count]=id_allAnnotsInOneCol 
  
  cat(paste0("loop ",count," is completed "))
  count=count+1
}

#====================================================




#get the first 2 columns of the result table
strain1strain2=t(combn(3979,2)) #now the strain ids are type integer


#get the 3rd column of the result table
#====================================================
annot_columns=cbind(id_allAnnotsInOneCol[strain1strain2[,1]],id_allAnnotsInOneCol[strain1strain2[,2]]) #get annotations for strain1 and strain2

start.time=Sys.time()
third_col=apply(annot_columns,1,FUN=function(annot1_annot2){
  annot1=strsplit(annot1_annot2[1],split=";")[[1]] #strsplit returns the result in the form of a list, so I have to convert back to vector
  annot2=strsplit(annot1_annot2[2],split=";")[[1]]
  paste(intersect(annot1,annot2),collapse=";")
})
end.time=Sys.time()
end.time-start.time #Time difference of 15.93562 mins
#====================================================



strain1_strain2_third_col=cbind(strain1strain2,third_col)
colnames(strain1_strain2_third_col)=c("strain1","strain2","co_annotation")

sup_pairwise_table=merge(as.data.frame(tab1),as.data.frame(strain1_strain2_third_col),by=c("strain1","strain2"))

##check if the merging doesn't fail because of orders of strain1 - strain2
##anyImcomplete(sup_pairwise_table) #no NAs are found. The merging should be ok


##manually check the some of the co-annotation results

#function to check coannotation of 2 strains
check_coannotation=function(strain1,strain2){ #strain1, strain2 are both integers
  strain1_annot=id_allAttributes[id_allAttributes$ids==strain1,c("Pwy","pcomplex","regulator","operon","kegg_modules")] %>% as.vector %>% removeNA %>% unique
  strain2_annot=id_allAttributes[id_allAttributes$ids==strain2,c("Pwy","pcomplex","regulator","operon","kegg_modules")] %>% as.vector %>% removeNA %>% unique
  
  return(intersect(strain1_annot,strain2_annot))
}


###1271	3896	0.87588503	1.0000000	"TSR-CPLX;TSR-GLUME;TSR-GLN;TSR-GLU;ECK120000017"
sup_pairwise_table[sup_pairwise_table$strain1==1271&sup_pairwise_table$strain2==3896,"co_annotation"]
check_coannotation(1271,3896)

###1377	1568	0.98708375	0.9992402	(no coannotation)
sup_pairwise_table[sup_pairwise_table$strain1==1377&sup_pairwise_table$strain2==1568,"co_annotation"]
check_coannotation(1377,1568)

###1635	1685	0.99625592	0.9998602	ECK120000540
sup_pairwise_table[sup_pairwise_table$strain1==1635&sup_pairwise_table$strain2==1685,"co_annotation"]
check_coannotation(1635,1685)


#replace strain1-strain2 with ECK. If no ECK, use EcoCyc identifier: 3341 3361 3370 (Nichols' strain IDs that don't have ECK) 

##bind ECK for strain1, strain2, respectively
id_ECK=id_allAttributes[,c("ids","ECK")] %>% unique
id_ECK$ids=as.numeric(id_ECK$ids)

temp1=id_ECK; names(temp1)[1]="strain1"
temp2=left_join(sup_pairwise_table,temp1,by="strain1")
temp1=id_ECK; names(temp1)[1]="strain2"
temp2=left_join(temp2,temp1,by="strain2")

##name the strains by EcocycID if there is no associated ECK
id_allAttributes$EcoCycID[id_allAttributes$ids=="3341"] #"G0-8894"
id_allAttributes$EcoCycID[id_allAttributes$ids=="3361"] #"G0-8906"
id_allAttributes$EcoCycID[id_allAttributes$ids=="3370"] #"G0-10202"

temp2$ECK.x[temp2$strain1=="3341"]="G0-8894"
temp2$ECK.x[temp2$strain1=="3361"]="G0-8906"
temp2$ECK.x[temp2$strain1=="3370"]="G0-10202"

##subset and rearrange temp2 to get updated sup_pairwise_table
temp2=temp2[,c("ECK.x","ECK.y","pcc","mi_ternary","co_annotation")]; names(temp2)[1:2]=c("strain1","strain2")

#round pcc and mi_ternary to reduce the file size
temp2$pcc=round(temp2$pcc,2)
temp2$mi_ternary=round(temp2$mi_ternary,2)


#save the result to Data/sup_pairwise_table.txt
##write.table(sup_pairwise_table,file="Data/sup_pairwise_table.txt",quote=F,sep="\t")






