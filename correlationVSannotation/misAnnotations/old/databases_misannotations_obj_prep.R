
load("Data/fake_pwyAnnot_list.RData")

dat=All_Data_NAimputed


#This is a completely random id_pwy 
if(!file.exists("Data/mis_annotation_random.RData")){  #If the file cannot be found, re-process
  
  ##Get the directory where the running script is at:
  ##https://stackoverflow.com/questions/47044068/get-the-path-of-current-script/47045368
  script_path=dirname(rstudioapi::getSourceEditorContext()$path) 
  source(paste(script_path,"/random_id_pwy_annot.R",sep=""))
  
  start.time = Sys.time()
  dat=All_Data_NAimputed
  attribute_list=attr_list(bad_id_pwy_annotations[,1],bad_id_pwy_annotations[,2])
  mis_annotation_random=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = pcc_dist)
  
  save(mis_annotation_random,file="Data/mis_annotation_random.RData")
  
  
  end.time = Sys.time()
  end.time - start.time #Time difference of 6.900174 mins
  
}



##For quantitative results

n=c(0,1,10,100,500,1000,1500,2000) #no. of misannotations. 0 stands for the original annotation set


#Mis-annotate by shuffling annotations
if(!file.exists("Data/mis_annotation_quantitative_list.RData")){
start.time = Sys.time()

mis_annotation_quantitative_list=list()
for(i in seq(n)){
  
  attribute_list=attr_list(fake_pwyAnnot_list[[i]][[2]],fake_pwyAnnot_list[[i]][[1]]) 
  
  name=paste("mis_annotation_",n[i],sep="")
  
  mis_annotation_quantitative_list[[name]]$pearson=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = pcc_dist)
  mis_annotation_quantitative_list[[name]]$spearman=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = spearman_dist)
  mis_annotation_quantitative_list[[name]]$euclidean=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = euclidean_dist)

}
save(mis_annotation_quantitative_list,file="Data/mis_annotation_quantitative_list.RData")
  
end.time = Sys.time()
end.time - start.time #Time difference of 3.14643 hours from  my pc

}


#Mis-annotate by leaving out annotations

load("Data/annotThrown_pwyAnnot.RData")

n=c(0,1,10,100,500,1000,1500,2000) #no. of misannotations. 0 stands for the original annotation set

if(!file.exists("Data/thrown_annotation_quantitative_list.RData")){
  start.time = Sys.time()
  
  thrown_annotation_quantitative_list=list()
  for(i in seq(n)){
    
    attribute_list=attr_list(annotThrown_pwyAnnot$ids,annotThrown_pwyAnnot[[i+1]]) 
    
    name=paste("thrown_annotation_",n[i],sep="")
    
    thrown_annotation_quantitative_list[[name]]$pearson=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = pcc_dist)
    thrown_annotation_quantitative_list[[name]]$spearman=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = spearman_dist)
    thrown_annotation_quantitative_list[[name]]$euclidean=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = euclidean_dist)
    
  }
  save(thrown_annotation_quantitative_list,file="Data/thrown_annotation_quantitative_list.RData")
  
  end.time = Sys.time()
  end.time - start.time 

}


#Some additional object I want to create for some testing

'

start.time = Sys.time()
indexNotNA=(1:dim(annotThrown_pwyAnnot)[1])[!is.na(annotThrown_pwyAnnot$Pwy)]

set.seed(6)
indexToBeThrown=sample(indexNotNA,1500)

new_annot=annotThrown_pwyAnnot$Pwy
new_annot[indexToBeThrown]=NA

attribute_list=attr_list(annotThrown_pwyAnnot$ids,new_annot)


diff_1500_original=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 5.044517 mins



start.time = Sys.time()
indexNotNA=(1:dim(annotThrown_pwyAnnot)[1])[!is.na(annotThrown_pwyAnnot$Pwy)]

set.seed(101)
indexToBeThrown=sample(indexNotNA,1500)

new_annot=annotThrown_pwyAnnot$Pwy
new_annot[indexToBeThrown]=NA

attribute_list=attr_list(annotThrown_pwyAnnot$ids,new_annot)


diff_1500=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = pcc_dist)
end.time = Sys.time()
end.time - start.time #Time difference of 5.058408 mins



##2.  A list that has original annotations + different no. of annotations left out
class(id_allAttributes)
id_pwyAnnot=id_allAttributes[,c("ids","Pwy")] %>% unique
dim(id_pwyAnnot)

n=c(1,1:20*100)

annotThrown_pwyAnnot_2=id_pwyAnnot
indexNotNA=(1:dim(annotThrown_pwyAnnot_2)[1])[!is.na(annotThrown_pwyAnnot_2$Pwy)]

for(i in seq(n)){

set.seed(i)
indexToBeThrown=sample(indexNotNA,n[i])

new_annot=annotThrown_pwyAnnot_2$Pwy
new_annot[indexToBeThrown]=NA

annotThrown_pwyAnnot_2[[paste("less_pwy_",i,sep="")]]=new_annot

#confirm the no. of annotations thrown out
original=annotThrown_pwyAnnot_2$Pwy
original[is.na(original)]="NA"

later=new_annot
later[is.na(later)]="NA"

sum(original!=later) %>% print
}






start.time = Sys.time()

thrown_annotation_quantitative_list_2=list()
for(i in seq(n)){

attribute_list=attr_list(annotThrown_pwyAnnot_2$ids,annotThrown_pwyAnnot_2[[i+1]]) 

name=paste("thrown_annotation_",n[i],sep="")

thrown_annotation_quantitative_list_2[[name]]$pearson=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = pcc_dist)
thrown_annotation_quantitative_list_2[[name]]$spearman=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = spearman_dist)
thrown_annotation_quantitative_list_2[[name]]$euclidean=dist_TF_cumsum(data=dat,attribute_list = attribute_list,dist_metric = euclidean_dist)

}
save(thrown_annotation_quantitative_list_2,file="Data/thrown_annotation_quantitative_list_2.RData")

end.time = Sys.time()
end.time - start.time #Time difference of 10.3031 hours

'



