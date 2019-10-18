#create mis-annotation sets for all annotation sets I have

#source the functions to be used
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste(currentDir,"mis_annotation_functions.R",sep="/"))



pwyMisAnnot=createMisAnnotSet("Pwy") #pwy
pcomplexMisAnnot=createMisAnnotSet("pcomplex") #pcomplex
operonMisAnnot=createMisAnnotSet("regulator") #operon
regulonMisAnnot=createMisAnnotSet("operon") #regulon
keggMisAnnot=createMisAnnotSet("kegg_modules") #kegg


misAnnotation_all=list(pwyMisAnnot=pwyMisAnnot,
                       pcomplexMisAnnot=pcomplexMisAnnot,
                       operonMisAnnot=operonMisAnnot,
                       regulonMisAnnot=regulonMisAnnot,
                       keggMisAnnot=keggMisAnnot
)

#save(misAnnotation_all,file="Data/misAnnotation_all.RData")







#create co-annotation tables
'
start.time=Sys.time()

strain1_strain2_coAnnotationTable_list_levels=list()
for(i in 1:5){ # go through different levels of mis-annotation (including the original)
  
  strain1_strain2_coAnnotationTable_list=list()
  for(j in 1:2){ #go through shuffled annotation set and removed annotation set
    
    for(k in seq(misAnnotation_all)){ #go through all 5 annotation sets 
      
      strain1_strain2_coAnnotation=pairwiseCoannotation(misAnnotation_all[[k]][[j]][[i]])  
      names(strain1_strain2_coAnnotation)[3]=c("pwy","pcomplex","operon","regulon","kegg")[k]
      
      if(k==1) {strain1_strain2_coAnnotationTable=strain1_strain2_coAnnotation; next} 
      
      strain1_strain2_coAnnotationTable=merge(strain1_strain2_coAnnotationTable,strain1_strain2_coAnnotation,by=c("Object 1","Object 2"))
      
      #######
      print(dim(strain1_strain2_coAnnotationTable)[1]) #should be 7914231 for every iteration
      #######
    }
    
    strain1_strain2_coAnnotationTable_list[[j]]=strain1_strain2_coAnnotationTable
    
  } 
  
  strain1_strain2_coAnnotationTable_list_levels[[i]]=strain1_strain2_coAnnotationTable_list 
  
}

end.time=Sys.time()
end.time-start.time #Time difference of 1.138508 hours


names(strain1_strain2_coAnnotationTable_list_levels)=c("0%","0.5%","5%","25%","50%")

for(i in 1:5){
    names(strain1_strain2_coAnnotationTable_list_levels[[i]])=c("shuffled","removed")
}
 




##reorganize the obj

###rename the coAnnotation columns
percentMis=c("0%", "0.5%", "5%", "25%", "50%")
type=c("shuff","remov")
annot=c("pwy","pcomplex","operon","regulon","kegg")
for(i in 1:5){
  for(j in 1:2){
    for(k in 1:5){
      names(strain1_strain2_coAnnotationTable_list_levels[[i]][[j]])[k+2]=paste(c(type[j],percentMis[i],annot[k]),
                                                                                collapse="_")
    }
  }
  
}

###bind all results into 1 dataframe
start.time=Sys.time()

for(i in 1:5){
  for(j in 1:2){ 
      
      if(i==1 & j==1 ){
        strain1_strain2_misCoannotation=strain1_strain2_coAnnotationTable_list_levels[[i]][[j]]
      }else{
        strain1_strain2_misCoannotation=merge(strain1_strain2_misCoannotation,strain1_strain2_coAnnotationTable_list_levels[[i]][[j]],
                                              by=c("Object 1","Object 2"))
      } 
      
    
  }
  
}

end.time=Sys.time()
end.time-start.time


str(strain1_strain2_misCoannotation)
names(strain1_strain2_misCoannotation)[1:2]=c("strain1","strain2")


#save(strain1_strain2_misCoannotation,file="Data/strain1_strain2_misCoannotation.RData")
'










