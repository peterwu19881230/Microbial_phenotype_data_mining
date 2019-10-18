#Repeat 25%_shuff multiple times to confirm the ROC-AUC result

#Object prep

#==========================================================================================



#source the functions to be used
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste(currentDir,"mis_annotation_functions.R",sep="/"))



createMisAnnotSet_modified=function(typeOfAnnotation){
  id_annot=id_allAttributes[,c("ids",typeOfAnnotation)] %>% unique
  
  shuffledAnnot_list=shuffleAnnotation(id_annot,fractionOfMisAnnot=rep(0.25,10))
  removedAnnot_list=removeAnnotation(id_annot,fractionOfMisAnnot=rep(0.25,10))
  
  return(list(shuffledAnnot_list=shuffledAnnot_list,
              removedAnnot_list=removedAnnot_list))
}


pwyMisAnnot=createMisAnnotSet_modified("Pwy") #pwy
pcomplexMisAnnot=createMisAnnotSet_modified("pcomplex") #pcomplex
operonMisAnnot=createMisAnnotSet_modified("regulator") #operon
regulonMisAnnot=createMisAnnotSet_modified("operon") #regulon
keggMisAnnot=createMisAnnotSet_modified("kegg_modules") #kegg


misAnnotation_all_shuff25percent=list(pwyMisAnnot=pwyMisAnnot,
                       pcomplexMisAnnot=pcomplexMisAnnot,
                       operonMisAnnot=operonMisAnnot,
                       regulonMisAnnot=regulonMisAnnot,
                       keggMisAnnot=keggMisAnnot)



#Function to create strain1_strain2_coAnnotation
#================================================================================================
pairwiseCoannotation=function(annot_id){
  
  attribute_list=attr_list(annot_id[[2]],annot_id[[1]])
  
  uniqueAttr_vec=unlist(attribute_list) %>% unique
  uniqueAttr_vec=uniqueAttr_vec[!is.na(uniqueAttr_vec)]
  
  
  idRow_attrCol=sapply(attribute_list,FUN = function(attribute){ #idRow_attrCol will be created as a matrix
    uniqueAttr_vec %in% attribute  #c("A","B","C") %in% NA will give FALSE, so no worries here
  }) %>% t
  
  
  ##This gets pairwise comparison of having the same annotations or not
  coAnnotationMatrix=idRow_attrCol %*% t(idRow_attrCol)
  coAnnotationMatrix=ifelse(coAnnotationMatrix>=1,1,0) #The matrix structure will be preserved even after using ifelse()
  
  colnames(coAnnotationMatrix)=rownames(coAnnotationMatrix) #Have to manually give the colnames (which are identical to the row names)
  
  
  coAnnotated_table=coAnnotationMatrix %>% as.dist %>% melt_dist #after this, the result becomes a dataframe
  names(coAnnotated_table)[3]="sameORnot"  #have to correct the colname here
  
  return(coAnnotated_table)
}
#================================================================================================





#create co-annotation tables
start.time=Sys.time()

strain1_strain2_coAnnotationTable_list_levels=list()
for(i in 1:11){ # go through different levels of mis-annotation (including the original)
  
  strain1_strain2_coAnnotationTable_list=list()
  for(j in 1:2){ #go through shuffled annotation set and removed annotation set
    
    for(k in seq(misAnnotation_all_shuff25percent)){ #go through all 5 annotation sets 
      
      strain1_strain2_coAnnotation=pairwiseCoannotation(misAnnotation_all_shuff25percent[[k]][[j]][[i]])  
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
end.time-start.time ##Time difference of 1.373838 hours





for(i in 1:11){
  names(strain1_strain2_coAnnotationTable_list_levels[[i]])=c("shuffled","removed")
}



percentMis="25%"
type=c("shuff","remov")
annot=c("pwy","pcomplex","operon","regulon","kegg")
for(i in 1:11){
  for(j in 1:2){
    for(k in 1:5){
      
      if(i==1){
        names(strain1_strain2_coAnnotationTable_list_levels[[i]][[j]])[k+2]=paste(c(type[j],"original",annot[k]),
                                                                                  collapse="_")
      }else{
        names(strain1_strain2_coAnnotationTable_list_levels[[i]][[j]])[k+2]=paste(c(type[j],percentMis,annot[k],i-1),
                                                                                  collapse="_")
      }
      
      
    }
  }
  
}

##(When run on my PC) R will complain not having enough memory at the following steps. What I did was I: save this obj -> restart my PC -> load this obj and continue
#save(strain1_strain2_coAnnotationTable_list_levels,file="Data/temp.RData")
#load("Data/temp.RData")


###bind all results into 1 dataframe
start.time=Sys.time()

for(i in 1:11){
  for(j in 1:2){ 
    
    if(i==1 & j==1 ){
      strain1_strain2_rep25percentMisAnnot=strain1_strain2_coAnnotationTable_list_levels[[i]][[j]]
    }else{
      strain1_strain2_rep25percentMisAnnot=merge(strain1_strain2_rep25percentMisAnnot,strain1_strain2_coAnnotationTable_list_levels[[i]][[j]],
                                            by=c("Object 1","Object 2"))
      
    } 
    
    print(names(strain1_strain2_coAnnotationTable_list_levels[[i]][[j]]))
    
  }
  
}

end.time=Sys.time()
end.time-start.time #Time difference of 58.02344 mins (on my Mac)


str(strain1_strain2_rep25percentMisAnnot)
names(strain1_strain2_rep25percentMisAnnot)[1:2]=c("strain1","strain2")


#save(strain1_strain2_rep25percentMisAnnot,file="Data/strain1_strain2_rep25percentMisAnnot.RData")
#==========================================================================================


load("Data/strain1_strain2_rep25percentMisAnnot.RData")
str(strain1_strain2_rep25percentMisAnnot)


##calculate cumulative sum of "in all 5 annotation sets"
load("Data/strain1strain2_allDistances.RData")
str(strain1strain2_allDistances)


strain1_strain2_rep25percentMisAnnot_allDistance=merge(strain1_strain2_rep25percentMisAnnot,strain1strain2_allDistances ,by=c("strain1","strain2"))
str(strain1_strain2_rep25percentMisAnnot_allDistance)


save(strain1_strain2_rep25percentMisAnnot_allDistance,file="strain1_strain2_rep25percentMisAnnot_allDistance.RData")



load("strain1_strain2_rep25percentMisAnnot_allDistance.RData")


##calculate metrics based on confusion matrix
percentMis=c("original",
             paste("25%_[a-z]{3,}_",1:10,"$",sep="")) 

type=c("shuff","remov")

df=strain1_strain2_rep25percentMisAnnot_allDistance

start.time=Sys.time()

col_names=c()
for(i in 1:11){
  for(j in 1:2){
    coAnnotColumns=df[,grepl(percentMis[i],names(df)) & grepl(type[j],names(df))]
    coAnnotColumn=ifelse(rowSums(coAnnotColumns)==5,1,0)
    
    
    if(i==1 & j==1){
      inAll5_column=coAnnotColumn
    }else{
      inAll5_column=cbind(inAll5_column,coAnnotColumn)
    }
    
    col_names=c(col_names,paste(c(percentMis[i],type[j]),collapse="_"))
    
  }
  
}

colnames(inAll5_column)=col_names


end.time=Sys.time()
end.time-start.time #Time difference of 1.302955 mins


inAll5_column_ordered=inAll5_column[order(strain1_strain2_rep25percentMisAnnot_allDistance$pcc),] 

apply(inAll5_column_ordered,2,sum) 



#plot ROC

cumsums=apply(inAll5_column_ordered,2,cumsum)


conf_metrics_list=list() #If I simply use this, there will be memory problem (even on my Mac): conf_metrics_list=apply(cumsums,2,confusionMatrix_metrics)
for(i in seq(dim(cumsums)[2])){ 
  conf_metrics_list[[i]]=confusionMatrix_metrics(cumsums[,i])
}
names(conf_metrics_list)=colnames(cumsums)


plot_ROC_25percent=function(samples,alpha,size,conf_metrics_list,legend_title="Mis-annotation set(s)"){
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  p=ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                   legend.text=element_text(size=15),
                                   axis.text.x = element_text(size=15),
                                   axis.text.y=element_text(size=15),
                                   axis.title=element_text(size=20))+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[1]]$specificity[samples],sensitivity=conf_metrics_list[[1]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[1])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[2]]$specificity[samples],sensitivity=conf_metrics_list[[2]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[2])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[3]]$specificity[samples],sensitivity=conf_metrics_list[[3]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[3])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[4]]$specificity[samples],sensitivity=conf_metrics_list[[4]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[4])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[5]]$specificity[samples],sensitivity=conf_metrics_list[[5]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[5])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[6]]$specificity[samples],sensitivity=conf_metrics_list[[6]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[6])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[7]]$specificity[samples],sensitivity=conf_metrics_list[[7]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[7])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[8]]$specificity[samples],sensitivity=conf_metrics_list[[8]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[8])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[9]]$specificity[samples],sensitivity=conf_metrics_list[[9]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[9])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[10]]$specificity[samples],sensitivity=conf_metrics_list[[10]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[10])),size=size)+
    geom_line(data=data.frame(FPR=1-conf_metrics_list[[11]]$specificity[samples],sensitivity=conf_metrics_list[[11]]$sensitivity[samples]),aes(FPR,sensitivity,color=names(conf_metrics_list[11])),size=size)+
    
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagonal line
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text',color=legend_title)
  
  return(p)
}


samples=seq(from=0,to=7914231,by=1989); samples[1]=1
roc_shuffle=plot_ROC_25percent(samples=samples,alpha=0.8,size=1,conf_metrics_list[grepl("shuff",names(conf_metrics_list))])
roc_remove=plot_ROC_25percent(samples=samples,alpha=0.8,size=1,conf_metrics_list[grepl("remov",names(conf_metrics_list))])

roc_shuffle 
roc_remove 




#Calculate AUC (Area under curve) using a package
library(pROC) #ref: https://stackoverflow.com/questions/4903092/calculate-auc-in-r



predictor=strain1_strain2_rep25percentMisAnnot_allDistance$pcc
for(i in (1:dim(inAll5_column)[2])[ apply(inAll5_column,2,sum)!=0 ]){  #only apply to columns that have at least one 1
  
  response=inAll5_column[,i]
  
  start.time = Sys.time()
  
  colnames(inAll5_column)[i] %>% print
  auc(response,predictor) %>% print #Syntax: auc(response,predictor)
  
  end.time = Sys.time()
  end.time - start.time #Time difference of 22.89565 secs
  
}

#Result:
'
[1] "original_shuff"
Area under the curve: 0.9131
[1] "original_remov"
Area under the curve: 0.9131
[1] "25%_[a-z]{3,}_1$_shuff"
Area under the curve: 0.9599
[1] "25%_[a-z]{3,}_1$_remov"
Area under the curve: 0.9166
[1] "25%_[a-z]{3,}_2$_shuff"
Area under the curve: 0.8813
[1] "25%_[a-z]{3,}_2$_remov"
Area under the curve: 0.9645
[1] "25%_[a-z]{3,}_3$_shuff"
Area under the curve: 0.9196
[1] "25%_[a-z]{3,}_3$_remov"
Area under the curve: 0.9477
[1] "25%_[a-z]{3,}_4$_shuff"
Area under the curve: 0.8193
[1] "25%_[a-z]{3,}_4$_remov"
Area under the curve: 0.9468
[1] "25%_[a-z]{3,}_5$_shuff"
Area under the curve: 0.8881
[1] "25%_[a-z]{3,}_5$_remov"
Area under the curve: 0.9752
[1] "25%_[a-z]{3,}_6$_shuff"
Area under the curve: 0.8858
[1] "25%_[a-z]{3,}_6$_remov"
Area under the curve: 0.8964
[1] "25%_[a-z]{3,}_7$_shuff"
Area under the curve: 0.9619
[1] "25%_[a-z]{3,}_7$_remov"
Area under the curve: 0.8683
[1] "25%_[a-z]{3,}_8$_shuff"
Area under the curve: 0.8491
[1] "25%_[a-z]{3,}_8$_remov"
Area under the curve: 0.9273
[1] "25%_[a-z]{3,}_9$_shuff"
Area under the curve: 0.9422
[1] "25%_[a-z]{3,}_9$_remov"
Area under the curve: 0.9256
[1] "25%_[a-z]{3,}_10$_shuff"
Area under the curve: 0.9086
[1] "25%_[a-z]{3,}_10$_remov"
Area under the curve: 0.761
'

