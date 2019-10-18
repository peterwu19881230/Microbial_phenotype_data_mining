
#All_Data_NAimputed & Ternary_Data_324cutff_NAremoved are used in the following algorithm

#Data prep
#================================================================================================
# strains_in_sub_profile=list()
# for(i in 1:3979){
#   strain=All_Data_NAimputed[i,]
#   print(i) #print the strain id (=index)
#   
#   severestCondIndex=which(strain==max(strain))
#   
#   strain_andRelatedStrains=which(Ternary_Data_324cutff_NAremoved[,severestCondIndex]!=0)
#   attributes(strain_andRelatedStrains)$names=NULL #remove the vector name (since it's identical to the vector value)
#   
#   strains_in_sub_profile=c(strains_in_sub_profile,list(strain_andRelatedStrains))
#   
#   }
# 
# 
# dir.create("Data/severestFitness_iterative")
# save(strains_in_sub_profile,file="Data/severestFitness_iterative/strains_in_sub_profile.RData")
#================================================================================================

load("Data/severestFitness_iterative/strains_in_sub_profile.RData")

##Look at the distribution of no. of strains for sub-profiles
strain_leng=sapply(strains_in_sub_profile,length)
hist(strain_leng)
summary(strain_leng)


##average pccs for each cluster of genes

calculate_avg_pcc=function(strain_cluster_df){
  pccs=melt_dist(as.dist(cor(t(strain_cluster_df))))[,3] #melt_dist is defined in functions.R
  abs_pccs=abs(pccs)
  return(mean(abs_pccs))
}


start.time=Sys.time()
avg_pcc_ForEachStrain=sapply(strains_in_sub_profile,FUN=function(strain_ids){
  strain_cluster_df=All_Data_NAimputed[strain_ids,]
  result=calculate_avg_pcc(strain_cluster_df)
  
  return(result)
})
end.time=Sys.time()
end.time-start.time #Time difference of 28.06992 secs


hist(avg_pcc_ForEachStrain)
summary(avg_pcc_ForEachStrain)


##coannotation enrichment for each cluster
all_strain_pairs=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2")]
calculate_average_coannotation=function(strains){ #input should be strain indices
  
  if(length(strains) %in% c(0,1) ){ #when there is only 1 significant phenotype under that condition
    coannotations_average_for_pairs=rep(NA,5)
    names(coannotations_average_for_pairs)=c("Pwy","pcomplex","operon","regulator","kegg_modules")
    return(coannotations_average_for_pairs) 
  }else{
    
    combs=t(combn(strains,2))
    
    TFindices=(all_strain_pairs[,1] %in% combs[,1]) & (all_strain_pairs[,2] %in% combs[,2]) 
    ##The above line is based on: t(combn()) always give the smaller number on the first column
    
    coannotations_df=strain1strain2_allAnnotations_allDistances[TFindices,c("Pwy","pcomplex","operon","regulator","kegg_modules")]
    
    coannotations_average_for_pairs=apply(coannotations_df,2,sum)/dim(combs)[1]
    
    return(coannotations_average_for_pairs) #coannotations_average_for_pairs is a named vector
    
    
  }
  
  
}  

#===================This block of code is too slow (But it should run. I have tested it) => jump to the next approach (start by conditions)
# library(pbapply)
# start.time=Sys.time()
# enrichedAnnotation_for_each_cluster=pbsapply(strains_in_sub_profile[1],calculate_average_coannotation)
# end.time=Sys.time()
# end.time-start.time
#====================


#Another approach: start with conditions and look for enrichment of coannotations

##Use the function defined above: calculate_average_coannotation()

##for positive phenotypes
positive_phenotype_matrix=(Ternary_Data_324cutff_NAremoved==1)

library(pbapply)
start.time=Sys.time()
coannotation_positive=pbapply(positive_phenotype_matrix,2,FUN=function(condition){
  strains=which(condition==T)
  result=calculate_average_coannotation(strains)
  return(result)
})
end.time=Sys.time()
end.time-start.time #Time difference of 1.794377 mins


##for negative phenotypes
negative_phenotype_matrix=(Ternary_Data_324cutff_NAremoved==-1)

start.time=Sys.time()
coannotation_negative=pbapply(negative_phenotype_matrix,2,FUN=function(condition){
  strains=which(condition==T)
  result=calculate_average_coannotation(strains)
  return(result)
})
end.time=Sys.time()
end.time-start.time #Time difference of 1.807359 mins


#negative control (random)
all_coannotation=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")]
random_coannotation=apply(all_coannotation,2,sum)/7914231






#########
co_annot_table=coannotation_negative #change this to coannotation_positive or coannotation_negative and get the other half of the result

result=list()
for(annot in c("Pwy","pcomplex","operon","regulator","kegg_modules")){
  
  
  enrichment=co_annot_table[annot,]/random_coannotation[annot]
  names(enrichment)=colnames(All_Data_NAimputed)
  
  result[[annot]]$summary=(enrichment) %>% summary 
  result[[annot]]$cond_enrichment_sorted=sort(enrichment,decreasing=T) #NA is ignored in sort()
  
}


table_first_thirty=lapply(result,FUN=function(annot_data){
  
  table_=round(annot_data$cond_enrichment_sorted[1:30],digits=2)
  data.frame(fold_enrichment=table_)
  
})
#########


#summary google sheet:
##https://docs.google.com/spreadsheets/d/19XIbsLOdql3wc1PamnpXw0oR47VGqcmogNmx9az3tEk/edit#gid=0





