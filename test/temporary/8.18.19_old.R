#(!) The code here should be executable but computing MI takes lots of time, as apposed to pcc. So I will recode.

#Generate corr VS annot graph based on these parameters:
#1. Ternary data using different FDRs
#2. different annotations
#3. differnt metrics calculated from confusion matrix


#specify parameters
#====================================================================================================
#1. load various data using different FDR
load("Data/ternary_data_various_FDR.RData") #This gives: 
##     Ternary_Data_324cutff_NAremoved_1_FDR,
##     Ternary_Data_324cutff_NAremoved_10_FDR,
##     Ternary_Data_324cutff_NAremoved_15_FDR,
##     Ternary_Data_324cutff_NAremoved_20_FDR,

##Ternary_Data_324cutff_NAremoved is already sourced from Nichols_preload.R. This one uses 5% FDR

#2.
annots=c("Pwy","pcomplex","operon","regulator","kegg_modules")

#3.
metrics=c("precision","sensitivity")
#====================================================================================================


#prep data
#====================================================================================================
my_mutual_info_dist=function(dat){
  mi=infotheo::mutinformation(t(dat) %>% as.data.frame) %>% natstobits
  mi_distObj=dist(mi)
  mi_distance=1-mi_distObj
  return(mi_distance)
}

start.time=Sys.time()



result_list=list()
i=1
for(dat in list(Ternary_Data_324cutff_NAremoved_1_FDR,
                Ternary_Data_324cutff_NAremoved,
                Ternary_Data_324cutff_NAremoved_10_FDR,
                Ternary_Data_324cutff_NAremoved_15_FDR,
                Ternary_Data_324cutff_NAremoved_20_FDR
)){
  
  for(annot in annots){
    
    id_annot=id_allAttributes[,c("ids",annot)] %>% unique()
    attribute_list=attr_list(id_annot[,1],id_annot[,2])
    
    result_list[[i]]=dist_TF_cumsum_matirxOperation(dat,attribute_list=attribute_list,dist_metric=my_mutual_info_dist)
    print(paste(i,"th result is done",sep=""))
    
    i=i+1
    
  }
  
}

end.time=Sys.time()
end.time-start.time

#(!)I might have to do any using a different way

#====================================================================================================


#get graphs
#====================================================================================================



#====================================================================================================



