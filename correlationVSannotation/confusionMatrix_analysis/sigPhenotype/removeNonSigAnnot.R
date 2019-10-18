#Plot the above result using a fixed pcc. Plot: Coverage VS. accuracy + Coverage VS. sensitivity




## Remove annotations that don't have any significant phenotypes by assigning NA, and return id_attr dataframe 
#------------------------------------------------------------------------------------
removeAnnot=function(attr){
  id_annot=id_allAttributes[,c("ids",attr)] %>% unique
  annotLeft=id_allAttributes[[attr]][id_allAttributes$AnySigFitness & !is.na(id_allAttributes[[attr]])] %>% unique
  allAnnot=id_annot[[attr]][!is.na(id_annot[[attr]])] %>% unique
  annotToBeRemoved=allAnnot[!(allAnnot %in% annotLeft)] 
  
  id_annotLeft=id_annot
  id_annotLeft[[attr]][id_annotLeft[[attr]] %in% annotToBeRemoved]=NA #Note: NA %in% "A" gives FALSE
  
  id_annotLeft=id_annotLeft %>% unique #This is a must-do because assigning new NAs will give duplicated rows
  
  return(id_annotLeft)
}

id_pwy=removeAnnot("Pwy")
id_pcomplex=removeAnnot("pcomplex")

#creating id_pwyANDpcomplex is more complicated:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
id_pwyANDpcomplex=inner_join(id_pwy,id_pcomplex,by="ids")
id_pwyANDpcomplex$pwyANDpcomplex=apply(id_pwyANDpcomplex[,c("Pwy","pcomplex")],1,FUN=function(Row) paste(Row,collapse=" & "))

id_pwyANDpcomplex$pwyANDpcomplex[is.na(id_pwyANDpcomplex$Pwy) | is.na(id_pwyANDpcomplex$pcomplex) ]=NA 
##if at least one type of annotation contains NA, I will assign NA for the intersection of pathway & pcomplex


id_pwyANDpcomplex=id_pwyANDpcomplex[,c("ids","pwyANDpcomplex")] %>% unique #unique() is a must-do because assigning new NAs will give duplicated rows
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#creating id_allAnnot is more complicated:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
id_regulon=removeAnnot("regulator")
id_operon=removeAnnot("operon")
id_kegg_module=removeAnnot("kegg_modules")

id_allAnnot=inner_join(id_pwyANDpcomplex,id_regulon,by="ids")
id_allAnnot=inner_join(id_allAnnot,id_operon,by="ids")
id_allAnnot=inner_join(id_allAnnot,id_kegg_module,by="ids")  

id_allAnnot$allAnnot=apply(id_allAnnot[,c("pwyANDpcomplex","regulator","operon","kegg_modules")],1,FUN=function(Row) paste(Row,collapse=" & "))


id_allAnnot$allAnnot[is.na(id_allAnnot$pwyANDpcomplex) | 
                       is.na(id_allAnnot$regulator) |
                       is.na(id_allAnnot$operon) |
                       is.na(id_allAnnot$kegg_modules)
                     ]=NA  ##if at least one type of annotation contains NA, I will assign NA for the intersection 

id_allAnnot=id_allAnnot[,c("ids","allAnnot")] %>% unique #unique() is a must-do because assigning new NAs will give duplicated rows
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------------------------------------------------------



