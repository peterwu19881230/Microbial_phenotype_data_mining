#response variables: 
#(a) super-pathway 
#(b) Don't really know how to categorize strains (strains-annotations are many-to-many) 
#=> Create series of annotated characters behind each strain and see how people in NCHU do it

##(a) These 2 lines show why super-pathway might not be a good response variable: too few annotated strains

##how many unique super pathways are there?
id_allAttributes$Pwy[grepl(pattern="super",x=tolower(id_allAttributes$Pwy))] %>% unique # 3 unique super pathways

##How many Nichols' IDs have been annotated with those super-pathways?
id_allAttributes$ids[grepl(pattern="super",x=tolower(id_allAttributes$Pwy))] %>% unique %>% length #35 (Not many)


#(b)
names(id_allAttributes) #These will be used as responses: "Pwy" "pcomplex" "regulator" "operon" "kegg_modules"
class(id_allAttributes$ids) #"character"

#Build a dataframe that each id gets all kinds of annotations 
id_usedAnnots_list=sapply(as.character(1:3979),FUN=function(strain){
  usedAttr=id_allAttributes[id_allAttributes$ids==strain,c("Pwy", "pcomplex", "regulator", "operon", "kegg_modules")] %>% unlist %>% unique
  usedAttr=usedAttr[!is.na(usedAttr)] #clean the NA value
})

#save(id_usedAnnots_list,file="Data/id_usedAnnots_list.RData")


#Create a dataframe format for NCHU

##find the strain that has the most annotations
sapply(id_usedAnnots_list,FUN = function(annots){
  length(annots)
}) %>% max #27


##extend the annotation vector so for each one length=27 ("" to make up for the blanks)
id_usedAnnots_list_flattened=sapply(id_usedAnnots_list,FUN=function(annots){ #sapply simplifies the result to a matrix
  if(length(annots)<27) annots=c(annots,rep("",27-length(annots)))
  return(annots)
})


id_usedAnnots_df=t(id_usedAnnots_list_flattened) 

dim(id_usedAnnots_df) #3979   27

phenotypes_annotations=cbind(All_Data_NAimputed,id_usedAnnots_df)
class(phenotypes_annotations)
dim(phenotypes_annotations) #3979 351 (324+27)



#save(phenotypes_annotations,file="Data/sourced/phenotypes_annotations.RData")

#write.table(phenotypes_annotations,file = "Data/phenotypes_annotations.txt",sep ="\t" ,quote = F)








