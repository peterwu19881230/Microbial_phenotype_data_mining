#UpSet graphing using Nichols'
##Refs:
##https://caleydo.org/tools/upset/
##https://github.com/hms-dbmi/UpSetR


#install.packages("hms-dbmi/UpSetR")
library(UpSetR)

#Test their example
#==================================================================================================
##https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html

# example of list input (list of named vectors)
listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), 
                  two = c(1, 2, 4, 5, 10), 
                  three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))

upset(fromList(listInput), order.by = "freq") #Note: intersection of one & three means: one & three but not two
#==================================================================================================


#Using Nichols' data
#==================================================================================================

## Data obj preparation
#------------------------------------------------------------------------------------------------------------
names(id_allAttributes)

#Here I haven't accounted strains that have 2 or more annotations from any annotation set

getAnnotatedID=function(attr){
  id_attr=id_allAttributes[,c("ids",attr)] %>% unique
  annotatedID=id_attr$ids[!is.na(id_attr[[paste(attr)]])] %>% unique
}

allID=as.character(1:3979)
Pwy=getAnnotatedID("Pwy")
kegg_modules=getAnnotatedID("kegg_modules")
pcomplex=getAnnotatedID("pcomplex")
regulator=getAnnotatedID("regulator")
operon=getAnnotatedID("operon")

sets=list(allID=allID,Pwy=Pwy,kegg_modules=kegg_modules,pcomplex=pcomplex,regulator=regulator,operon=operon)

names(sets)=c("All 3979 strains","EcoCyc Pathway","KEGG Modules","EcoCyc Protein complex","Regulon","Operon")
str(sets)
#------------------------------------------------------------------------------------------------------------


## Plot the result

upset(fromList(sets),nset=6,text.scale=2.5) #When nset is not specified, it only shows the 5 sets with larger number of ids (Weird default!)
#upset(fromList(sets),nset=6,text.scale=2.5,empty.intersections = "on")


#(!) the bottom left number "4000" will always get cut. Haven't found a solution
dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
pdf(file=paste(dir_of_workingScript,"UpSet_allAnnotationSets",sep="/"),width=7*2,height=7, onefile=FALSE)
##Note: if onefile=FALSE is not there, there will be an additional empty page (I don't know why)
##Solution Ref: https://stackoverflow.com/questions/12481267/in-r-how-to-prevent-blank-page-in-pdf-when-using-gridbase-to-embed-subplot-insi

upset(fromList(sets),nset=6,text.scale=2.5)
dev.off()

#==================================================================================================
