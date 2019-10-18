#1. Find the indices of strains in Nichols' that were duplicated
#2. Do any of them have any significant phenotypes?
#3. Are they annotated in any of the following databases?: 


#1
x=id_allAttributes[,c("ids","Original Name")] %>% unique
dim(x)
x$"Original Name" %>% duplicated %>% sum #12
id=x$ids[duplicated(x$"Original Name",fromLast=T)|duplicated(x$"Original Name",fromLast=F)]



#2
##Unfortunately most of them have significant phenotypes => Shouldn't be removed from the analysis
dat=Ternary_Data_324cutff_NAremoved[rownames(Ternary_Data_324cutff_NAremoved) %in% id,]
apply(dat,1,FUN=function(strain){
  sig=ifelse(strain==1|strain==-1,1,0)
  sum(sig)
})



#3
##Unfortunately they are well annotated and can't be removed
id_allAttributes[id,c("Pwy","pcomplex","regulator","operon","kegg_modules")]

