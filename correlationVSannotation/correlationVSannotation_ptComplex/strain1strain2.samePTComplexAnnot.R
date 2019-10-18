#Goal: Create a molten table: strain1 - strain2 - at least 1 same ptcomplex annotation or not (obj name: strain1strain2.sameAnnot)

# Prepare the annotation table 
##just use this defined obj from proteinComplex.R: ids.originalName.pcomplex

# Use the defined obj: sort.pcc.sql.NoIdent
id_in_complex=unique(ids.originalName.pcomplex$id) #total 618 unique ids, 618/3979 ~ 15% genes are involved in protein complexes

start.time=Sys.time()
bi.ptcomplex.annot=apply(sort.pcc.sql.NoIdent,1,FUN=function(row){
  strain1=row[1]
  strain2=row[2]
  if((strain1 %in% ids.originalName.pcomplex$id) & (strain2 %in% ids.originalName.pcomplex$id)){ #if both of them are annotated with protein complexes
    howmany=ids.originalName.pcomplex$pcomplex[ids.originalName.pcomplex$id==strain1] %in% ids.originalName.pcomplex$pcomplex[ids.originalName.pcomplex$id==strain2]
    if(sum(howmany)>=1){
    return(1)  
    }else{0}
    
  }else{0}
}
)
end.time=Sys.time()
end.time-start.time #Time difference of 3.284059 mins
#Why did strain1 %in% strain1 return 0.02..?


strain1strain2.pcc.samePTcomplex_Annot=cbind(sort.pcc.sql.NoIdent,bi.ptcomplex.annot)
save(strain1strain2.pcc.samePTcomplex_Annot,file="Data/strain1strain2.pcc.samePTcomplex_Annot.RData")


##Table that contains only genes involved in protein complexes
TFvec=(strain1strain2.pcc.samePTcomplex_Annot$strain1 %in% id_in_complex) & (strain1strain2.pcc.samePTcomplex_Annot$strain2 %in% id_in_complex)
inptcomplex_strain1strain2.pcc.samePTcomplex_Annot=strain1strain2.pcc.samePTcomplex_Annot[TFvec,]
save(inptcomplex_strain1strain2.pcc.samePTcomplex_Annot,file="Data/inptcomplex_strain1strain2.pcc.samePTcomplex_Annot.RData")



