names(id_allAttributes)

tab=id_allAttributes[,c("ids","Pwy","pcomplex","regulator","operon","kegg_modules")] %>% unique
dim(tab)



TF=apply(tab,1,FUN=function(row){
  (is.na(row) %>% sum)==5 #strains with no annotations from any of the sources have 5 NAs
})

strains=tab[TF,]$ids %>% unique #28 strains that are not anywhere annotated

## "657"  "159"  "3131" "3952" "658"  "168"  "367"  "2833" "2863" "3275" "3286" "3295"
## "3304" "3308" "3314" "3329" "3332" "3334" "3342" "3344" "3345" "3353" "3356" "3359"
## "3363" "3370" "3934" "3371"

