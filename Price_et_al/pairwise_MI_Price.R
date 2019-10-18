load("Data/ternary_Nich_Price.RData")

##Ref: https://cran.r-project.org/web/packages/infotheo/infotheo.pdf
start.time = Sys.time()
ternary_dat=ternary_Nich_Price[,325:486] #only use the price data here
price_mi_ternary=( 1- mutinformation(t(ternary_dat) %>% as.data.frame) %>% natstobits ) %>% melt_dist 
end.time = Sys.time()
end.time-start.time #Time difference of 22.15182 mins on my PC




names(price_mi_ternary)=c("strain1","strain2","mi_ternary")
price_mi_ternary$strain1=as.integer(price_mi_ternary$strain1)
price_mi_ternary$strain2=as.integer(price_mi_ternary$strain2)
price_mi_ternary=left_join(price_mi_ternary,strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy")],by = c("strain1", "strain2"))
##The above have merging problem (produce lots of NA for pwy) because of the order of strain1, strain2

##swap some strain1 and strain2
price_mi_ternary[is.na(price_mi_ternary$Pwy),c("strain1","strain2")]=price_mi_ternary[is.na(price_mi_ternary$Pwy),c("strain2","strain1")]
##remerge
price_mi_ternary=left_join(price_mi_ternary[,c("strain1","strain2","mi_ternary")],strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy")],by = c("strain1", "strain2"))

anyNA(price_mi_ternary$Pwy) #all NAs are gone



dat_to_join=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2")]
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])==5 )

dat_to_join$sameORnot=as.numeric(TF)

price_mi_ternary_allAnnot=left_join(price_mi_ternary[,c("strain1","strain2","mi_ternary")],dat_to_join,by = c("strain1", "strain2"))



price_mi_ternary_allAnnot=price_mi_ternary_allAnnot[order(price_mi_ternary_allAnnot$mi_ternary),]
price_mi_ternary_allAnnot$cumsum=cumsum(price_mi_ternary_allAnnot$sameORnot)


##save(price_mi_ternary_allAnnot,file="Data/price_mi_ternary_allAnnot.RData")
