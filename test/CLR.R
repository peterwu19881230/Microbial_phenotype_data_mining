#Context Likelihood or (in the paper it's "of") Relatedness Network
##https://rdrr.io/cran/parmigene/man/clr.html
##paper: (shown as ref in the above url)

#install.packages("parmigene")
library(parmigene)

#test example
mat <- matrix(rnorm(1000), nrow=10)
mi  <- knnmi.all(mat)
grn <- clr(mi)



#test on Nichols' data

#generate grn
##============================================
#mat <- All_Data_NAimputed

#start.time = Sys.time()
#mi  <- knnmi.all(mat)
#grn <- clr(mi)
#end.time = Sys.time()
#end.time - start.time  #Time difference of 16.21768 mins
##============================================
#save(grn,file="Data/grn.RData")
load("Data/grn.RData")

strain1_strain2_grn=grn %>% as.dist %>% meltANDsort_dist(decreasing=T)
names(strain1_strain2_grn)=c("strain1","strain2","grn")



strain1strain2coPwyAnnot=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy")]
strain1strain2coPwyAnnot$strain1=as.character(strain1strain2coPwyAnnot$strain1)
strain1strain2coPwyAnnot$strain2=as.character(strain1strain2coPwyAnnot$strain2)

grn_df_coPwyAnnot=left_join(strain1_strain2_grn,strain1strain2coPwyAnnot,by= c("strain1","strain2"))


grn_df_coPwyAnnot$cumsum=cumsum(grn_df_coPwyAnnot$Pwy)
names(grn_df_coPwyAnnot)[4]="sameORnot"
str(grn_df_coPwyAnnot)


#pcc based cumsum
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcc")]
Pwy_pcc=cumsum(df$Pwy[order(df$pcc)])

#spearman based cumsum
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","spearman")]
Pwy_spearman=cumsum(df$Pwy[order(df$spearman)])

#plot
samples=seq(from=0,to=7914231,by=173); samples[1]=1; 
limit=3000 #limit=7914231
x=grn_df_coPwyAnnot

plot(samples[samples<=limit],x$cumsum[samples][samples<=limit],type="l")
lines(samples[samples<=limit],Pwy_pcc[samples][samples<=limit],col="cyan")
lines(samples[samples<=limit],Pwy_spearman[samples][samples<=limit],col="pink")
abline(a=0,b=1,col="blue") #positive control
abline(a=0,b=sum(x$sameORnot)/length(x$sameORnot),col="red")#negative control

#Conclusion: CLR works similarly as spearman does







