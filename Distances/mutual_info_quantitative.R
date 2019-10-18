#pairwise MI for quantitative data

#Test on the whole dataset
##start.time = Sys.time()
##mi_all=cminjk(t(All_Data_NAimputed))
##end.time = Sys.time()
##end.time - start.time  #Time difference of 27.39689 mins
##save(mi_all,file="Data/mi_all.RData")
load("Data/mi_all.RData")
mi_all=mi_all/log(2) #convert nats to bits
class(mi_all)
dim(mi_all)


mi_all_dist=mi_all %>% as.dist
mi_df=meltANDsort_dist(mi_all_dist,decreasing = T)
names(mi_df)=c("strain1","strain2","mi") #correct the name of the columns
str(mi_df)
hist(mi_df$mi)

summary(mi_df$mi)
summary(1-strain1strain2_allAnnotations_allDistances$mi) #this is the mi for ternary data

strain1strain2coPwyAnnot=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy")]
strain1strain2coPwyAnnot$strain1=as.character(strain1strain2coPwyAnnot$strain1)
strain1strain2coPwyAnnot$strain2=as.character(strain1strain2coPwyAnnot$strain2)
  
mi_df_coPwyAnnot=left_join(mi_df,strain1strain2coPwyAnnot,by= c("strain1","strain2"))
mi_df_coPwyAnnot$cumsum=cumsum(mi_df_coPwyAnnot$Pwy)
names(mi_df_coPwyAnnot)[4]="sameORnot"
str(mi_df_coPwyAnnot)



#pcc based cumsum
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcc")]
Pwy_pcc=cumsum(df$Pwy[order(df$pcc)])

#spearman based cumsum
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","spearman")]
Pwy_spearman=cumsum(df$Pwy[order(df$spearman)])

#plot
samples=seq(from=0,to=7914231,by=173); samples[1]=1; 
limit=3000 #limit=7914231
x=mi_df_coPwyAnnot

plot(samples[samples<=limit],x$cumsum[samples][samples<=limit],type="l")
lines(samples[samples<=limit],Pwy_pcc[samples][samples<=limit],col="cyan")
lines(samples[samples<=limit],Pwy_spearman[samples][samples<=limit],col="pink")
abline(a=0,b=1,col="blue") #positive control
abline(a=0,b=sum(x$sameORnot)/length(x$sameORnot),col="red")#negative control

#Conclusion: Good news! MI works the best for quantitative data
