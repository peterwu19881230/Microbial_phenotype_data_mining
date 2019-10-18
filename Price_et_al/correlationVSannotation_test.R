#correlation VS annotation test

#Condition mapping from Curtis is here:
##https://docs.google.com/spreadsheets/d/1ivzQBRTA_02KkYyzyGx-Yjtak45FAmYIvWRmNfLyjWY/edit#gid=0


##Name mapping for Deutschbauers' to Nichols'
#===============================================================================================
#Place to download the data: http://genomics.lbl.gov/supplemental/bigfit/
new_dat=read.table("Data/fit_logratios_good.tab.txt",quote="",fill=T,sep="\t",header = T)
new_dat_scoresOnly=new_dat[,-(1:4)]


new_dat$sysName=as.character(new_dat$sysName)
#names(new_dat)
names(new_dat)[2]="bNumber"
duplicated(new_dat$bNumber) %>% sum #No duplicates in Deutschbauer's data (Awesome!)

bNumber_inDeutschbauer=cbind(new_dat$bNumber,T) %>% as.data.frame(stringsAsFactors=F)
names(bNumber_inDeutschbauer)=c("bNumber","in_Deutschbauer_s")

temp=left_join(id_allAttributes,bNumber_inDeutschbauer,by = "bNumber")
temp$in_Deutschbauer_s[is.na(temp$in_Deutschbauer_s)]=F

temp2=unique(temp[,c("ids","in_Deutschbauer_s")])

#sum(temp2$in_Deutschbauer_s %>% as.logical) #3525 strains mapped to Nichols
#===============================================================================================




##The id for Nichols' / bNumber for Deutschbauer's
strain_to_merge=temp[as.logical(temp$in_Deutschbauer_s)==T,c("ids","bNumber")] %>% unique
##Nichols' has duplicates => some of the Deutschbauer's will be duplicated

##get the indices for subtracting Deutschbauer's data
index=c()
count=0
for(bNumber in strain_to_merge$bNumber){
  count=count+1
  index[count]=which(new_dat$bNumber==bNumber)
}

Nich_Deut=cbind(All_Data_NAimputed[strain_to_merge$ids %>% as.numeric,],
                new_dat_scoresOnly[index,])

Deut=new_dat_scoresOnly[index,] #This only contains strains that are also in Nichols'
rownames(Deut)=rownames(Nich_Deut)
#str(Deut)

Nich_Price_quantitative=Nich_Deut

#save(Nich_Price_quantitative,file="Data/Nich_Price_quantitative.RData")


##Do pairwise PCCs look the same as those in Nichols (for the mapped strains)?
#===============================================================================================

Nich=All_Data_NAimputed[strain_to_merge$ids %>% as.numeric,]
Deut=new_dat_scoresOnly[index,]; rownames(Deut)=rownames(Nich)

identical(rownames(Deut),rownames(Nich)) #TRUE

Nich_pairwisePCC_Table=cor(t(Nich)) %>% as.dist %>% meltANDsort_dist #the result is a dataframe
names(Nich_pairwisePCC_Table)[3]="Distance_Nich"

Deut_pairwisePCC_Table=cor(t(Deut)) %>% as.dist %>% meltANDsort_dist #the result is a dataframe
names(Deut_pairwisePCC_Table)[3]="Distance_Deut"

pcc_comparison_table=left_join(Nich_pairwisePCC_Table,Deut_pairwisePCC_Table)

dim(pcc_comparison_table)

cor(pcc_comparison_table$Distance_Nich,pcc_comparison_table$Distance_Deut) #0.02946715 #Very bad correlation
cor(pcc_comparison_table$Distance_Nich,pcc_comparison_table$Distance_Deut,method="spearman") #0.009994566 #Very bad correlation

#The next thing would be to make them into categorical variables and then take various correlation (distance metrics)

#===============================================================================================



#Correlation VS Annotation experiment
#===============================================================================================


##original_attribute_list=attr_list(ECK.Pathway$ids,ECK.Pathway$Data)
##attribute_list=original_attribute_list[strain_to_merge$ids]

##dat=Nich_Deut
##start.time = Sys.time()
##Nich_Deut_AorrAnnotPwy_pcc=dist_TF_cumsum_matirxOperation(dat,attribute_list,pcc_dist) #pcc_dist is defined in functions.R
##end.time = Sys.time()
##end.time - start.time #Time difference of 35.36343 secs
#save(Nich_Deut_AorrAnnotPwy_pcc,file="Data/Nich_Deut_AorrAnnotPwy_pcc.RData")
load("Data/Nich_Deut_AorrAnnotPwy_pcc.RData")


##dat=Deut
##start.time = Sys.time()
##Deut_AorrAnnotPwy_pcc=dist_TF_cumsum_matirxOperation(dat,attribute_list,pcc_dist) #pcc_dist is defined in functions.R
##end.time = Sys.time()
##end.time - start.time #Time difference of 23.28664 secs
##save(Deut_AorrAnnotPwy_pcc,file="Data/Deut_AorrAnnotPwy_pcc.RData")
load("Data/Deut_AorrAnnotPwy_pcc.RData")


#control (If I only use data from Nichols)
##dat=All_Data_NAimputed[strain_to_merge$ids %>% as.numeric,]
##start.time = Sys.time()
##Nich_AorrAnnotPwy_pcc=dist_TF_cumsum_matirxOperation(dat,attribute_list,pcc_dist) #pcc_dist is defined in functions.R
##end.time = Sys.time()
##end.time - start.time #Time difference of 28.16337 secs
##save(Nich_AorrAnnotPwy_pcc,file="Data/Nich_AorrAnnotPwy_pcc.RData")
load("Data/Nich_AorrAnnotPwy_pcc.RData")

#plot
##7914231 -> 6211050

samples=seq(from=0,to=6211050,by=47); samples[1]=1
limit=3000 #If limit=6211050, it's clear Deut alone is much better
x=Nich_Deut_AorrAnnotPwy_pcc
x2=Deut_AorrAnnotPwy_pcc
x3=Nich_AorrAnnotPwy_pcc
plot(samples[samples<=limit],x$cumsum[samples][samples<=limit],type="l")
lines(samples[samples<=limit],x2$cumsum[samples][samples<=limit],col="cyan")
lines(samples[samples<=limit],x3$cumsum[samples][samples<=limit],col="green")
abline(a=0,b=1,col="blue") #positive control
abline(a=0,b=sum(x$sameORnot)/length(x$sameORnot),col="red")#negative control
#===============================================================================================



#Correlation VS Annotation experiment (Ternary)    
#===============================================================================================
#(!) I never figured out how Price et al. decided their significant cutoffs (again not clear). I filtered every fitness here based on |fitness|>=3
#ternary_dat=TernaryConvert(as.matrix(dat),thresh=3)
#(updated 7/22/2019) The authors replied. I have generated the ternary files based on: Price_ternary.R


load("Data/sourced/ternary_Price.RData")


#Use the indices defined on the code above to subset the strains and merge the ternary data of Nichols' and Price's. The resulting data only have strains that are in both Nichols' and Price's
ternary_Nich_Price=cbind(Ternary_Data_324cutff_NAremoved[strain_to_merge$ids %>% as.numeric,],
                        ternary_Price[index,])

#save(ternary_Nich_Price,file="Data/ternary_Nich_Price.RData")


##Ref: https://cran.r-project.org/web/packages/infotheo/infotheo.pdf
start.time = Sys.time()
ternary_dat=ternary_Nich_Deut
NichDeut_mi_ternary=( 1- mutinformation(t(ternary_dat) %>% as.data.frame) %>% natstobits ) %>% melt_dist 
end.time = Sys.time()
end.time-start.time #Time difference of 1.109382 hours

names(NichDeut_mi_ternary)=c("strain1","strain2","ternary_mi")
###save(NichDeut_mi_ternary,file="Data/NichDeut_mi_ternary.RData")


NichDeut_mi_ternary$strain1=as.integer(NichDeut_mi_ternary$strain1)
NichDeut_mi_ternary$strain2=as.integer(NichDeut_mi_ternary$strain2)
NichDeut_mi_ternary=left_join(NichDeut_mi_ternary,strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy")],by = c("strain1", "strain2"))
##The above have merging problem (produce lots of NA for pwy) because of the order of strain1, strain2

##swap some strain1 and strain2
NichDeut_mi_ternary[is.na(NichDeut_mi_ternary$Pwy),c("strain1","strain2")]=NichDeut_mi_ternary[is.na(NichDeut_mi_ternary$Pwy),c("strain2","strain1")]
##remerge
NichDeut_mi_ternary=left_join(NichDeut_mi_ternary[,c("strain1","strain2","mi_ternary")],strain1strain2_allAnnotations_allDistances[,c("strain1","strain2","Pwy")],by = c("strain1", "strain2"))

anyNA(NichDeut_mi_ternary$Pwy) #all NAs are gone


NichDeut_mi_ternary=NichDeut_mi_ternary[order(NichDeut_mi_ternary$mi_ternary),]
NichDeut_mi_ternary$cumsum=cumsum(NichDeut_mi_ternary$Pwy)

NichDeut_mi_ternary_pwy=NichDeut_mi_ternary

###save(NichDeut_mi_ternary_pwy,file="Data/NichDeut_mi_ternary_pwy.RData")
load("Data/NichDeut_mi_ternary_pwy.RData")


#plot
##7914231 -> 6211050

samples=seq(from=0,to=6211050,by=47); samples[1]=1
limit=3000 #If limit=6211050, it's clear Deut alone is much better
x=Nich_Deut_AorrAnnotPwy_pcc
x2=Deut_AorrAnnotPwy_pcc
x3=Nich_AorrAnnotPwy_pcc
x4=NichDeut_mi_ternary_pwy
plot(samples[samples<=limit],x$cumsum[samples][samples<=limit],type="l")
lines(samples[samples<=limit],x2$cumsum[samples][samples<=limit],col="cyan")
lines(samples[samples<=limit],x3$cumsum[samples][samples<=limit],col="green")
lines(samples[samples<=limit],x4$cumsum[samples][samples<=limit],col="brown")
abline(a=0,b=1,col="blue") #positive control
abline(a=0,b=sum(x$sameORnot)/length(x$sameORnot),col="red")#negative control

##I need tp get Nich only ternary to compare. qualitative from 2 datasets is close but not better than Nich's quantitative
##Also it would be worth doing quantitative MI for 2 datasets
#===============================================================================================







#Correlation VS Annotation experiment (Ternary with all 5 annotation sets)    
#===============================================================================================
load("Data/NichDeut_mi_ternary_pwy.RData") #this is the data that use only pwy annotation, but I can still use the strain1-strain2-MI to do the following experiment



dat_to_join=strain1strain2_allAnnotations_allDistances[,c("strain1","strain2")]
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])==5 )

dat_to_join$sameORnot=as.numeric(TF)

NichDeut_mi_ternary_allAnnot=left_join(NichDeut_mi_ternary_pwy[,c("strain1","strain2","mi_ternary")],dat_to_join,by = c("strain1", "strain2"))



NichDeut_mi_ternary_allAnnot=NichDeut_mi_ternary_allAnnot[order(NichDeut_mi_ternary_allAnnot$mi_ternary),]
NichDeut_mi_ternary_allAnnot$cumsum=cumsum(NichDeut_mi_ternary_allAnnot$sameORnot)


##save(NichDeut_mi_ternary_allAnnot,file="Data/NichDeut_mi_ternary_allAnnot.RData")
#===============================================================================================

