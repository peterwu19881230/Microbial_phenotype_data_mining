#Mutual information example from Nichols'


#Dr.Cai's suggestion for MI:
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3929353/
##https://www.mathworks.com/matlabcentral/fileexchange/30998-kernel-estimate-for-conditional-mutual-information
##https://cran.r-project.org/web/packages/mpmi/mpmi.pdf 

#========================================================================================================================================================
X=All_Data_NAimputed[1,] %>% as.numeric
Y=All_Data_NAimputed[2,] %>% as.numeric

cor(X,Y) #0.3041736


#About kernal density: https://en.wikipedia.org/wiki/Kernel_density_estimation ( down below is says R implements kernal density by density() )

#histogram with density
hist(X,breaks = 50,freq=F) #ref: https://www.r-bloggers.com/histogram-density-plot-combo-in-r/
lines(density(X),col="blue") 

hist(Y,breaks = 50,freq=F) #ref: https://www.r-bloggers.com/histogram-density-plot-combo-in-r/
lines(density(Y),col="blue") 


#density only
plot(density(X)) #ref: https://www.statmethods.net/graphs/density.html
plot(density(Y))


f_X=approxfun(density(X)) #ref: https://stackoverflow.com/questions/43570344/calculate-probability-from-density-function
curve(f_X(x),xlim=c(-7,7)) #verify that the function is correct (should look almost identical to the kernal density plot of X)


#Mutual informtaion 

#This one doesn't work -> by viewing the definition of Entropy()  I think this only works for categorical data. It doesn't work for continuous data
#Ref: https://www.rdocumentation.org/packages/DescTools/versions/0.99.19/topics/Entropy
#library("DescTools")
#Entropy(X)+Entropy(Y)-Entropy(X,Y)
#MutInf(X,Y) #The default unit uses bit(base for log()=2)

H_X=-sum( f_X(X)*log(f_X(X)) ) #default base for log() is e (2.71828...)

#I should try to figure out how to get this function: p(x,y)

#This one works: https://cran.r-project.org/web/packages/mpmi/mpmi.pdf 

#install.packages("mpmi")
library(mpmi)

#Test on a pair
v1=All_Data_NAimputed[3061,] %>% as.numeric
v2=All_Data_NAimputed[3107,] %>% as.numeric

cor(v1,v2)
mi=cminjk.pw(v1,v2) 
##They didn't say anything about the unit in their doc (##https://cran.r-project.org/web/packages/mpmi/mpmi.pdf) => In the following code I have figured out
mi


v1_ternary=Ternary_Data_324cutff_NAremoved[3061,]
v2_ternary=Ternary_Data_324cutff_NAremoved[3107,]

mi_ternary=dminjk.pw(v1_ternary,v2_ternary)
mi_ternary #0.1377335



library("DescTools")
mi_ternary_bit=MutInf(v1_ternary,v2_ternary) #0.1987074

mi_ternary/mi_ternary_bit #This confirms that dminjk.pw() uses the unit "nat"
##=> I therefore assume that other functions for continuous data also use "nat"
#========================================================================================================================================================




#What about mutual information for binary data?
#========================================================================================================================================================

#Note: functions used below are defined in the function directory


#Test

##test cases:

x1=c(1,1,1,1,1)
x2=c(1,1,1,1,1)

x3=c(0,0,1,0,0)
x4=c(0,0,1,0,0)

x9=c(0,1,1,0,0)
x10=c(0,1,1,0,0)

x11=c(0,1,1,1,0)
x12=c(0,1,1,1,0)

x13=c(0,1,1,1,1)
x14=c(0,1,1,1,1)

x5=c(0,0,0,0,0)
x6=c(0,0,0,0,0)

x7=c(0,0,1,0,0)
x8=c(0,0,-1,0,0)

x15=c(1,1,1,1,1)
x16=c(0,0,0,0,0)

x17=c(1,1,1,1,1)
x18=c(-1,-1,-1,-1,-1)


##test
library("DescTools") #my implementation matches MutInf() in this library
mutualInfo_binary(x1,x2) #0 -> In our case we want it to have perfect score (=1)
MutInf(x1,x2) #0

mutualInfo_binary(x3,x4) #0.7219281
MutInf(x3,x4) #0.7219281

mutualInfo_binary(x5,x6) #0
MutInf(x5,x6) #0

mutualInfo_binary(x7,x8) #0.7219281
MutInf(x7,x8) #0.7219281


Ternary_Data_324cutff_condCollapsed[1:2,1:20]
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[1,1:20],Ternary_Data_324cutff_condCollapsed[2,1:20]) #0.09040651
MutInf(Ternary_Data_324cutff_condCollapsed[1,1:20],Ternary_Data_324cutff_condCollapsed[2,1:20])

#For a strain pair that has PCC ~ 0.96 (highest PCC)
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[3061,],Ternary_Data_324cutff_condCollapsed[3107,]) #0.428811
#PCC ~ 0.8
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[2515,],Ternary_Data_324cutff_condCollapsed[3002,]) #0.428811
#PCC ~ 0.7
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[1211,],Ternary_Data_324cutff_condCollapsed[1495,]) #0.3603876
#PCC ~ 0.65
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[2336,],Ternary_Data_324cutff_condCollapsed[2431,]) #0.1032528
#PCC ~ 0.6
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[1699,],Ternary_Data_324cutff_condCollapsed[1742,]) #0
#PCC ~ 0.5
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[1425,],Ternary_Data_324cutff_condCollapsed[1637,]) #0
#PCC ~ 0.4
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[758,],Ternary_Data_324cutff_condCollapsed[1746,]) #0
#PCC ~ 0.3
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[908,],Ternary_Data_324cutff_condCollapsed[1335,]) #0
#PCC ~ 0 (lowest PCC)
mutualInfo_binary(Ternary_Data_324cutff_condCollapsed[398,],Ternary_Data_324cutff_condCollapsed[2571,]) #0



#If No. of conditions are fixed and I change the ratio of 1s, what will the trend be?

mutualInfos=numeric(324)
for(i in 1:324){
  identical_seq=c(rep(1,i),rep(0,324-i))
  mutualInfos[i]=mutualInfo_binary(identical_seq,identical_seq)
}

plot(1:324,mutualInfos)

#maximum:
mutualInfo_binary(c(rep(1,162),rep(0,162)),c(rep(1,162),rep(0,162)))




#PCC ~ 0.96 (highest PCC)
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[3061,],Ternary_Data_324cutff_condCollapsed[3107,]) #0.3079791
#PCC ~ 0.8
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[2515,],Ternary_Data_324cutff_condCollapsed[3002,]) #0.3079791
#PCC ~ 0.7
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[1211,],Ternary_Data_324cutff_condCollapsed[1495,]) #0.2849059
#PCC ~ 0.65
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[2336,],Ternary_Data_324cutff_condCollapsed[2431,]) #0.0920689
#PCC ~ 0.6
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[1699,],Ternary_Data_324cutff_condCollapsed[1742,]) #0
#PCC ~ 0.5
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[1425,],Ternary_Data_324cutff_condCollapsed[1637,]) #0
#PCC ~ 0.4
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[758,],Ternary_Data_324cutff_condCollapsed[1746,]) #0
#PCC ~ 0.3
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[908,],Ternary_Data_324cutff_condCollapsed[1335,]) #0
#PCC ~ 0 (lowest PCC)
mutualInfo_binary_weightedForNichols(Ternary_Data_324cutff_condCollapsed[398,],Ternary_Data_324cutff_condCollapsed[2571,]) #0


##Some notes about mutual info:
###When there are lots of (1,1) or (-1,-1) for 2 strains, the MI will also be low. However, in Nichols' Ternary data is is a sparse matrix => This problem doesn't exist
###See the ratio of 1 or -1 for each strain:
strainSigPhenotype=apply(Ternary_Data_324cutff,1,FUN=function(strain){
  ifelse(strain %in% c(1,-1),1,0) %>% sum(na.rm=T)
})
summary(strainSigPhenotype) #max no. of phenotypes for a strain = 119, which is 119/324 (36.7%) of total conditions => if not >50, the above mentioned problem shouldn't bother
hist(main="Histogram + Kernal Density of no. of significant phenotypes for each strain",xlab="No. of significant phenotypes",strainSigPhenotype,freq=F,breaks=70,ylim=c(0,0.5))
lines(density(strainSigPhenotype),col="blue",add=T)
#========================================================================================================================================================








