#The following is used to test and debug dist_TF_cumsum_old
#===============================================================================

##A toy size
dat=matrix(rnorm(30),nrow=3)
rownames(dat)=c("1","2","3")

a=c("1","1","2","2","2","3") #name
b=c("A","B","A","B","C","D") #attribute

attribute_list=list()
for(i in 1:length(rownames(dat))){
  TFvector=(a==rownames(dat)[i])
  attribute_list[[i]]=b[TFvector]
}


##Another toy size
dat=matrix(rnorm(30),nrow=5)
rownames(dat)=c("1","2","3","4","5")

a=c("1","1","2","2","2","3","4","5") #name
b=c("A","B","A","B","C","D",NA,NA) #attribute

attribute_list=list()
for(i in 1:length(rownames(dat))){
  TFvector=(a==rownames(dat)[i])
  attribute_list[[i]]=b[TFvector]
}


##Some real data
dat=All_Data_NAimputed[1:100,]

a=ECK.Pathway$ids[ECK.Pathway$ids %in% 1:100]
b=ECK.Pathway$Data[ECK.Pathway$ids %in% 1:100]

attribute_list=list()
for(i in 1:length(rownames(dat))){
  TFvector=(a==rownames(dat)[i])
  attribute_list[[i]]=b[TFvector]
}

start.time = Sys.time()
listOfResults=dist_TF_cumsum_old(dat,attribute_list)
end.time = Sys.time()
end.time - start.time


##Plot them
samples=seq(from=0,to=7914231,by=1989); samples[1]=1 #change samples to 1:10000 to get the first 10000 results
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[2]]$cumsum[samples]),aes(no,CumSum,color="PCC_squared"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[3]]$cumsum[samples]),aes(no,CumSum,color="euclidean"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[4]]$cumsum[samples]),aes(no,CumSum,color="maximum"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[5]]$cumsum[samples]),aes(no,CumSum,color="manhattan"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[6]]$cumsum[samples]),aes(no,CumSum,color="canberra"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[7]]$cumsum[samples]),aes(no,CumSum,color="binary"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[8]]$cumsum[samples]),aes(no,CumSum,color="minkowski"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[9]]$cumsum[samples]),aes(no,CumSum,color="spearman"))+  
  geom_abline(intercept=,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same pwy/ptcomplex",aesthetic='custom text')+
  scale_color_discrete(name="Distance")


# Test one_attr:

## Test
a=c(1,1,2,2,2,3) #name
b=c("A","B","A","B","C","D") #attribute
a_b=data.frame(a=a,b=b,stringsAsFactors = F)
one_attr(a,b)

##Test using data with NA
a=c(1,1,2,2,2,3,4,5) #name
b=c("A","B","A","B","C","D",NA,NA) #attribute
a_b=data.frame(a=a,b=b,stringsAsFactors = F)
one_attr(a,b)


## Test using real data
dat=ECK.Pathway
dat$Data[dat$Data=="Not in any pathway"]=NA

start=Sys.time()
vec=one_attr(dat$ids[1:1000],dat$Data[1:1000])
end=Sys.time()
end-start #Time difference of 10.56773 secs