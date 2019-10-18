#Debug bottomUpExp()

set.seed(101)
ids=sample(1:1000,30)

temp=hclust(pcc_dist(All_Data[ids,]))
temp$labels=as.numeric(temp$labels)
plot(temp)

#create a fake table for experiment
tags=c()
for(i in 1:10){
  tag=sample(1:30,1)
  print(tag)
  tags=c(tags,rep(tag,3))
}
tags=as.character(tags)
tab=data.frame(ids,tags,stringsAsFactors = F)

bottomUpExp(temp,tab) #This doesn't work 

elements.hclust(temp,10) #This works
fractionBottomUp(temp,c(1,2,3,4,5,6)) #This works

#Guess: labels are not ordered like 1,2,3,4.....n. Instead, they are like sample(1:k,n)
#Problem fixed by modifying functions.R

