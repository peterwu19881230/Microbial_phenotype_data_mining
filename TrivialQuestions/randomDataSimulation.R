#Fake data to do simulation: similar stuff like what's done for Nichols'

set.seed(102)
fakeData=matrix(rnorm(3979*324),ncol=324)
fakeStrainCC=cor(t(fakeData))

unique_fakeStrainCC=fakeStrainCC %>% as.dist %>% as.numeric
unique_cor_strains=cor_strains %>% as.dist %>% as.numeric

#Empiral CDF according to my understanding of this: https://en.wikipedia.org/wiki/Empirical_distribution_function
samples=c(1,(1:39)*202929) # Have to sample. Otherwise it's too big and will crash
sort_unique_cor_strains=sort(unique_cor_strains)
probability=(1:length(sort_unique_cor_strains))/length(sort_unique_cor_strains)
plot(sort_unique_cor_strains[samples],probability[samples])
abline(h=0.95,col="cyan")
abline(h=0.99,col="blue")
curve(pnorm(x,mean=mean(sort_unique_cor_strains),sd=sd(sort_unique_cor_strains)),from=-1,to=1,lty="dotted",add=T)

## Residual plot
resid=(1:length(sort_unique_cor_strains))/length(sort_unique_cor_strains)-
      pnorm(sort_unique_cor_strains,mean=mean(sort_unique_cor_strains),sd=sd(sort_unique_cor_strains))
plot(1:length(sort_unique_cor_strains),resid)
abline(h=0)


#Some histograms

hist(unique_fakeStrainCC,breaks=100)
#qqnorm(unique_fakeStrainCC) <-Why does this take so long?

hist(unique_cor_strains,breaks=100)
#qqnorm(unique_cor_strains) <-Why does this take so long?


#Nichols' strainCC to compared with fake PCC
FakeData=data.frame(fakeFitness=fakeData %>% t %>% as.vector)
NicholsData=data.frame(fitness=All_Data %>% t %>% as.vector)
FakeData_pcc=data.frame(fakeCC=unique_fakeStrainCC)
NicholsData_pcc=data.frame(CC=unique_cor_strains)

ggplot()+  
  geom_density(data=FakeData,alpha=.2, fill="blue", aes(x=fakeFitness))+
  geom_density(data=NicholsData,alpha=.2, fill="#FF6666", aes(x=fitness))
  

ggplot()+  
  #geom_histogram(aes(y=..density..),binwidth=.5,colour="black", fill="white")+ # Histogram with density instead of count on y-axis
  geom_density(data=FakeData_pcc,alpha=.2, fill="blue", aes(x=fakeCC))+
  geom_density(data=NicholsData_pcc,alpha=.2, fill="#FF6666", aes(x=CC))
  



#The 48-gene mock pathway looks pretty normal
ggplot()+geom_density(data=data.frame(cc26=randomPCC[[26]]),alpha=.1, fill="#FF6666", aes(x=cc26))


#Distribution of all pcc + distribution of mock pathways 
ggplot()+ 
  geom_density(data=NicholsData_pcc,alpha=.2, fill="#FF6666", aes(x=CC))+
  geom_density(data=data.frame(cc1=randomPCC[[1]]),alpha=.01, aes(x=cc1))+
  geom_density(data=data.frame(cc2=randomPCC[[2]]),alpha=.01, aes(x=cc2))+
  geom_density(data=data.frame(cc3=randomPCC[[3]]),alpha=.01, aes(x=cc3))+
  geom_density(data=data.frame(cc4=randomPCC[[4]]),alpha=.01, aes(x=cc4))+
  geom_density(data=data.frame(cc5=randomPCC[[5]]),alpha=.01, aes(x=cc5))+
  geom_density(data=data.frame(cc6=randomPCC[[6]]),alpha=.01, aes(x=cc6))+
  geom_density(data=data.frame(cc7=randomPCC[[7]]),alpha=.01, aes(x=cc7))+
  geom_density(data=data.frame(cc8=randomPCC[[8]]),alpha=.01, aes(x=cc8))+
  geom_density(data=data.frame(cc9=randomPCC[[9]]),alpha=.01, aes(x=cc9))+
  geom_density(data=data.frame(cc10=randomPCC[[10]]),alpha=.01, aes(x=cc10))+
  geom_density(data=data.frame(cc11=randomPCC[[11]]),alpha=.01, aes(x=cc11))+
  geom_density(data=data.frame(cc12=randomPCC[[12]]),alpha=.01, aes(x=cc12))+
  geom_density(data=data.frame(cc13=randomPCC[[13]]),alpha=.01, aes(x=cc13))+
  geom_density(data=data.frame(cc14=randomPCC[[14]]),alpha=.01, aes(x=cc14))+
  geom_density(data=data.frame(cc15=randomPCC[[15]]),alpha=.01, aes(x=cc15))+
  geom_density(data=data.frame(cc16=randomPCC[[16]]),alpha=.01, aes(x=cc16))+
  geom_density(data=data.frame(cc17=randomPCC[[17]]),alpha=.01, aes(x=cc17))+
  geom_density(data=data.frame(cc18=randomPCC[[18]]),alpha=.01, aes(x=cc18))+
  geom_density(data=data.frame(cc19=randomPCC[[19]]),alpha=.01, aes(x=cc19))+
  geom_density(data=data.frame(cc20=randomPCC[[20]]),alpha=.01, aes(x=cc20))+
  geom_density(data=data.frame(cc21=randomPCC[[21]]),alpha=.01, aes(x=cc21))+
  geom_density(data=data.frame(cc22=randomPCC[[22]]),alpha=.01, aes(x=cc22))+
  geom_density(data=data.frame(cc23=randomPCC[[23]]),alpha=.01, aes(x=cc23))+
  geom_density(data=data.frame(cc24=randomPCC[[24]]),alpha=.01, aes(x=cc24))+
  geom_density(data=data.frame(cc25=randomPCC[[25]]),alpha=.01, aes(x=cc25))+
  geom_density(data=data.frame(cc26=randomPCC[[26]]),alpha=.01, aes(x=cc26))


#Distribution of sample means
n=c(3:22,27,28,31,41,47,48)

pdf("Data/randomPCC_i.pdf")

for(i in 1:length(n)){
par(mfrow=c(1,3))
hist(main=paste("Histogram of ",n[i]," gene\n mock pathway",sep=""),randomPCC[[i]],breaks=100,freq=F)
curve(dnorm(x,mean=mean(unique_cor_strains),sd=sd(unique_cor_strains)/sqrt(n[i])),col="blue",lwd=2,add= T) #lwd is the thickness of the curve
##Note: this curve is according to the Central Limit Theorem => Distribution of mean pcc of n gene pathway

##qqplot with ~N
qqnorm(randomPCC[[i]])

##Residual plot (Using empirical CDF)
resid=(1:length(randomPCC[[i]]))/length(randomPCC[[i]])-
  pnorm(sort(randomPCC[[i]]),mean=mean(randomPCC[[i]]),sd=sd(randomPCC[[i]]))
plot(main = "Residual plot\n (Using empirical CDF)",1:length(randomPCC[[i]]),resid,pch=20,cex=0.2) #pch specifies the type of points (circle, triangle, square)
abline(h=0)
}
dev.off()

#Confidence for each n gene mock pathways
i=1:length(n)

## 95%
qnorm(1-0.05/2,mean=mean(unique_cor_strains),sd=sd(unique_cor_strains)/sqrt(n[i]))

## 99%
qnorm(1-0.01/2,mean=mean(unique_cor_strains),sd=sd(unique_cor_strains)/sqrt(n[i]))




