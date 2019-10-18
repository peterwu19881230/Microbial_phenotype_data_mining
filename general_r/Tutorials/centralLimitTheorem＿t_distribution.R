#Practice of Central Limit theorem
#Ref: https://en.wikipedia.org/wiki/Central_limit_theorem



#Distribution of mean of X distribution (no matter X is normal or not) ~ N( mu, sigma(X)/sqrt(n) ), where n = sample size (p.s. sample size != No. of sample)
#Usually sigma is unknown, so we would use (x_bar - mu_h0)/(sd/sqrt(n)) ~ t(df=n-1)

#Here I assume X is the population
X=rnorm(100000) # ~ Z
hist(X)

mean=numeric(1000)
for(i in 1:1000){
  mean[i]=sample(X,100) %>% mean
}

hist(mean,freq=F)
curve(dnorm(x,mean=0,sd=1/sqrt(100)),col="blue",lwd=2,add= T)


#Here I assume X2 is the population
X2=runif(100000) #uniform distribution. Has mean=0.5
hist(X2)

mean=numeric(1000)
sd=numeric(1000)
for(i in 1:1000){
  values=sample(X2,100)
  mean[i]=values %>% mean
  sd[i]=values %>% sd
}

hist((mean-0.5)/(sd/sqrt(100)),freq=F) # 0.5 is the mean for H0
curve(dt(x,df=100-1),col="blue",lwd=2,add= T)


#Here I assume X3 is the population
X3=rpois(100000,lambda=3) #poisson distribution. Has mean=3
hist(X3)

mean=numeric(1000)
sd=numeric(1000)
for(i in 1:1000){
  values=sample(X3,100)
  mean[i]=values %>% mean
  sd[i]=values %>% sd
}

hist((mean-3)/(sd/sqrt(100)),freq=F) # 3 is the mean for H0
curve(dt(x,df=100-1),col="blue",lwd=2,add= T)



#t-distribution
set.seed(102)
mean=numeric(1000)
sd=numeric(1000)
for(i in 1:1000){
  values=rnorm(100) #sample from Z
  mean[i]=values %>% mean
  sd[i]=values %>% sd
}

mean_h0=0

hist((mean-mean_h0)/(sd/sqrt(100)),freq=F) # 0 is the mean for H0
curve(dt(x,df=100-1),col="blue",lwd=2,add= T)


set.seed(102)
mean=numeric(1000)
sd=numeric(1000)
for(i in 1:1000){
  values=runif(100) #sample from uniform distribution
  mean[i]=values %>% mean
  sd[i]=values %>% sd
}

mean_h0=0.5

hist((mean-mean_h0)/(sd/sqrt(100)),freq=F) # 0.5 is the mean for H0
curve(dt(x,df=100-1),col="blue",lwd=2,add= T)


set.seed(102)
mean=numeric(100000)
sd=numeric(100000)
for(i in 1:100000){
  values=runif(3) 
  mean[i]=values %>% mean
  sd[i]=values %>% sd
}

mean_h0=0.5

hist((mean-mean_h0)/(sd/sqrt(100)),freq=F,breaks=seq(from=-3300,to=2000,by=0.01),xlim=c(-10,10)) # 3 is the mean for H0
curve(dt(x,df=3-1),col="blue",lwd=2,add= T)

##t-distribution with different df compared with ~N
curve(dt(x,df=2),col=1,xlim=c(-5,5),n=500,ylim=c(0,0.4)) #n= no. of points plotted at x axis (default=101 <=resolution is not good enough)
curve(dt(x,df=6),col=2,xlim=c(-5,5),add=T,n=500)
curve(dt(x,df=10),col=3,xlim=c(-5,5),add=T,n=500)

curve(dt(x,df=Inf),col=4,xlim=c(-5,5),add=T,n=500) #When df=Inf it is ~Z
curve(dnorm(x),col=5,xlim=c(-5,5),add=T,n=500) # ~Z





