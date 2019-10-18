#Analysis of randomPCC_3, random_absPCC_3 made by cor_3_on_ada.R


unique_cor_strains=as.numeric(as.dist(cor_strains))
abs_unique_cor_strains=abs(unique_cor_strains)

#Note: InrandomPCC and random_absPCC are matrix where 1st row is the sample mean, 2nd row is the sample standard deviation

dim(randomPCC_3)
hist(main="Histogram of 3 gene mock pathway",randomPCC_3[1,],freq=F,breaks=200)
curve(dnorm(x,mean=mean(randomPCC_3[1,]),sd=sd(randomPCC_3[1,])),col="blue",add=T) #It seems to follow a ~N using

mean_h0=0
t=(randomPCC_3[1,]-mean_h0)/(randomPCC_3[2,]/sqrt(3))
hist(t,freq=F,breaks=seq(from=-500,to=600,by=0.1),xlim=c(-10,10)) #xlim is to limit the region the whole distribution is being plotted. Otherwise, the range will be too broad and the plot would be difficult to visualize
curve(dt(x,df=2),col="blue",add=T)



dim(random_absPCC_3)
mean_h0=mean(random_absPCC_3[1,])
t=(random_absPCC_3[1,]-mean_h0)/(random_absPCC_3[2,]/sqrt(3))
hist(t,freq=F,breaks=seq(from=-1500,to=400,by=0.1),xlim=c(-10,10)) #xlim is to limit the region the whole distribution is being plotted. Otherwise, the range will be too broad and the plot would be difficult to visualize
curve(dt(x,df=2),col="blue",add=T) #This isn't a t distribution at all!




#when n=48
hist(randomPCC[[26]],freq=F)
curve(dnorm(x,mean=mean_h0,sd=sd(unique_cor_strains)/sqrt(48)),col="blue",add=T) #Here I assume sd(unique_cor_strains) can be used to replace sigma (which is unknown)



hist(randomPCC[[1]],freq=F)
curve(dnorm(x,mean=mean_h0,sd=sd(unique_cor_strains)/sqrt(3)),col="blue",add=T) #Here I assume sd(unique_cor_strains) can be used to replace sigma (which is unknown)

