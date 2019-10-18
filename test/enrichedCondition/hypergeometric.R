#Draw a probability mass function of a hypergeometric distribution

N=3979*324
success=15833

values=c(rep(1,success),rep(0,N-success))
values=sample(values)

noOfDraw=3979

nOfSuccess=numeric(10^4)
for(i in 1:10^4){ #this takes about 10 sec
  nOfSuccess[i]=sample(values,noOfDraw) %>% sum
  
}

#modify by using the function
 

hist(nOfSuccess,breaks=100,freq=F)
curve(dhyper(x,m=success,n=N-success,k=noOfDraw),from=0,to=100,col="blue",add=T) #dhyper does something weird here (~by Dr. Karmakar)
#Why do the graphs look so different?
