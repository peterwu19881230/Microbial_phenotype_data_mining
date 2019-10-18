#How to plot distributions using basic graph and ggplot

#Normal distribution is used as an example:
set.seed(1)
dat=rnorm(1000)


##Basic graph
hist(dat) #histogram 

hist(dat,freq = F) #density plot
curve(dnorm(x),col="blue",add=T) 
##->append the curve. Any function (including self-defined) can be used as long as the first argument = x


##ggplot version
x=data.frame(dat=dat) #(I hate this so much) ggplot only deals with data frames as input


ggplot(data = x)+geom_histogram(aes(x=dat)) #historgram

ggplot(data = x)+geom_density(aes(x=dat)) #Density plot


ggplot()+ #histogram + density plot (have to use density for histogram)
geom_histogram(data = x,aes(x=dat,y=..density..))+ #have to use density to overlay with the curve below
geom_density(aes(x=dat),colour="blue")










