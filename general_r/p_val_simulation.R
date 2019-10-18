# Do p-values from sampling out of a standard normal distribution follow a uniform distribution?


#Refs:
#https://stats.stackexchange.com/questions/10613/why-are-p-values-uniformly-distributed-under-the-null-hypothesis
#https://www.cyclismo.org/tutorial/R/pValues.html


#Simulation 
p_val=c()
run_times=10000

set.seed(101)
for(i in 1:run_times){
  x_bar=rnorm(100) %>% mean
  z=(x_bar-0)/(1/sqrt(100))
  p_val[i]=2*pnorm(-abs(z))
}

hist(p_val)