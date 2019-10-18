#Test tSNE algorithm using Nichols' data
#Ref: https://www.analyticsvidhya.com/blog/2017/01/t-sne-implementation-r-python/

#install.packages("Rtsne")
library("Rtsne")

train<- All_Data_NAimputed

## Executing the algorithm on curated data
exeTimeTsne<- system.time(tsne <- Rtsne(train, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)) #Fitting performed in 23.40 seconds.

exeTimeTsne

## Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels="o") #Doesn't show any separable patterns. Why?


