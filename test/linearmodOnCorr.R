#regress pcc, spearman, mi... on co-annotation (1 or 0)
names(strain1strain2_allAnnotations_allDistances)

#try pathway first
result=lm( Pwy ~ pcc+spearman+mi,data=strain1strain2_allAnnotations_allDistances) #note these are distances
summary(result)

#next will be to use probability=beta1*pearson+beta2*spearman+beta3*mi 

##Verification with Nichols (should be very good because fitting is based on this)
pcc=strain1strain2_allAnnotations_allDistances$pcc
spearman=strain1strain2_allAnnotations_allDistances$spearman
mi=strain1strain2_allAnnotations_allDistances$mi
prob=-0.0313067*pcc+0.0281014*spearman+0.0519444*mi



#plot (compare pearson with newDistanceY)
samples=seq(from=0,to=7914231,by=1989); samples[1]=1
limit=5000

annot=strain1strain2_allAnnotations_allDistances$Pwy

cumsum_prob=annot[order(prob,decreasing = T)] %>% cumsum 

cumsum1=annot[order(strain1strain2_allAnnotations_allDistances$pcc)] %>% cumsum
cumsum2=annot[order(strain1strain2_allAnnotations_allDistances$spearman)] %>% cumsum
cumsum3=annot[order(strain1strain2_allAnnotations_allDistances$mi)] %>% cumsum

plot(samples[samples<=limit],cumsum_prob[samples][samples<=limit],type="l")
lines(samples[samples<=limit],cumsum1[samples][samples<=limit],col="cyan")
lines(samples[samples<=limit],cumsum2[samples][samples<=limit],col="pink")
lines(samples[samples<=limit],cumsum3[samples][samples<=limit],col="green")
abline(a=0,b=1,col="blue") #positive control
abline(a=0,b=max(cumsum1)/length(cumsum1),col="red")#negative control


##Verification with Price's (If it's good. Woohoo.)
