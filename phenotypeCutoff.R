#Goal: determine which scores should be the cutoffs for phenotypes in Nichols' (Previously we thought were -3 and 3)
#Understanding FDR: https://eranraviv.com/understanding-false-discovery-rate/

#I did part I and part II according to Nichols' (Described in the supplementary session: Phenotype Analysis). However it doesn't seem to be right
#=============================
#Part I: negative phenotypes
#Part II: positive phenotypes
#=============================

#Part III: If FDR is defined using 2 sides:

#Part IV: Verification

#Part V: FDR using each of the 324 condition 

#part I
allFitness=as.numeric(as.matrix(All_Data))
allFitness=allFitness[!is.na(allFitness)]

##Try different p
numbers=seq(from=0.025,to=0.0001,by=-0.0001)
for(p in numbers){
FDR=p/(sum(allFitness<=qnorm(p))/length(allFitness))
  if(FDR<=0.05){
    print(p)
    print(qnorm(p)) #The fitness score where it has 5% FDR
  }
}

###Verification
0.0004/(sum(allFitness<=qnorm(0.0004))/length(allFitness)) #FDR=0.0439987



##Try different q (q=qnorm(p), p=pnorm(q))
numbers=seq(from=-2.5,to=-4,by=-0.001)
n=1
for(q in numbers){
  FDR=pnorm(q)/(sum(allFitness<=q)/length(allFitness))
  if(FDR<=0.05){
    n=n+1
    print(pnorm(q))
    print(q) #The fitness score where it has 5% FDR
    if(n==50)break
    }
}

###Verification
pnorm(-3.307)/(sum(allFitness<=-3.307)/length(allFitness)) #FDR=0.04997598



##Determine the cutoff using real scores
numbers=sort(allFitness[allFitness<=-3],decreasing = T) #Use scores<=-3 to speed up the process
n=1
for(q in numbers){
  FDR=pnorm(q)/(sum(allFitness<=q)/length(allFitness))
  if(FDR<=0.05){
    n=n+1
    print(pnorm(q))
    print(q) #The fitness score where it has 5% FDR
    if(n==50)break
  }
}

###Verification
pnorm(-3.306904)/(sum(allFitness<=-3.306904)/length(allFitness)) #FDR=0.04998899

###This way the No. of negative phenotypes would be:
length(allFitness[allFitness<=-3.306904]) #12143


#part II
numbers=sort(allFitness[allFitness>=3],decreasing = F) #Use scores>=3 to speed up the process
n=1
for(q in numbers){
  FDR=(1-pnorm(q))/(sum(allFitness>=q)/length(allFitness))
  if(FDR<=0.05){
    n=n+1
    print(1-pnorm(q))
    print(q) #The fitness score where it has 5% FDR
    if(n==50)break
  }
}

###Verification
(1-pnorm(3.824592))/(sum(allFitness>=3.824592)/length(allFitness)) #FDR=0.04993414

###This way the No. of positive phenotypes would be:
length(allFitness[allFitness>=3.824592]) #1688


### Calculate the claimed ratio of negative phenotypes: positive phenotypes (=4:1)
12143/1688 #7.19372 => Doesn't match the paper's ratio


#Part III: If FDR is defined using 2 sides:
Zcutoff=seq(from=3.45,to=3.5,by=(3.5-3.45)/1000 )#3.45 is just a value that I picked up.How do I solve the Zcutoff mathematically?

n=1
for(q in Zcutoff){
  scores=sort(allFitness[allFitness<=-q | allFitness>=q],decreasing = F) #Use scores>=3 to speed up the process
  #FDR
  FDR=2*(1-pnorm(q))/(length(scores)/length(allFitness))
  if(FDR<=0.05){
    n=n+1
    print(q) #The closest: 3.463
    print(FDR) #The closest correspoinding to q=3.463 =>0.04999856
    #if(n==50)break
  }
  
}

 
## Calculate the claimed ratio of negative phenotypes: positive phenotypes (~4:1)
length(allFitness[allFitness<=-3.463])/length(allFitness[allFitness>=3.463]) #3.683243: 1 => close to 4:1


##Conclusion: The cutoff should be -3.463<=fitness score<=3.463


#Part IV: Verification 

##The paper says 49% of the strains tested had at least 1 phenotype if using FDR=5%
new.binary=All_Data
new.binary=ifelse(new.binary<=-3.463 | new.binary>=3.463,1,0)

###No. of Phenotypes for each strain under fitness score cutoff=3.463
noOfPhenotypes=apply(new.binary,1,sum,na.rm=T)

sum(noOfPhenotypes>=1)/3979 # 0.5001257 => close to the reported 49% (1957/3979)

### total number of phenotypes
sum(new.binary,na.rm=T) #13750 => Close. reported was 13497

### What is "phenotypes per screen"? Verified with Dr. Siegele => per screen means per condition
noOfPhenotypes.cond=apply(new.binary,2,sum,na.rm=T)
max(noOfPhenotypes.cond) # 146 => The reported value is 173
hist(noOfPhenotypes.cond) # similar to their fig. S2

##The paper says 5% FDR or less was calculated for each phenotype

###histogram of fitness scores and the H0 curve (~Z)
hist(allFitness,breaks = 500,freq=F) 
curve(dnorm(x),xlim=c(-40,10),col="blue",add=T,n=10000)


#Calculate FDR for each fitness score (This would take estimated 1hr on my desktop. I haven't run) => Make sure this is the right way
##In the paper they say FDR for each fitness scores were calculated. That's why I did this seemingly redundant calculation to verify.

'
start.time=Sys.time()

total=length(allFitness)
FDR=numeric(total)

n=1
for(q in allFitness){ #Iterate through all fitness scores can find me all possible fractions

  if(q>=0) leng=sum(allFitness>=q) 
  if(q<0)  leng=sum(allFitness<=q)
  
  FDR[n]=(1-pnorm(abs(q)))/(leng/total)
  n=n+1
}

end.time=Sys.time()
end.time-start.time
save(FDR,file="Data/sourced/FDR.RData") 
'



#Part V: FDR using each of the 324 condition
start.time=Sys.time()

FDR_cond=list()

total=3979 
for(i in 1:324){
  
  n=1
  dat=All_Data_NAimputed[,i] #Some conditions have NAs, but they constitute only a small amount. I think ignoring them would be fine.
  fdr=numeric(total)
  for(q in dat){ #Iterate through all fitness scores can find me all possible fractions
    
    if(q>=0) leng=sum(dat>=q) 
    if(q<0)  leng=sum(dat<=q)
    
    fdr[n]=(1-pnorm(abs(q)))/(leng/total)
    n=n+1
  }
  
  FDR_cond[[i]]=cbind(score=dat,FDR=fdr) %>% as.data.frame
}

end.time=Sys.time()
end.time-start.time
#save(FDR_cond,file="Data/sourced/FDR_cond.RData") 

##summary statistics:

min_forEachCond=sapply(FDR_cond,FUN=function(q_fdr){
q_fdr$score[q_fdr$FDR<=0.05] %>% abs %>% min   
})

hist(min_forEachCond) # Have to double check this
#save(min_forEachCond,file="Data/sourced/min_forEachCond.RData") 



#Using FDRs calculated by each condition instead of by all scores:
#Ternary_Data_324cutff is used (generated in phnotypeData.R)


##I got more significant phenotypes out of this ( sum(Ternary_Data_NAnotimputed %in% c(1,-1),na.rm=T) => this gets me 13750 )
sum(Ternary_Data_324cutff %in% c(1,-1),na.rm=T) #15833


## Check other statistics:

### Calculate the claimed ratio of negative phenotypes: positive phenotypes (=4:1)
sum(Ternary_Data_324cutff %in% -1,na.rm=T)/sum(Ternary_Data_324cutff %in% 1,na.rm=T) # 3.617381


###The paper says 49% of the strains tested had at least 1 phenotype if using FDR=5%
noOfPhnotypes=sum(apply(Ternary_Data_324cutff,1,FUN=function(row){
  if(sum(row %in% c(-1,1))>=1){
    return(1)
  }else return(0)
}))

noOfPhnotypes/3979 #55.5%. Dr. Siegele and I both think it's ok

###phenotypes per screen (condition)
noOfPhenotypes.cond=apply(Ternary_Data_324cutff,2,FUN=function(col){
  sig=ifelse(col %in% c(-1,1),1,0)
  return(sum(sig))
})
max(noOfPhenotypes.cond) # 194 => The reported value is 173
hist(noOfPhenotypes.cond,breaks=20) # Not very different from fig. S2





