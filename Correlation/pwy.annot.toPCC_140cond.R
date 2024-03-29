#Goal: Bind pathway annotations to each pair of genes under a pcc


# Prepare the annotation table (change the "Not in any pathway" annotation to NA. Otherwise the genes that are not in any pathway are going to be grouped as if they are in the same pathway)
new.ECK.Pathway=ECK.Pathway
new.ECK.Pathway$Data=sapply(new.ECK.Pathway$Data,FUN=function(Data){
  ifelse(Data=="Not in any pathway",NA,Data)
})

# Get table for pcc and convert all pcc using abs()
'
abs.sort.nonRedun.140cond.strainCC.sql.NoIdent=sort.nonRedun.140cond.strainCC.sql.NoIdent
abs.sort.nonRedun.140cond.strainCC.sql.NoIdent$Pearson.Correlation.Coefficient=abs(sort.nonRedun.140cond.strainCC.sql.NoIdent$Pearson.Correlation.Coefficient)
save(abs.sort.nonRedun.140cond.strainCC.sql.NoIdent,file="Data/abs.sort.nonRedun.140cond.strainCC.sql.NoIdent.RData")
'

# Get the pwy annotations corresponding to the pcc table
##The following runs about 2 hr on my PC when using 2 cores (2 Rstudio sessions).
'
strain.table=abs.sort.nonRedun.140cond.strainCC.sql.NoIdent
start.time=Sys.time() #Take ? min to finish on my mac
strain1.annot.allPCC.140cond=sapply(strain.table$strain1,FUN = function(strain){ 
  index=which(as.numeric(new.ECK.Pathway$ids)==strain)
  return(new.ECK.Pathway$Data[index])
}
)
end.time=Sys.time()
end.time-start.time
save(strain1.annot.allPCC.140cond,file="Data/strain1.annot.allPCC.140cond.RData")


start.time=Sys.time() 
strain2.annot.allPCC.140cond=sapply(strain.table$strain2,FUN = function(strain){ 
  index=which(as.numeric(new.ECK.Pathway$ids)==strain)
  return(new.ECK.Pathway$Data[index])
}
)
end.time=Sys.time()
end.time-start.time
save(strain2.annot.allPCC.140cond,file="Data/strain2.annot.allPCC.140cond.RData")

#For concenience I just saved them as .RData for loading in case they are accidentally removed
'



# The following uses the dataset that only contains genes annotated in pathways
## Prepare the annotation table (remove the ones that are not in any pathway)
pwy.annotation.table=ECK.Pathway[ECK.Pathway$Data!="Not in any pathway",1:2]

## How many unique gene ids are there? (How many genes are annotated in at least 1 pathway?)
length(unique(pwy.annotation.table$ids)) #893

# Get the pcc table
uniq.id=as.numeric(unique(pwy.annotation.table$ids))
pwy.pcc.140cond=sort.nonRedun.140cond.strainCC.sql.NoIdent[(sort.nonRedun.140cond.strainCC.sql.NoIdent$strain1 %in% uniq.id ) & (sort.nonRedun.140cond.strainCC.sql.NoIdent$strain2 %in% uniq.id) ,]

# Use abs() and sort
pwy.pcc.140cond$Pearson.Correlation.Coefficient=abs(pwy.pcc.140cond$Pearson.Correlation.Coefficient)
pwy.pcc.140cond=pwy.pcc.140cond[order(pwy.pcc.140cond$Pearson.Correlation.Coefficient,decreasing=T),]



## Why does No. of unique id=892 in pwy.pcc.140cond's strain1 and strain2 => This is not wrong. Take a look at the note in Nichols_correlation.R for sort.pcc.sql.NoIdent

# Get the pwy annotations corresponding to the pcc table
strain1.annot.140cond=sapply(pwy.pcc.140cond$strain1,FUN = function(strain){ #This takes like 1 min
  index=which(as.numeric(ECK.Pathway$ids)==strain)
  return(ECK.Pathway$Data[index])
}
)
strain2.annot.140cond=sapply(pwy.pcc.140cond$strain2,FUN = function(strain){ #This takes like 1 min
  index=which(as.numeric(ECK.Pathway$ids)==strain)
  return(ECK.Pathway$Data[index])
}
)


# Get the binary vector (true when at least 1 annotation for strain1 matches strain2's) for each pairwise pcc. Bind to pwy.pcc
bi.pwy.annot=c()
for(i in 1:length(strain1.annot.140cond)){
  bi.pwy.annot[i]=ifelse(sum(strain1.annot.140cond[[i]] %in% strain2.annot.140cond[[i]])>=1,1,0)
}

pwy.pcc.140cond$"At least 1 same annotation"=bi.pwy.annot
