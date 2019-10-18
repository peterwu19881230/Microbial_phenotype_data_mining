#Goal: create mean(mean(radom strainCC from j=3,4,5...)). j represents No. of genes in pathways

##Get the values of j
Pathways=ECK.Pathway$Data[ECK.Pathway$Data!="Not in any pathway"] ##extract pathway frequency and remove the "Not in any pathway"
j=table(Pathways) %>% sort %>% as.numeric %>% unique 
j=j[j!=1 & j!=2] #remove 1 gene and 2 gene pathways

##Sample j strains from All_Data. Average the pairwise pcc.
##This takes about 10 sec.
'
random.strainCC=c()
k=0
for(i in j){ ##j=(3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 20 21 22 26 27 28 31 41 48)
  k=k+1
  for(n in 1:100){ ##repeat 100 times for sampling. 
  indices=sample(1:3979,i)
  mat=All_Data[indices,]
  random.strainCC[k]=cor(t(mat),use="pairwise.complete.obs", method="pearson") %>% mean 
  }
} 
'
#result: random.strainCC=0.31839590 0.36290688 0.22348155 0.14289477 0.13858540 0.12503695 0.19356493 0.10724866 0.09479564 0.08776620
#0.06965361 0.10059420 0.08588408 0.06719162 0.07673919 0.05604816 0.04632447 0.04610556 0.04656969 0.03294944
#0.03824613 0.03749123 0.04389637 0.02093455 0.02074612


###names(random.strainCC)=corresponding j's values
names(random.strainCC)=as.character(j)


##Call the average random values by: random.strainCC["j"]














