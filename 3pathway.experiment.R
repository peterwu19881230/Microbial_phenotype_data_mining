##3 pathways experiment (total 17 genes involved)
dta=ECK.Pathway[order(ECK.Pathway$Data),]
sub.dta=dta[1:17,]  ## Pwy1: 1CMET2-PWY  Pwy2: 2PHENDEG-PWY Pwy3: 4AMINOBUTMETAB-PWY

temp=All_Data[as.numeric(sub.dta$ids),]


sig.3.pwy=All_Data[as.numeric(sub.dta$ids),]

##Use this to change labels to pathways
rownames(sig.3.pwy)=c("pw1-1","pw1-2","pw1-3","pw1-4","pw1-5","pw1-6","pw1-7","pw1-8","pw1-9","pw1-10",
                      "pw2-1","pw2-2",
                      "pw3-1","pw3-2","pw3-3","pw3-4","pw3-5")


library(dplyr)
dis=sig.3.pwy %>% pcc_dist ##pcc.dist is a function defined by me


hc=dis %>% hclust 
plot(hc)

input=data.frame(id=1:17,label=c("pw1","pw1","pw1","pw1","pw1","pw1","pw1","pw1","pw1","pw1",
                                 "pw2","pw2",
                                 "pw3","pw3","pw3","pw3","pw3"))


Three.pwy.bottom=bottomUpExp(hc,input)


##https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html
library(dendextend)
dd=as.dendrogram(hc)
dd %>% set("branches_k_color", k = 3) %>% plot(main = "3 pathways")

## A demo of how to set colors for a specific label
hc %>% as.dendrogram %>% set("by_labels_branches_col", value = c("pw1-4","pw1-10")) %>% plot


#My idea to determine the optimum No. of clusters
probs=c()
for(i in 1:(dim(sub.dta)[1]-1)){
  k=i
  tree=cutree(hc, k = k) 
  
  
  #Swap the values and labels of "tree"
  temp=names(tree)
  names(temp)=tree
  tree=temp
  rm(temp)
  
  
  ##put the genes in the same cluster into the same sublist of a list
  hc.k.list=list()
  for(j in 1:(names(tree) %>% unique %>% length)){
    hc.k.list[[j]]=tree[names(tree)==as.character(j)] 
    
    ##Use the calculated random strainCC values
    
  }
  
  
  ###dis %>% as.numeric %>% mean
}


