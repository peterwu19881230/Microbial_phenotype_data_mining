#Box plots for pathways (Repeat Curtis' avg pcc experiment)


##Make pathway_ids (another format of ECK.Pathway_table but: 1. doesn't contain genes not involved in pathways 
in_pwy_table=ECK.Pathway_table[!(is.na(ECK.Pathway_table$PathwayID)),]
in_pwy_table$ids=as.numeric(in_pwy_table$ids)
##(!)In this experiment the pathway annotations shouldn't be limited to what was used in Nichols' fig.S1, so I use all of them

id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique
id_pwy=id_pwy[!is.na(id_pwy$Pwy),]

pathway_ids=attr_list(id_pwy$Pwy,id_pwy$ids)




cors_for_pwys=data.frame()
count=1
for(i in 1:length(pathway_ids)){
  if(length(pathway_ids[[i]]) %in% c(1,2)){ #There are 7 pwy that has only 1 gene, 10 pwy that has only 2 genes (They are pathways that have at least 3 genes, but not all of them could be found in Nichols')
    print("The following only has 1 or 2 gene:")
    print(pathway_ids[i])
    next
  }
  comb=t(combn(pathway_ids[[i]],2))
  print(dim(comb)[1])
  
  pcc=c()
  for(j in 1:dim(comb)[1]){
    pcc[j]=cor_strains[comb[j,1],comb[j,2]] 
  }
  cors_for_pwys=rbind(cors_for_pwys,c(names(pathway_ids[i]),mean(pcc),sd(pcc),length(pathway_ids[[i]])),stringsAsFactors=F)
  
  count=count+1
}
names(cors_for_pwys)=c("pwy","avg_pcc","std_pcc","no_gene_used")
cors_for_pwys$avg_pcc=as.numeric(cors_for_pwys$avg_pcc) 
cors_for_pwys$std_pcc=as.numeric(cors_for_pwys$std_pcc)
cors_for_pwys$no_gene_used=as.numeric(cors_for_pwys$no_gene_used) 

##Add the column for p-values
#attach(cors_for_pwys)
#cors_for_pwys$p_value=2*pt(-abs( avg_pcc/(std_pcc/sqrt(no_gene_used)) ),df=no_gene_used-1)

#save(cors_for_pwys,file="Data/sourced/cors_for_pwys.RData")



##Do absolute value
abs_cors_for_pwys=data.frame()
abs_cors_for_pwys_all=list()
count=1
for(i in 1:length(pathway_ids)){
  if(length(pathway_ids[[i]]) %in% c(1,2)){ #There are 7 pwy that has only 1 gene, 10 pwy that has only 2 genes (They are pathways that have at least 3 genes, but not all of them could be found in Nichols')
    print("The following only has 1 or 2 gene:")
    print(pathway_ids[i])
    next
  }
  comb=t(combn(pathway_ids[[i]],2))
  print(dim(comb)[1])
  
  abs_pcc=c()
  for(j in 1:dim(comb)[1]){
    abs_pcc[j]=cor_strains[comb[j,1],comb[j,2]] %>% abs
  }
  abs_cors_for_pwys=rbind(abs_cors_for_pwys,c(names(pathway_ids[i]),mean(abs_pcc),sd(abs_pcc),length(pathway_ids[[i]])),stringsAsFactors=F)
  
  abs_cors_for_pwys_all[[count]]=abs_pcc
  names(abs_cors_for_pwys_all)[count]=names(pathway_ids[i])
  
  count=count+1
}
names(abs_cors_for_pwys)=c("pwy","avg_pcc","std_pcc","no_gene_used")
abs_cors_for_pwys$avg_pcc=as.numeric(abs_cors_for_pwys$avg_pcc) 
abs_cors_for_pwys$std_pcc=as.numeric(abs_cors_for_pwys$std_pcc)
abs_cors_for_pwys$no_gene_used=as.numeric(abs_cors_for_pwys$no_gene_used) 

##Add the column for p-values
#attach(abs_cors_for_pwys)
#abs_cors_for_pwys$p_value=2*pt(-abs( (avg_pcc-mean(abs_unique_cor_strains))/(std_pcc/sqrt(no_gene_used)) ),df=no_gene_used-1)
##(!)Not sure if I should use mean() as mu_h0
##(!)Not sure whether I should use no. of genes used or no. of total genes (some genes are not in Nichols') in that pathway

#sum(p_value<0.05)/280
#save(abs_cors_for_pwys,file="Data/sourced/abs_cors_for_pwys.RData")



# boxplot for abs_pcc of 280 pathways

##old way
boxplot(abs_cors_for_pwys_all)


##ggplot way
df.abs_cors_for_pwys_all=data.frame(pwy=rep(names(abs_cors_for_pwys_all),sapply(abs_cors_for_pwys_all,length)),abs_pcc=unlist(abs_cors_for_pwys_all))
##Ref: https://stackoverflow.com/questions/12639501/convert-list-to-data-frame-while-keeping-list-element-names  

random_expectation=mean(cor_strains %>% as.dist %>% abs %>% mean)

ggplot(df.abs_cors_for_pwys_all,aes(x=reorder(pwy,-abs_pcc,FUN=median),y=abs_pcc))+theme_minimal()+#How to reorder bars based on bar values: https://sebastiansauer.github.io/ordering-bars/
  scale_y_continuous(limits = c(0, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 5))+
  geom_boxplot()+
  geom_hline(yintercept = random_expectation,linetype="dashed",colour="blue")





#boxplot for n=c(3:22,27,28,31,41,47,48)
#I am currently using this: no matter how many genes are used, 
# pathways of same No. of genes defined by EcoCyc should be put into 1 box 

#Use n gene pathway according to Pathways.col -> box plot for n=1,2,....max(n) gene pathway in Pathways.col
unique_n=ECK.Pathway_table$NoGeneInPathway_byPathwaysCol %>% as.character %>% as.numeric %>% unique %>% sort ##sort() removes the NA
##Notes about converting factor to numeric: https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information

cors_for_nGenePath=list()
count2=1
for(i in unique_n[-1]){ #-1 is to remove 1 gene pathway 
  id_for_nGenePath=ECK.Pathway_table$ids[!is.na(ECK.Pathway_table$NoGeneInPathway_byPathwaysCol)&ECK.Pathway_table$NoGeneInPathway_byPathwaysCol==i]
  
  unique_Path_forNoGenes=ECK.Pathway_table[as.numeric(id_for_nGenePath),]$PathwayID  %>% unique
  unique_Path_forNoGenes=unique_Path_forNoGenes[-1] #remove the NA
  
  count1=1
  cors_for_nGenePath[[count2]]=numeric(0) #Reluctantly I have to later grow a vector 
  for(pathway in unique_Path_forNoGenes){
    id_for_path=ECK.Pathway_table$ids[which(ECK.Pathway_table$PathwayID==pathway)]; if(length(id_for_path)<2) next #skip if there are not at least 2 genes found in Nichols'
    pair=id_for_path %>% as.numeric %>% combn(m=2) %>% t
    
    for(j in 1:dim(pair)[1]){
      cors_for_nGenePath[[count2]][count1]=cor_strains[pair[j,1],pair[j,2]] %>% abs
      count1=count1+1
    }
    
  }
  count2=count2+1
}
names(cors_for_nGenePath)=unique_n[-1] #-1 is to remove 1 gene pathway


boxplot(cors_for_nGenePath,main="Box plot for n gene pathways",xlab="No. of genes in the pathway (defined in Pathways.col)",ylab="|PCC|")
random_abs_pcc=cor_strains %>% as.dist %>% as.numeric %>% abs %>% mean
abline(h=random_abs_pcc,lty=2,col="blue")

#mean for n gene pathways
sapply(cors_for_nGenePath,mean)