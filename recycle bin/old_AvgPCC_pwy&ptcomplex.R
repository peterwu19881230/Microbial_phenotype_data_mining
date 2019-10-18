#Goal: Determine the significance of avg pcc of pwys and ptcomplexes. 

#Part I: Pwys (Curtis has done part of this)
#Part II: Ptcomplexes
#Part III: modules from KEGG 


#Pwys (Curtis has done part of this)
##Note: Here I use just pwys in Nichols' fig.S1

##Make pathway_ids (another format of ECK.Pathway_table but: 1. doesn't contain genes not involved in pathways 2. contains only genes in Nichls' fig. S1)
in_pwy_table=ECK.Pathway_table[!(is.na(ECK.Pathway_table$PathwayID)),]
in_pwy_table$ids=as.numeric(in_pwy_table$ids)
##(!)In this experiment the pathway annotations shouldn't be limited to what was used in Nichols' fig.S1, so I use all of them

unique_pwy_annot=unique(in_pwy_table$PathwayID)

pathway_ids=list()
for(i in 1:length(unique_pwy_annot)){
  TFvec=(in_pwy_table$PathwayID==unique_pwy_annot[i])
  ids=in_pwy_table$ids[TFvec]
  pathway_ids[[i]]=ids
}
names(pathway_ids)=unique_pwy_annot

##Remove genes that are not in fig. S1
s1.pwy[!(s1.pwy %in% unique_pwy_annot)] ## 2 pathways being removed: "COA-PWY" "ARABCAT-PWY"
pathway_ids=pathway_ids[unique_pwy_annot %in% s1.pwy] #s1.pwy (contains 299 pwys) is from figS1.R (297 pathways. Note that from Pathways.col there are 299 pathways used if 1 and 2 gene pathways are removed from Pathways.col)
#save(pathway_ids,file="Data/pathway_ids.RData")


#The following gathers all avg from genes in pwy
cor_in_pwy=c()
for(i in 1:length(pathway_ids)){
  if(length(pathway_ids[[i]])==1){ #There are 7 pwy that has only 1 gene
    print("The following only has 1 gene:")
    print(pathway_ids[i])
    next
  }
  comb=t(combn(pathway_ids[[i]],2))
  print(dim(comb)[1])
  
  for(j in 1:dim(comb)[1]){
    cor_in_pwy=c(cor_in_pwy,cor_strains[comb[j,1],comb[j,2]])
  }
}


#Distribution of mean pccs
hist(cor_in_pwy,main="Distribution of mean pccs of 290 pathways",xlab="Pearson Correlation Coefficient")
summary(cor_in_pwy)



##Violin plot => after updating to R 3.5 the following block of code returns error: Error: package ‘sm’ could not be loaded

#library(vioplot) #https://cran.r-project.org/web/packages/vioplot/vioplot.pdf
##A reference (not sure if it's the case for all the violin plots): https://datavizcatalogue.com/methods/violin_plot.html
#vioplot(cor_in_pwy,sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient,names=c("Gene pairs in pathways","All gene pairs")) #In pathways here means in pathways used in fig.S1
#mtext("Violin plot of pcc in pwy and for all",side=3,1,0,cex=1.5,col="black") 
#mtext("Pearson correlation coefficient",side=2,2,0,cex=1,col="black") 


#Negative control (average of all pcc)
mean(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient) ##0.0007193512

#Part II: Ptcomplexes
##I used inptcomplex_strain1strain2.pcc.samePTcomplex_Annot defined in another script
pcc_in_ptcomplex=inptcomplex_strain1strain2.pcc.samePTcomplex_Annot$Pearson.Correlation.Coefficient[
  inptcomplex_strain1strain2.pcc.samePTcomplex_Annot$bi.ptcomplex.annot==1
]

#Distribution of mean pccs
hist(pcc_in_ptcomplex,main="Distribution of mean pccs of 206 ptcomplexes",xlab="Pearson Correlation Coefficient")
summary(pcc_in_ptcomplex)




#Part III:  modules from KEGG => For this part I only do abs_pcc because the non-abs doesn't seem to be useful
abs_pcc_in_KEGGmod=1-keggInModules_quantitative[[1]]$pcc[keggInModules_quantitative[[1]]$sameORnot==1] 
# |pcc| = 1- pcc_distance=> have to use non-abs first

#Distribution of mean pccs
hist(abs_pcc_in_KEGGmod,main="Distribution of mean pccs of modules",xlab="Pearson Correlation Coefficient")
summary(abs_pcc_in_KEGGmod)




##Violin plot (The shape is a little different than using vioplot() but overall they look similar. Is it because ggplot() gives higher resolution?)
df1 = data.frame(pwy=abs(cor_in_pwy))
df2 = data.frame(ptcom=abs(pcc_in_ptcomplex))
df3=data.frame(kegg=abs_pcc_in_KEGGmod)
df4=data.frame(all=abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient))


ggplot() +
  geom_violin(data = df1,aes("Gene pairs in pathways",pwy)) +
  geom_boxplot(data = df1,aes("Gene pairs in pathways",pwy),width=0.1,outlier.shape = NA)+
  geom_violin(data = df2,aes("Gene pairs in protein complexes",ptcom))+
  geom_boxplot(data = df2,aes("Gene pairs in protein complexes",ptcom),width=0.1,outlier.shape = NA)+
  geom_violin(data = df3,aes("Gene pairs in KEGG modules",kegg))+
  geom_boxplot(data = df3,aes("Gene pairs in KEGG modules",kegg),width=0.1,outlier.shape = NA)+
  geom_violin(data = df4,aes("All gene pairs",all)) +
  geom_boxplot(data = df4,aes("All gene pairs",all),width=0.1,outlier.shape = NA)+
  scale_x_discrete("",limits=c("Gene pairs in pathways","Gene pairs in protein complexes","Gene pairs in KEGG modules","All gene pairs"))+ 
  #I want the x axis to be empty. And if I don't use this, the order is not right
  scale_y_continuous("Abs of Pearson Correlation Coefficient")



##My old way of doing Violin plot
'
vioplot(pcc_in_ptcomplex,sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient,names=c("Gene pairs in protein complexes","All gene pairs")) #Gene pairs in protein complexes means: average pcc of genes in 206 ptcomplexes
mtext("Violin plot of pcc in ptcomplex and for all",side=3,1,0,cex=1.5,col="black") 
mtext("Pearson correlation coefficient",side=2,2,0,cex=1,col="black") 

#Violin plot for Pathways + Protein complexes + All genes
vioplot(cor_in_pwy,pcc_in_ptcomplex,sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient,names=c("Gene pairs in pathways","Gene pairs in protein complexes","All gene pairs")) 
mtext("Violin plot of pcc in pwy, ptcomplex and for all",side=3,1,0,cex=1.5,col="black") 
mtext("Pearson correlation coefficient",side=2,2,0,cex=1,col="black")

#(After taking absolute value) Violin plot for Pathways + Protein complexes + All genes
vioplot(abs(cor_in_pwy),abs(pcc_in_ptcomplex),abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient),names=c("Gene pairs in pathways","Gene pairs in protein complexes","All gene pairs")) 
mtext("Violin plot of pcc in pwy, ptcomplex and for all",side=3,1,0,cex=1.5,col="black") 
mtext("Pearson correlation coefficient",side=2,2,0,cex=1,col="black") 
'

mean(abs(cor_in_pwy))
mean(abs(pcc_in_ptcomplex))
mean(abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient))

#Suggested by Dr. Cai, I can use Mann-Whitney U test. And I am thinking that I might want to use 1-sided test
#Ref: https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-the-differences-between-one-tailed-and-two-tailed-tests/

## An example:
##wilcox.test(c(10,10,10,10,10),c(1,2,3,4,5),alternative="greater")

wilcox.test(abs(cor_in_pwy),abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient),alternative="greater") #p-value < 2.2e-16
wilcox.test(abs(pcc_in_ptcomplex),abs(sort.pcc.sql.NoIdent$Pearson.Correlation.Coefficient),alternative="greater") #p-value < 2.2e-16


#Repeat Curtis' avg pcc experiment
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
#Should I do this for each box?: no matter how many genes are used, 
# pathways of same No. of genes defined by EcoCyc should be put into 1 box 

#Have to ask Dr. Hu about this: use n gene pathway according to Pathways.col -> box plot for n=1,2,....max(n) gene pathway in Pathways.col
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


#boxplot
#Note: up to this point there are only 280 pathways being used
ordered_cors_for_pwys=abs_cors_for_pwys[order(abs_cors_for_pwys$no_gene_used,-abs_cors_for_pwys$avg_pcc),]
##Ref about using order(a,-1): https://stackoverflow.com/questions/7793295/how-to-order-a-data-frame-by-one-descending-and-one-ascending-column
##Note: I feel that the part about rev() is BSing. rev() doesn't really help (I wonder if the person really had run it)

#have to turn no_gene_used into a factor in order to get proper color filling
ordered_cors_for_pwys$no_gene_used=as.factor(ordered_cors_for_pwys$no_gene_used)
ordered_cors_for_pwys$pwy=as.factor(ordered_cors_for_pwys$pwy)

#Lost the way I plot this using different colors for genes in n gene pathways, but I guess it's not important anymore 
tab=left_join(df.abs_cors_for_pwys_all,ordered_cors_for_pwys[,c("pwy","no_gene_used")],by="pwy")
tab=arrange(tab,as.numeric(no_gene_used),-abs_pcc)
tab$pwy=factor(tab$pwy,levels=unique(tab$pwy))#Do this to prevent automatic x axis reordering

##ref: https://stackoverflow.com/questions/43877663/order-multiple-variables-in-ggplot2
##ref: http://www.sthda.com/english/wiki/ggplot2-facet-split-a-plot-into-a-matrix-of-panels
##Also, I am separating them into 2 lines

ggplot(tab[tab$no_gene_used %in% c(3,4,5,6,7),], aes(x=pwy, y=abs_pcc,fill=no_gene_used)) + 
  #ggtitle("Average |PCC| for gene pairs in EcoCyc Pathways")+
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)


ggplot(tab[!tab$no_gene_used %in% c(3,4,5,6,7),], aes(x=pwy, y=abs_pcc,fill=no_gene_used)) + 
  theme(plot.title= element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),legend.position="none",
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20))+
  facet_grid(. ~ no_gene_used, scales="free_x", space="free_x")+
  geom_boxplot()+ylab("")+xlab("")+
  #ylab("|PCC|")+xlab("EcoCyc Pathway Annotations")+
  geom_hline(yintercept = random_expectation,linetype="longdash",colour="black",size=1,alpha=0.5)




