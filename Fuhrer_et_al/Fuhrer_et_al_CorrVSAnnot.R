pos=read.table("Data/genomwide_metabolomic/zscore_pos.tsv")
neg=read.table("Data/genomwide_metabolomic/zscore_neg.tsv")


all=cbind(t(pos),t(neg))
str(all)

gene_names=read.xlsx("Data/genomwide_metabolomic/sample_id_zscore.xls",sheetName="Sheet1",header = F)
str(gene_names)
sum(duplicated(gene_names$associated_gene_names)) #no duplications in gene names
names(gene_names)="associated_gene_names"


#assign Nichols' id to Fuhrer's data (I can remap Fuhrer's id to all the annotation sets instead of using Nichols' id to help, but I am not sure it's worth it)
geneName_NicholsID=id_allAttributes[,c("associated_gene_names","ids")] %>% unique
str(geneName_NicholsID)

df=left_join(gene_names,geneName_NicholsID,by="associated_gene_names")
str(df) #increased no. of ID (3807 -> 3828) because of duplications in Nichols

df=df[!is.na(df$ids),]#remove NA (remove genes that are in Fuhrer's but not in Nichols)
str(df) #3453 genes remaining

#Remove duplicates
df_clean=df[! (duplicated(df$associated_gene_names,fromLast = T) | duplicated(df$associated_gene_names,fromLast = F)), ]
str(df_clean)

#See what those duplicates are
df_removed=df[ (duplicated(df$associated_gene_names,fromLast = T) | duplicated(df$associated_gene_names,fromLast = F)), ]
str(df_removed) 
##there should only be 12*2 duplicated strains removed, but here it shows 42 
##-> There are genes with different original names but same gene names in Nichols -> For simplification, I think it's ok to remove all of them


Index=which(gene_names$associated_gene_names %in% df_clean$associated_gene_names  )

clean_phenoDat=all[Index,]

temp=merge(data.frame(associated_gene_names=gene_names$associated_gene_names[Index]),df_clean,by="associated_gene_names")
identical(temp$associated_gene_names,gene_names$associated_gene_names[Index]) #make sure the above merging doesn't change the order of the first column

rownames(clean_phenoDat)=temp$ids

str(clean_phenoDat)

dat=clean_phenoDat
id_pwy=id_allAttributes[,c("ids","Pwy")] %>% unique()
id_pwy=id_pwy[which(id_pwy$ids %in% df_clean$ids),] #remove ids that are not in dat

attribute_list=attr_list(id_pwy$ids,id_pwy$Pwy)


#start.time = Sys.time()
#Fuhrer_test=dist_TF_cumsum_matirxOperation(data=dat,attribute_list=attribute_list,dist_metric=pcc_dist)
#end.time = Sys.time()
#end.time - start.time #Time difference of 8.424094 mins
#save(Fuhrer_test,file="Data/Fuhrer_test.RData")
load("Data/Fuhrer_test.RData")


Fuhrer_test_metrics=confusionMatrix_metrics(cumSums=Fuhrer_test$cumsum)



plot_ROC=function(samples,alpha,size,df){
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(FPR=1-df$specificity[samples],sensitivity=df$sensitivity[samples]),aes(FPR,sensitivity),size=size)+
    geom_abline(slope=1,intercept=0,linetype="dashed")+ #This is the diagnal line
    
    
    labs(x ="FPR (1-specificity)",y="TPR (sensitivity)",aesthetic='custom text')+
    scale_x_continuous(expand = c(0,0),labels=comma)
}


samples=seq(from=0,to=5815755,by=1137); samples[1]=1
plot_ROC(samples=samples,alpha=0.5,size=1,df=Fuhrer_test_metrics)


#Calculate AUC (Area under curve) using a package
library(pROC) #ref: https://stackoverflow.com/questions/4903092/calculate-auc-in-r

predictor=Fuhrer_test$Distance

#pwy
response=Fuhrer_test$sameORnot

start.time = Sys.time()
pwyAUC=auc(response,predictor) #Syntax: auc(response,predictor)
end.time = Sys.time()
end.time - start.time #Time difference of 29.62732 secs

pwyAUC #Area under the curve: 0.5362

#Note: the eaiser thing to redo the above (and finish all combinations of annotation sets) is to calculate each distance for this dataset and map to strain1strain2.... obj for Nichols


#different distances:

#pcc
#spearman
#mi
#qualitative mi
#mhd










