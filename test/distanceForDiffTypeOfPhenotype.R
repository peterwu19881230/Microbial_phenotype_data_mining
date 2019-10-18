# Quantitative phenotypic profiles after removing non-significant phenotypes

#Phenotype data used: sigQuantitative_Data_324cutoff

#remove strains that have no significant phenotypes
TF=apply(sigQuantitative_Data_324cutoff,1,FUN=function(row){
  sum(row!=0)>0
})

sum(TF) #2210

dat=sigQuantitative_Data_324cutoff[TF,]

##PCC based distance
pcc= ( 1-abs(cor(t(dat))) ) %>% as.dist %>% melt_dist

##Spearman based distance
spearman=( 1-abs(cor(t(dat),method="spearman")) ) %>% as.dist %>% melt_dist

##Euclidean distance
euclidean=get_dist(dat,method="euclidean") %>% melt_dist


##In strain1strain2_allDistances.R I have demonstrated that the order of the row for each distance should be the same,
##so here I will just bind them without row name matching
distances__sigQuantitative_Data_324cutoff=data.frame(strain1=as.integer(pcc$`Object 1`),
                                                     strain2=as.integer(pcc$`Object 2`),
                                                     pcc=pcc$Distance,
                                                     spearman=spearman[,3],
                                                     euclidean=euclidean[,3])








# Quantitative phenotypic profiles after collapsing conditions (most significant phenotypes are used among different concentrations of conditions and if there's a conflict,there will be no phenotype (fitness = 0) )

#Phenotype data used: sigQuantitative_Data_324cutoff_condCollapsed

TF=apply(sigQuantitative_Data_324cutoff_condCollapsed,1,FUN=function(row){ 
  sum(row!=0)>0 
}) #This TF vector is different than the previous one because after collapsing conditions, 
   #some strains that have at least 1 significant phenotypes will not have any (I remove conflicting phenotypes)
   # from sum(TF) I know I only have to remove 1 additional strain which has such conflict

sum(TF) #2209

dat=sigQuantitative_Data_324cutoff_condCollapsed[TF,]

##PCC
pcc= ( 1-abs(cor(t(dat))) ) %>% as.dist %>% melt_dist

##Spearman
spearman=( 1-abs(cor(t(dat),method="spearman")) ) %>% as.dist %>% melt_dist

##Euclidean
euclidean=get_dist(dat,method="euclidean") %>% melt_dist


##In strain1strain2_allDistances.R I have demonstrated that the order of the row for each distance should be the same,
##so here I will just bind them without row name matching


distances__sigQuantitative_Data_324cutoff_condCollapsed=data.frame(strain1=as.integer(pcc$`Object 1`),
                                                     strain2=as.integer(pcc$`Object 2`),
                                                     pcc=pcc$Distance,
                                                     spearman=spearman[,3],
                                                     euclidean=euclidean[,3])
















