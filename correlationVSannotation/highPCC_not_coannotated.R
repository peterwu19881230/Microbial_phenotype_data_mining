#list of high pcc strain pairs in Nichols that are not co-annotated

#Not co-annotated
anyCoAnnotation=ifelse(rowSums(sigPheno_strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules")])>0,
                       T,F)

noCoannotation=sigPheno_strain1strain2_allAnnotations_allDistances[!anyCoAnnotation,]
temp=noCoannotation[order(noCoannotation$pcc),c("strain1","strain2","pcc")] #pcc here means pcc based distance (1-|PCC|)

id_geneName=(id_allAttributes[,c("ids","associated_gene_names")] %>% unique)[order(as.integer(id_geneName$ids)),] 

temp$strain1=(id_geneName$associated_gene_names %>% as.character)[temp$strain1]
temp$strain2=(id_geneName$associated_gene_names %>% as.character)[temp$strain2]
temp$pcc=1-temp$pcc #PCC based distance is changed to |PCC|
result=temp
head(result)                                                                              
dim(result)
(result$pcc>0.6) %>% sum #2034
(result$pcc>0.7) %>% sum #1320
(result$pcc>0.8) %>% sum #693
(result$pcc>0.9) %>% sum #152

#Everything
temp=sigPheno_strain1strain2_allAnnotations_allDistances[order(sigPheno_strain1strain2_allAnnotations_allDistances$pcc),c("strain1","strain2","pcc")] #pcc here means pcc based distance (1-|PCC|)
id_geneName=(id_allAttributes[,c("ids","associated_gene_names")] %>% unique)[order(as.integer(id_geneName$ids)),] 

temp$strain1=(id_geneName$associated_gene_names %>% as.character)[temp$strain1]
temp$strain2=(id_geneName$associated_gene_names %>% as.character)[temp$strain2]
temp$pcc=1-temp$pcc #PCC based distance is changed to |PCC|
result_everything=temp



#plot |PCC| against no. of pairs not co-annotated + everything
result_0.6=result[result$pcc>=0.6,]
result_everything_0.6=result_everything[result_everything$pcc>=0.6,]

size=1
alpha=0.3
ggplot()+theme_minimal()+
  geom_line(data = result_0.6,aes(x=pcc,y=1:length(pcc),color="Not co-annotated"),size=size)+
  geom_line(data = result_everything_0.6,aes(x=pcc,y=1:length(pcc),color="All"),size=size)+
  xlim(0.6,1)+scale_x_reverse()+
  scale_colour_manual("",values=c("Not co-annotated"=alpha("blue",alpha),"All"=alpha("black",alpha)))+
  theme(axis.text = element_text(size=15),
        legend.text=element_text(size=15))+
  labs(x="",y="")











