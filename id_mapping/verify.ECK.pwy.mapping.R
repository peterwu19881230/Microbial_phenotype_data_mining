#Goal: how many pathways are there in ECK-pathways mapping?

uniq=unique(EcoCycID.Pwys$`Pathway ID`) # 424

j=0
for(i in 1:length(uniq)){
  indices=which(EcoCycID.Pwys$`Pathway ID` %in% unique[i]) 
  if(length(indices)>2) j=j+1
} #This gives you 300

#If the following pathway is dropped, there are 299, the exact No. of pathways that Nichols et al. had
sum(EcoCycID.Pwys$`Pathway Name`=="tRNA charging") ## 109 genes in this pathway


