#Average PCC for each treatment. For treatments, I don't think I should use absolute value for PCC

avg_std=matrix(numeric(114*2),nrow=114)
for(i in 1:length(uniqueChemIndex)){
  
  index=uniqueChemIndex[[i]]
  
  ##Filter out treatment with only 1 concentration
  if(length(index)==1) {avg_std[i,]=NA; next}
  
  #I use All_Data_NAimputed here. Only 2234/1289196 (0.17%) are NA => probably doesn't matter
  avg=cor(All_Data_NAimputed[,index]) %>% as.dist %>% mean
  std=cor(All_Data_NAimputed[,index]) %>% as.dist %>% sd
  
  avg_std[i,]=c(avg,std)
  
}

rownames(avg_std)=names(uniqueChemIndex)
colnames(avg_std)=c("avg","std")
avg_std=as.data.frame(avg_std)


#some stats:

#simple histogram
hist(avg_std[,1])
summary(avg_std[,1])

#ggplot version (with standard deviation as error bars)
filtered_data=avg_std[is.na(avg_std$avg)==F,] #Filter NA first

ggplot(filtered_data,aes(x=avg))+
  geom_histogram(color="black", fill="white")




#Barplot for each treatment
filtered_data=avg_std
filtered_data$avg[is.na(filtered_data$avg)]=0 #give 0 for those NA so they are shown as nothing on the barplot

filtered_data$std[is.na(filtered_data$std)]=0 
#If there are only 2 conditions, there is only 1 PCC => no std
#It should be NA. I gave 0 so no error bars are shown for these treatments

filtered_data=filtered_data[order(filtered_data$avg,decreasing = T),] #sort the data based on avg PCC

##clean one name because it has "\n" =>for it to be properly displayed on the bar graph
which(rownames(filtered_data)=='Calcofluor (F3543    
Fluorescent Brightener 28) ')

rownames(filtered_data)[95]='Calcofluor (F3543 Fluorescent Brightener 28)'

ggplot(filtered_data,aes(x=reorder(rownames(filtered_data),-avg),y=avg))+theme_minimal()+#How to reorder bars based on bar values: https://sebastiansauer.github.io/ordering-bars/
  scale_y_continuous(limits = c(0, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_bar(stat="identity",color="black",fill="#CCE5FF",width = 0.5)+ #Read this to understand: stat="identity": https://www.rdocumentation.org/packages/ggplot2/versions/1.0.1/topics/geom_bar
  geom_errorbar(aes(ymin=avg-1*filtered_data$std,ymax=avg+filtered_data$std),width=0.2)


##Average PCC across all 324 condition pairs
average_all=cor_conditions %>% as.dist %>% mean #0.02503471
sort.ConditionPCC.sql.NoIdent$Distance %>% mean #0.02503471 (This is just to doulbe check that the above line works)



