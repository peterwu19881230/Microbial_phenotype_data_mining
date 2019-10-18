#Look at fitness scores for auxotrophs under minimal media

##auxotrophs: https://docs.google.com/spreadsheets/d/1lGMs6E-u8-KOrn-_7FBYkdqvwAjzIAdu/edit#gid=255348281
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDir)
auxotroph=read.table("auxotrophs.txt") %>% as.matrix %>% as.character #note: this .txt file is copied-pasted from the above gsheet

#get Nichols' id for those auxotrophs
auxotroph_id=sapply(auxotroph,FUN=function(original_name){
  index=which(ECK_1st_table$`Original Name`==original_name)
  
  return(ECK_1st_table$ids[index] %>% as.numeric)
  
}) #the result is a named vector


#get fitness scores of the auxotrophs under minimal media
TF=( names(uniqueChemIndex) %in% c("NH4Cl (MOPS)","Iron excess-FeSO4","Iron starvation-FeSO4","Acetate (M9)",
                                   "Glucosamine (M9)","Glucose (M9)","Glycerol (M9)","Maltose (M9)","N-acetyl Glucosamine","Succinate (M9)") )
minimal_uniq=names(uniqueChemIndex)[TF]

cond_indices_minimal=uniqueChemIndex[minimal_uniq] %>% as.numeric

auxotroph_minimal_fitness=All_Data[auxotroph_id,cond_indices_minimal] %>% as.matrix %>% as.numeric


##Result by analysis
png(filename=paste(currentDir,"auxotroph_minimal_fitness.png",sep="/"))
hist(auxotroph_minimal_fitness,breaks=100,xlab='fitness',main="Distribution of fitness of 102 auxotrophs under 10 minimal media")
dev.off()

summary(fitness_auxotroph_minimal)


