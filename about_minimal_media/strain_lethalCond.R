#look at fitness scores for CE genes under lethal conditions


currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)


#CE genes & condition: https://docs.google.com/spreadsheets/d/0BybqJiJOk5KZZXk5Nk1ERFd1T2c/edit#gid=69631418
strain_lethalCond=read.table(paste(currentDir,"strain_lethalCond.txt",sep="/"),sep="\t") #this .txt file is copied-pasted from the above gsheet

strain_lethalCond[,1] %>% unique %>% length #no. of unique CE genes

CE_original_name=strain_lethalCond[,1] 

tab_to_join=ECK_1st_table[,c("Original Name","ids")]
names(tab_to_join)=c("Original_Name","ids")

CE_original_name_id=left_join(data.frame('Original_Name'=CE_original_name),tab_to_join,by="Original_Name") %>% unique 
##=>Note: some names match to more than 1 id (this is correct. Nichols has duplicated strains)



##get the indices for those conditions in All_Data
cond=read.csv("Data/allData.csv",header=F)[1,-1] %>% as.character
trimmed_cond=base::trimws(cond) #trim conditions to get rid of white spaces at the end of the condition string to avoid further matching problems

cond_indices=match(base::trimws(strain_lethalCond[,3]),trimmed_cond) #I also have to trim condition strings from what Debby gave me

strain_lethalCond_cond_indices=cbind(strain_lethalCond,cond_indices)
names(strain_lethalCond_cond_indices)=c("Original_Name","cond","original_cond","cond_indices")

temp=left_join(CE_original_name_id,strain_lethalCond_cond_indices[,c(1,4)],by="Original_Name")

auxotroph_ids__lethal_conds=temp[,-1]
auxotroph_ids__lethal_conds$ids=as.integer(auxotroph_ids__lethal_conds$ids)
auxotroph_ids__lethal_conds$cond_indices=as.integer(auxotroph_ids__lethal_conds$cond_indices)


dat=as.matrix(All_Data)
rownames(dat)=as.character(1:3979)
colnames(dat)=as.character(1:324)
melted_dat=melt(dat)
names(melted_dat)=c("ids","cond_indices","fitness")

auxotroph_ids__lethal_conds__fitness=left_join(auxotroph_ids__lethal_conds,melted_dat,by=c("ids","cond_indices"))



##Result by analysis
png(filename=paste(currentDir,"auxotroph_ids__lethal_conds__fitness.png",sep="/"))
hist(auxotroph_ids__lethal_conds__fitness$fitness,breaks=30,xlab="fitness",main="Distribution of fitness scores for lethal phenotypes")
dev.off()

summary(auxotroph_ids__lethal_conds__fitness$fitness)





