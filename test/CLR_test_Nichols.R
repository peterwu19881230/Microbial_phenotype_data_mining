#Context Likelyhood Relatedness tested on Nichols' data
#ref: https://rdrr.io/bioc/minet/man/clr.html

start.time = Sys.time()
mi_all=build.mim(t(All_Data_NAimputed)) #the default estimator is spearman, but I am not sure it's correct (or acceptable)
end.time = Sys.time()
end.time - start.time 


##Calculate similarity based on CLR:
net = clr(mi_all)

#CLR for all genes
all_pairwise_clr=net %>% as.dist %>% as.vector
summary(all_pairwise_clr) #min=0, max=63.60
mean(all_pairwise_clr)


#CLR for gcv pathway:
load("Data/ECK_1st_table.RData")
columns=GetIndices(c("ECK2898-GCVP","ECK2900-GCVT","ECK2899-GCVH","ECK2548-GLYA"),ECK_1st_table$sorted_ECK_missing_gene_names_added)  
pairwise_clr=net[columns,columns] %>% as.dist %>% as.vector
mean(pairwise_clr) #better than random (population)

#CLR for acrA, acrB, marA, marB, marR (example from Nichols' fig. 4)
##index for acrA, acrB, marA, marB, marR, respectively: 3965, 1278, 2775, 2432, 2693
columns=c(3965, 1278, 2775, 2432, 2693)
pairwise_clr=net[columns,columns] %>% as.dist %>% as.vector
mean(pairwise_clr) #better than random (population)



#gene pairs that have the 2 highest CLR also have high PCC

#Nichols id =c(1950,1987) 
cor(t(All_Data_NAimputed[c(1950,1987),])) #0.93

#Nichols id =c(1950,1987) 
cor(t(All_Data_NAimputed[c(3065,3142),])) #0.81




