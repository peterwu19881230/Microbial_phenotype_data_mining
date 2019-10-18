#Regenerate h-clustering in Nichols' (confirming fig. 1:)

str(rownames(All_Data_NAimputed))

#cluster strains
hclust_original=All_Data_NAimputed %>% pcc_dist %>% hclust(method="complete")

#cluster conditions
hclust_original_condition=t(All_Data_NAimputed) %>% pcc_dist %>% hclust(method="complete")
plot(hclust_original_condition)


cond_in_fig1=c("BILE.0.5....UNSPECIFIED","BILE.1.0....UNSPECIFIED","BILE.2.0....UNSPECIFIED","DEOXYCHOLATE.2.0....UNSPECIFIED","DEOXYCHOLATE.0.5....UNSPECIFIED","DEOXYCHOLATE.0.1....UNSPECIFIED",
              "SDS0.5..EDTA0.5..","SDS1.0..EDTA0.5..","SDS0.5..EDTA0.1..",
              "SDS.0.5....UNSPECIFIED","SDS.1.0....UNSPECIFIED","SDS.2.0....UNSPECIFIED","SDS.3.0....UNSPECIFIED","SDS.4.0....UNSPECIFIED",
              "BENZALKONIUM.10...UNSPECIFIED","BENZALKONIUM.25...UNSPECIFIED",
              "DIBUCAINE.0.4...UNSPECIFIED","DIBUCAINE.0.8...UNSPECIFIED","DIBUCAINE.1.2...UNSPECIFIED",
              "NOVOBIOCIN.30...UNSPECIFIED",
              "TRICLOSAN.0.05...UNSPECIFIED",
              "ACRIFLAVINE.10...UNSPECIFIED",
              "ETHIDIUMBROMIDE.50...UNSPECIFIED",
              "ACRIFLAVINE.2...UNSPECIFIED",
              "ETHIDIUMBROMIDE.10...UNSPECIFIED","ETHIDIUMBROMIDE.2...UNSPECIFIED",
              "PROPIDIUMIODIDE.20...UNSPECIFIED","PROPIDIUMIODIDE.50...UNSPECIFIED",
              "MINOCYCLINE.0.2...UNSPECIFIED","MINOCYCLINE.0.5...UNSPECIFIED","MINOCYCLINE.1.0...UNSPECIFIED",
              "PUROMYCIN.25...UNSPECIFIED","PUROMYCIN.5...UNSPECIFIED",
              "PYOCYANIN.10.0...UNSPECIFIED",
              "MITOMYCINC.0.1...UNSPECIFIED",
              "STREPTONIGRIN.0.5...UNSPECIFIED"
)

as.dendrogram(hclust_original_condition) %>% set("by_labels_branches_col", value = cond_in_fig1) %>% plot(main = "Conditions")
##Not exactly the same as Nichols' says in fig.1, but most conditions cluster together 
##(Reason: might be because that I imputed the data or there's slight differences in ways I and Nasos treat pairwise pcc)


#fig. C upper strains
ECKs=c("ECK3622","ECK3617","ECK3618","ECK3620","ECK3621","ECK3610","ECK3609","ECK0223","ECK3042")

#Get the id I created for Nichols' ECK
ids=lapply(ECKs,FUN=function(ECK){
  id_allAttributes$ids[id_allAttributes$ECK==ECK] %>% unique
}) %>% unlist

ids=as.numeric(ids)


hclust_original$labels=1:3979 #set the label for this h-clust. I don't know why hclust() doesn't use the original rownames in All_Data_NAimputed
as.dendrogram(hclust_original) %>% set("by_labels_branches_col", value = ids) %>% plot(main = "Strains")
#I see 2 major branches and they are far away from each other

#verify with heatmap:
dat=All_Data_NAimputed[ids,cond_in_fig1]
Heatmap(dat,col=colorRamp2(c(-10,0,10),c("green","black","red")),cluster_rows=FALSE,cluster_columns=FALSE)
##Note: grey area in fig1 are where they have NA. They look black here because I have imputed the data.


#fig. C middle strains
ECKs=c("ECK2272","ECK2281","ECK2273","ECK2282","ECK2274","ECK2276","ECK2270","ECK2278","ECK2271","ECK2279")

#Get the id I created for Nichols' ECK
ids=lapply(ECKs,FUN=function(ECK){
  id_allAttributes$ids[id_allAttributes$ECK==ECK] %>% unique
}) %>% unlist

ids=as.numeric(ids)

hclust_original$labels=1:3979 #set the label for this h-clust. I don't know why hclust() doesn't use the original rownames in All_Data_NAimputed
as.dendrogram(hclust_original) %>% set("by_labels_branches_col", value = ids) %>% plot(main = "Strains")
#I see one major cluster. They do cluster together


#verify with heatmap:
dat=All_Data_NAimputed[ids,cond_in_fig1]
Heatmap(dat,col=colorRamp2(c(-10,0,10),c("green","black","red")),cluster_rows=FALSE,cluster_columns=FALSE)
##Looks identical by my eyes.

#Lower one only contains 2 strains and the conditions are all differnt. Reluctant to do it.




