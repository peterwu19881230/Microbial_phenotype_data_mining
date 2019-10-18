#Does genes in TCA cycle cluster together

## Get the gene ids in TCA cycle
id=id_allAttributes$ids[id_allAttributes$Pwy=="TCA"&!is.na(id_allAttributes$Pwy)] %>% unique

## Draw the dendrogram from the tree and color the wiltypes (Haven't found a way to remove the labels)
strainCC=1-abs(cor(t(All_Data),use="pairwise.complete.obs"))
class(strainCC); dim(strainCC)
hclust_pcc_complete=strainCC %>% as.dist %>% hclust

hclust_pcc_complete %>% as.dendrogram %>% set("by_labels_branches_col", value = id) %>% plot(raise.dendrogram(-100))


plot(hclust_pcc_complete %>% as.dendrogram %>% set("by_labels_branches_col", value = id))
#plot(raise.dendrogram(hclust_pcc_complete %>% as.dendrogram %>% set("by_labels_branches_col", value = id),-0.2))


#A way to raise (or move down) the branches. positive values raise, negative move down
#plot(dend15)
#plot(raise.dendrogram(dend15,-0.5))
