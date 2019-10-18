#2 gene pathway analysis


#Find gene pairs in pathways where there are only 2 genes involved (including 2 gene pathways + pathways that only have 2 genes in Nichols')
df=id_allAttributes[,c("ids","Pwy")] %>% unique

tab=table(df$Pwy)
twoGenePwy=names(tab[tab==2])



ids=df$ids[df$Pwy %in% twoGenePwy] %>% unique 
## Whether using unique() or not makes a difference
## Supposedly there are genes that are involved in multiple 2-gene pathways or involved in mutiple pathways where only 2 
## of Nichols' genes are involved







df2=strain1strain2_allAnnotations_allDistances
pcc=(1-df2$pcc)[( df2$strain1 %in% as.numeric(ids) ) & 
            ( df2$strain2 %in% as.numeric(ids) ) &
              df2$Pwy==1]

hist(pcc,main="|PCC| in 2 gene pathways (or pathways where only 2 genes in Nichols' are involved)",xlab="|PCC|")
text(x=0.6,y=60,labels=paste(names(summary(pcc)),collapse=" - "))
text(x=0.6,y=50,labels=paste(round(summary(pcc),3),collapse=" - "))






