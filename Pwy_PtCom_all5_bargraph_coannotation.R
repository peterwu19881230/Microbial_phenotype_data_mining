## Data obj preparation
#------------------------------------------------------------------------------------------------------------
noOfGenes=function(df,annotLogical){
  df[annotLogical,c("strain1","strain2")] %>% as.matrix %>% as.vector %>% unique %>% length
}

df=strain1strain2_allAnnotations_allDistances

Pwy=noOfGenes(df,as.logical(df$Pwy))
pcomplex=noOfGenes(df,as.logical(df$pcomplex))
operon=noOfGenes(df,as.logical(df$operon))
regulator=noOfGenes(df,as.logical(df$regulator))
kegg_modules=noOfGenes(df,as.logical(df$kegg_modules))

PwyANDPcomplex=noOfGenes(df,as.logical(df$Pwy & df$pcomplex))
all5=noOfGenes(df,as.logical(df$Pwy & df$pcomplex & df$operon & df$regulator & df$kegg_modules))


sets=c(Pwy=Pwy,pcomplex=pcomplex,operon=operon,regulator=regulator,kegg_modules=kegg_modules,PwyANDPcomplex=PwyANDPcomplex,all5=all5)

names(sets)=c("Pathway","Protein complex","Operon","Regulon","KEGG modules","Pathway & Protein complex","All annotation sets")
#------------------------------------------------------------------------------------------------------------





#A bargraph containing no. of genes involved in co-annotated pairs
#================================================================================================================
name=names(sets)

inputDF=data.frame(name=factor(name,levels=name),
                   count=sets
)

Cols=c("#56B4E9","#F3518A","#24D7B7","#267A6B","#3533A7") 


p1=ggplot(data=inputDF[1:5,],aes(x=name,y=count,fill=name))+geom_bar(stat="identity")+
  theme_minimal()+
  theme(text=element_text(size=20),axis.text.x=element_blank())+
  labs(x="",y="No. of co-annotated genes",aesthetic='custom text',fill="Annotation set")+  
  scale_fill_manual(values=Cols)+
  theme(legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylim(0, 4000)


Cols=c("#BE70EA","#09D38A") 

p2=ggplot(data=inputDF[6:7,],aes(x=name,y=count,fill=name))+geom_bar(stat="identity")+
  theme_minimal()+
  theme(text=element_text(size=20),axis.text.x=element_blank())+
  labs(x="",y="",aesthetic='custom text',fill="Annotation set")+  #y label is specified in the above plot
  scale_fill_manual(values=Cols)+
  theme(legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylim(0, 4000)



grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")) #ref: https://gist.github.com/tomhopper/faa24797bb44addeba79

#================================================================================================================




#I can also do a bargraph containing no. of "genes pairs" involved in co-annotated pairs












