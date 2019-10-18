## Data obj preparation
#------------------------------------------------------------------------------------------------------------
names(id_allAttributes)

#Here I haven't accounted strains that have 2 or more annotations from any annotation set

getAnnotatedID=function(attr){
  id_attr=id_allAttributes[,c("ids",attr)] %>% unique
  annotatedID=id_attr$ids[!is.na(id_attr[[paste(attr)]])] %>% unique
}

allID=as.character(1:3979)
Pwy=getAnnotatedID("Pwy")
kegg_modules=getAnnotatedID("kegg_modules")

##For no. of complex. I don't want to show the monomer annotations because they are by definition not complexes 
##(for homo-oligomer I think it's ok to let them stay there because 1. they are by definition complexes 2. they don't cause problem in determining co-annotated pairs)
attr="pcomplex"
isMonomer=grepl("MONOMER",id_allAttributes$pcomplex)
id_attr=id_allAttributes[!isMonomer,c("ids",attr)] %>% unique
pcomplex=id_attr$ids[!is.na(id_attr[[paste(attr)]])] %>% unique

regulator=getAnnotatedID("regulator")
operon=getAnnotatedID("operon")


sets=list(allID=allID,Pwy=Pwy,pcomplex=pcomplex,operon=operon,regulator=regulator,kegg_modules=kegg_modules)

names(sets)=c("All 3979 strains","Pathway","Protein complex","Operon","Regulon","Modules")
str(sets)
#------------------------------------------------------------------------------------------------------------


#graph things based on individual no. of annotated ids in each annotation
#================================================================================================================
name=names(sets)

inputDF=data.frame(name=factor(name,levels=name),
                   count=sapply(sets,length)
)

Cols=c("grey","#56B4E9","#F3518A","#24D7B7","#267A6B","#3533A7") 

p1=ggplot(data=inputDF,aes(x=name,y=count,fill=name))+geom_bar(stat="identity")+
  theme_minimal()+
  theme(text=element_text(size=20),axis.text.x=element_blank())+
  labs(x="",y="No. of annotated genes",aesthetic='custom text',fill="Annotation set")+  
  scale_fill_manual(values=Cols)+
  theme(legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylim(0, 4000)
  

#================================================================================================================


#graph things based on the way I use combinations of annotations
#================================================================================================================
PwyANDPcomplex=intersect(Pwy,pcomplex)
str(PwyANDPcomplex)

all5=Pwy %>% intersect(kegg_modules) %>% 
  intersect(pcomplex) %>%
  intersect(regulator) %>%
  intersect(operon) 
str(all5)

name=c("All 3979 strains",
       "Pathway",
       "Protein complex",
       "Pathway & Protein complex",
       "All annotation sets")

inputDF=data.frame(name=factor(name,levels=name),
                   count=c(length(allID),length(Pwy),length(pcomplex),length(PwyANDPcomplex),length(all5))
                   )


Cols=c("grey","#56B4E9","#F3518A","#BE70EA","#09D38A") 
 

p2=ggplot(data=inputDF,aes(x=name,y=count,fill=name))+geom_bar(stat="identity")+
  theme_minimal()+
  theme(text=element_text(size=20),axis.text.x=element_blank())+
  labs(x="",y="No. of annotated genes",aesthetic='custom text',fill="Annotation set(s)")+  
  scale_fill_manual(values=Cols)+
  theme(legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylim(0, 4000)
#================================================================================================================

#graph only combinations of annotation sets
#================================================================================================================
PwyANDPcomplex=intersect(Pwy,pcomplex)
str(PwyANDPcomplex)

all5=Pwy %>% intersect(kegg_modules) %>% 
  intersect(pcomplex) %>%
  intersect(regulator) %>%
  intersect(operon) 
str(all5)

name=c("Pathway & Protein complex",
       "All annotation sets")

inputDF=data.frame(name=factor(name,levels=name),
                   count=c(length(PwyANDPcomplex),length(all5))
)


Cols=c("#BE70EA","#09D38A") 


p3=ggplot(data=inputDF,aes(x=name,y=count,fill=name))+geom_bar(stat="identity")+
  theme_minimal()+
  theme(text=element_text(size=20),axis.text.x=element_blank())+
  labs(x="",y="No. of annotated genes",aesthetic='custom text',fill="Annotation sets")+  
  scale_fill_manual(values=Cols)+
  theme(legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ylim(0, 4000)
#================================================================================================================

grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p3), size = "last")) #ref: https://gist.github.com/tomhopper/faa24797bb44addeba79


