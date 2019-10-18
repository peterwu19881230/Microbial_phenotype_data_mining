#Construct strain1strain2_allAnnotations on ADA (copy paste the following code to ADA)

#dependency: id_allAttributes
##scp /Users/peterwu/TAMU\ Google\ Drive/OMP/NicholsDataMining_Data/Data/sourced/id_allAttributes.RData peterwu19881230@ada.tamu.edu:/home/peterwu19881230/

##ssh peterwu19881230@ada.tamu.edu
##module load R/3.5.0-iomkl-2017b-recommended-mt
##R


#1st and 2nd column for strain pairs
strain1strain2_allAnnotations=as.data.frame(t(combn(1:3979,2)))
names(strain1strain2_allAnnotations)=c("strain1","strain2")


#Add columns for various annotations: Pwy, pcomplex, regulator, operon, kegg_modules, GO
library(parallel)
no_cores <- detectCores()-1 
cl <- makeCluster(no_cores)


attrs=c("Pwy", "pcomplex", "regulator", "operon", "kegg_modules")





for(i in 1:length(attrs)){
  
  start.time=Sys.time()
  
  attr=attrs[i]
  clusterExport(cl,c("id_allAttributes","attr"))
  
  result=parApply(cl=cl,X=strain1strain2_allAnnotations,MARGIN=1,FUN=function(strain1strain2){ 
    
    id1=strain1strain2[1]
    id2=strain1strain2[2]
    
    annot1=unique(id_allAttributes[id_allAttributes$ids==id1,attr])
    annot2=unique(id_allAttributes[id_allAttributes$ids==id2,attr])
    print(annot1)
    print(annot2)    
    ifelse( sum( annot1 %in% annot2)>=1 & 
              !(anyNA(annot1)) & 
              !(anyNA(annot2)),
            1,0)
    
  })
  
  end.time=Sys.time()
  end.time-start.time
  
  save(result,file=paste("strain1strain2_annotCol_",attr,".RData",sep=""))
}

#scp peterwu19881230@ada.tamu.edu:/home/peterwu19881230/strain1strain2_annotCol_* Data/

load("Data/strain1strain2_annotCol_Pwy.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,Pwy=result)
load("Data/strain1strain2_annotCol_pcomplex.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,pcomplex=result)
load("Data/strain1strain2_annotCol_operon.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,operon=result)
load("Data/strain1strain2_annotCol_regulator.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,regulator=result)
load("Data/strain1strain2_annotCol_kegg_modules.RData")
strain1strain2_allAnnotations=cbind(strain1strain2_allAnnotations,kegg_modules=result)

##Load GO Data....



#type check
class(strain1strain2_allAnnotations)
apply(strain1strain2_allAnnotations,2,class)

#save(strain1strain2_allAnnotations,file="Data/strain1strain2_allAnnotations.RData")