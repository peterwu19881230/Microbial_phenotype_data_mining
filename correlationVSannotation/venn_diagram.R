#Ref: https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html

# install.packages('venneuler')
library(venneuler)

#An example
#vd <- venneuler(c(A=30, B=30, C=110, "A&B"=10, "A&C"=20, "B&C"=10 ,"A&B&C"=10))
#plot(vd)

names(id_allAttributes)


#Here I haven't accounted strains that have 2 or more annotations from any annotation set

getAnnotatedID=function(attr){
  id_attr=id_allAttributes[,c("ids",attr)] %>% unique
  annotatedID=id_attr$ids[!is.na(id_attr[[paste(attr)]])] %>% unique
}

allID=as.character(1:3979)
Pwy=getAnnotatedID("Pwy")
kegg_modules=getAnnotatedID("kegg_modules")
pcomplex=getAnnotatedID("pcomplex")
regulator=getAnnotatedID("regulator")
operon=getAnnotatedID("operon")

sets=list(allID=allID,Pwy=Pwy,kegg_modules=kegg_modules,pcomplex=pcomplex,regulator=regulator,operon=operon)

names(sets)=c("All 3979 strains","EcoCyc Pathway","KEGG Modules","EcoCyc Protein complex","Regulon","Operon")



get_combinations=function(sets){ #I have verified that this function works with sets variable that has length from 2~5
  
  length_forVen=c()
  for(i in 1:length(sets)){
    
    if(i==1){
      length_forVen=c(length_forVen,sapply(sets,FUN=length))
      names(length_forVen)=names(sets)
    }
    
    if(i==2){
      
      comb=combn(length(sets),i)
      
      for(j in 1:dim(comb)[2]){
        
        index1=comb[1,j]
        index2=comb[2,j]
        
        intersectedID=intersect(sets[[index1]],sets[[index2]]) %>% length 
        
        length_forVen=c(length_forVen,intersectedID)
        names(length_forVen)[length(length_forVen)]=paste(names(sets[index1]),names(sets[index2]),sep="&")
      }
      
    }
    
    if(i==3){
      comb=combn(length(sets),i)
      
      for(j in 1:dim(comb)[2]){
        
        index1=comb[1,j]
        index2=comb[2,j]
        index3=comb[3,j]
        
        intersectedID=intersect(sets[[index1]],sets[[index2]]) 
        intersectedID=intersect(intersectedID,sets[[index3]]) %>% length 
        
        length_forVen=c(length_forVen,intersectedID)
        names(length_forVen)[length(length_forVen)]=paste(names(sets[index1]),names(sets[index2]),names(sets[index3]),sep="&")
      }
    }
    
    if(i==4){
      comb=combn(length(sets),i)
      
      for(j in 1:dim(comb)[2]){
        
        index1=comb[1,j]
        index2=comb[2,j]
        index3=comb[3,j]
        index4=comb[4,j]
        
        intersectedID=intersect(sets[[index1]],sets[[index2]]) 
        intersectedID=intersect(intersectedID,sets[[index3]]) 
        intersectedID=intersect(intersectedID,sets[[index4]]) %>% length
        
        length_forVen=c(length_forVen,intersectedID)
        names(length_forVen)[length(length_forVen)]=paste(names(sets[index1]),names(sets[index2]),names(sets[index3]),names(sets[index4]),sep="&")
      }
      
    } 
    
    if(i==5){
      comb=combn(length(sets),i)
      
      for(j in 1:dim(comb)[2]){
        
        index1=comb[1,j]
        index2=comb[2,j]
        index3=comb[3,j]
        index4=comb[4,j]
        index5=comb[5,j]
        
        intersectedID=intersect(sets[[index1]],sets[[index2]]) 
        intersectedID=intersect(intersectedID,sets[[index3]])
        intersectedID=intersect(intersectedID,sets[[index4]])
        intersectedID=intersect(intersectedID,sets[[index5]]) %>% length
        
        length_forVen=c(length_forVen,intersectedID)
        names(length_forVen)[length(length_forVen)]=paste(names(sets[index1]),names(sets[index2]),names(sets[index3]),names(sets[index4]),names(sets[index5]),sep="&")
      }
      
    } 
    
    if(i==6){
      intersectedID=intersect(sets[[1]],sets[[2]]) 
      intersectedID=intersect(intersectedID,sets[[3]]) 
      intersectedID=intersect(intersectedID,sets[[4]]) 
      intersectedID=intersect(intersectedID,sets[[5]]) 
      intersectedID=intersect(intersectedID,sets[[6]]) %>% length
      
      length_forVen=c(length_forVen,intersectedID)
      names(length_forVen)[length(length_forVen)]=paste(names(sets[1]),names(sets[2]),names(sets[3]),names(sets[4]),names(sets[5]),names(sets[6]),sep="&")
    }
      
  }
  
  return(length_forVen)
  }
  
  
  
vec_forVen=get_combinations(sets[1:6]) 
#including the "All 3979 strains" will have this problem: Other annotation sets are not plotted within "All 3979 strains"
#There are some more problems. Eg. Regulon should be within operon


vd = venneuler(vec_forVen)
#By googling + looking into the structure of vd, I wansn't able to change the color of the labels

plot(vd) 

#pdf(file="venDiagram_allAnnotationSets")
#plot(vd)
#dev.off()

#What's the result of accounting strains that have 2 or more annotations from every annotation set?








