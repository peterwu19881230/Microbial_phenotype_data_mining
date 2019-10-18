#Goal:Load commonly used stuff for the project of Nichols'

#Set stringsAsFactors = FALSE so things won't be automatically converted to factor
options(stringsAsFactors = FALSE) 


#Load commonly used packages. If not installed, install.
#========================================================================
if(!require(pacman)){
  install.packages("pacman")
  library(pacman)
}
#Tidyverse: https://www.tidyverse.org/packages/
#Note: BiocManager has to be installed before its downstream packages
pacman::p_load(tidyverse,BiocManager,readxl,factoextra,pheatmap,ComplexHeatmap,XML,reshape2,googlesheets,googledrive,
               circlize,ape,dendextend,sqldf,RMySQL,AnnotationHub,infotheo,foreach,doParallel,parallel,
               stringi,purrr,ggfortify) #Note ggplot2 is loaded when tidyverse or factoextra is loaded 

if(!require(GOSemSim)){
  BiocManager::install("GOSemSim")
}

if(!require(org.EcK12.eg.db)){
  BiocManager::install("org.EcK12.eg.db")
}

#========================================================================



#Load data and other useful files:
##source("UnitTesting/Datasets.R") ##optional



##Source all the self-defined functions
functions=dir("functions/")
paths=paste("functions/",functions,sep="")
for(func in paths){
  source(func)
}

functions=dir("general_r/functions/")
paths=paste("general_r/functions/",functions,sep="")
for(func in paths){
  source(func)
}

##Load useful data in "Data" directory


#Lazy loading ref: https://www.r-bloggers.com/lazy-load-with-archivist/
Data=dir("Data/sourced")
paths=paste("Data/sourced/",Data[!( grepl(".rdx",Data,fixed = T) | grepl(".rdb",Data,fixed = T)) ],sep="")
for(path in paths){
  if(class(try(lazyLoad(path)))=="try-error"){
    lazyLoad = local({load(path); 
      environment()})
    tools:::makeLazyLoadDB(lazyLoad, path)
    lazyLoad(path)
  }else lazyLoad(path)
  
  print(paste(path," has been loaded"))
}




##Source other useful files
load("id_mapping/genes.dat.RData")






