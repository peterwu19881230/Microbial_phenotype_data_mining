#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("AnnotationDbi", version = "3.8")

##function to separate CC, MF BP: AnnotationDbi::Term()


##2 functions to creae subset of the obj GOannotation (GOannotation obj was created by strain1_strain2_allGO.R)

#========================================================================
filter_annot=function(annot1_annot2,branch){
  
  annot1=annot1_annot2[[1]]
  annot2=annot1_annot2[[2]]
  
  is_branch_1=(AnnotationDbi::Ontology(annot1)==branch)
  is_branch_2=(AnnotationDbi::Ontology(annot2)==branch)
  
  new_annot1=annot1[is_branch_1]; if(identical(new_annot1,character(0))) new_annot1=as.character(NA)  
  new_annot2=annot2[is_branch_2]; if(identical(new_annot2,character(0))) new_annot2=as.character(NA)  
  
  return(list(new_annot1,new_annot2))
}
#========================================================================

filter_GO_class=function(GOannotation,branch=""){
  
  if(branch=="") stop("No branch specified!")
  
  lapply(GOannotation,FUN=filter_annot,branch=branch)
}
#========================================================================

#========================================================================

filter_GO_class_parallel=function(GOannotation,branch=""){
  
  if(branch=="") stop("No branch specified!")
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  
  result=parLapply(cl=cl,GOannotation,fun=filter_annot,branch=branch)
  stopCluster(cl)
  return(result)
}
#========================================================================



#Test if the functions work
'
GOannotation_test=GOannotation[1:3]

filter_GO_class(GOannotation_test,branch="BP")
filter_GO_class(GOannotation_test,branch="CC")
filter_GO_class(GOannotation_test,branch="MF")

filter_GO_class(GOannotation_test,branch="nonsense")

BP=unlist(filter_GO_class(GOannotation_test,branch="BP"))
BP=BP[!is.na(BP)]
AnnotationDbi::Ontology(BP)

MF=unlist(filter_GO_class(GOannotation_test,branch="MF"))
MF=MF[!is.na(MF)]
AnnotationDbi::Ontology(MF)

CC=unlist(filter_GO_class(GOannotation_test,branch="CC"))
CC=CC[!is.na(CC)]
AnnotationDbi::Ontology(CC)


##Make sure the parallelized function works
identical(filter_GO_class(GOannotation_test,branch="BP"),filter_GO_class_parallel(GOannotation_test,branch="BP"))
identical(filter_GO_class(GOannotation_test,branch="MF"),filter_GO_class_parallel(GOannotation_test,branch="MF"))
identical(filter_GO_class(GOannotation_test,branch="CC"),filter_GO_class_parallel(GOannotation_test,branch="CC"))
'
#========================================================================


if(!file.exists("Data/sourced/GOannotation_separate.RData")){
  start.time=Sys.time()
  
  GOannotation_BP=filter_GO_class_parallel(GOannotation,branch="BP")
  
  end.time=Sys.time()
  end.time-start.time #Time difference of 2.521393 hours
  
  start.time=Sys.time()
  GOannotation_MF=filter_GO_class_parallel(GOannotation,branch="MF") 
  end.time=Sys.time()
  end.time-start.time
  
  start.time=Sys.time()
  GOannotation_CC=filter_GO_class_parallel(GOannotation,branch="CC") 
  end.time=Sys.time()
  end.time-start.time
  
  save(GOannotation_BP,GOannotation_MF,GOannotation_CC,file="Data/sourced/GOannotation_separate.RData")
}


#Get distance from MF, BP, CC,respectively
if(!file.exists("Data/Wang_pairwise_similarity_separate.RData")){
  
  ##parallelized code to get distances from Wang method 
  #BiocManager::install("AnnotationHub", version = "3.8")
  library(AnnotationHub)
  library(GOSemSim)
  
  ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)
  ECK_GO_MF=godata('org.EcK12.eg.db',ont="MF",computeIC = F)
  ECK_GO_CC=godata('org.EcK12.eg.db',ont="CC",computeIC = F)
  
  
  Wang_pairwise_similarity_separate=list()
  
  ##BP====================================================================================
  start.time=Sys.time()
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  
  
  #Export packages to each cluster (https://stackoverflow.com/questions/33761123/r-how-can-i-export-methods-provided-by-a-package-to-a-psock-cluster)
  clusterEvalQ(cl,c(library(AnnotationHub),library(GOSemSim),library(parallel))) # I tried to pass more than 1 library and it worked
  
  #Export required variables to each cluster (http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)
  clusterExport(cl=cl, c("GOannotation_BP","ECK_GO_BP")) # I tried to pass more than 1 variable and it worked
  
  Wang_pairwise_similarity_separate[["BP"]]=parLapply(cl=cl,X=GOannotation_BP,fun=function(annots){
    mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_BP, measure="Wang",combine="BMA")
  }
  ) 
  
  stopCluster(cl)
  
  
  end.time=Sys.time()
  end.time-start.time #Time difference of 8.071771 mins
  
  
  
  ##====================================================================================
  
  
  ##MF====================================================================================
  start.time=Sys.time()
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  
  
  #Export packages to each cluster (https://stackoverflow.com/questions/33761123/r-how-can-i-export-methods-provided-by-a-package-to-a-psock-cluster)
  clusterEvalQ(cl,c(library(AnnotationHub),library(GOSemSim),library(parallel))) # I tried to pass more than 1 library and it worked
  
  #Export required variables to each cluster (http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)
  clusterExport(cl=cl, c("GOannotation_MF","ECK_GO_MF")) # I tried to pass more than 1 variable and it worked
  
  Wang_pairwise_similarity_separate[["MF"]]=parLapply(cl=cl,X=GOannotation_MF,fun=function(annots){
    mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_MF, measure="Wang",combine="BMA")
  }
  ) 
  
  stopCluster(cl)
  
  
  end.time=Sys.time()
  end.time-start.time #Time difference of 5.200617 mins
  
  
  
  ##====================================================================================
  
  ##CC====================================================================================
  start.time=Sys.time()
  
  library(parallel)
  no_cores <- detectCores()-1
  cl <- makeCluster(no_cores)
  
  
  #Export packages to each cluster (https://stackoverflow.com/questions/33761123/r-how-can-i-export-methods-provided-by-a-package-to-a-psock-cluster)
  clusterEvalQ(cl,c(library(AnnotationHub),library(GOSemSim),library(parallel))) # I tried to pass more than 1 library and it worked
  
  #Export required variables to each cluster (http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)
  clusterExport(cl=cl, c("GOannotation_CC","ECK_GO_CC")) # I tried to pass more than 1 variable and it worked
  
  Wang_pairwise_similarity_separate[["CC"]]=parLapply(cl=cl,X=GOannotation_CC,fun=function(annots){
    mgoSim(annots[[1]],annots[[2]], semData=ECK_GO_CC, measure="Wang",combine="BMA")
  }
  ) 
  
  stopCluster(cl)
  
  
  end.time=Sys.time()
  end.time-start.time #Time difference of 13.63786 mins
  
  
  
  ##====================================================================================
  
  
  save(Wang_pairwise_similarity_separate,file="Data/sourced/Wang_pairwise_similarity_separate.RData")
  
  
}






