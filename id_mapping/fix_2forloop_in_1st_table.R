indices<-c()
###2 for() is too slow. Maybe use sapply? => Changed the first for() into sapply. No significant difference

for(j in 1:3979){
  
  
  
  for(i in 1:4494){
    IsFound<-sum(geneNames_genes.col[i,1:5] %in% id_ECKs_CorrectedECKs_AssociatedGeneNames$associated_gene_names[j])
    
    if(IsFound==1){
      indices<-c(indices,i) 
      break
    }
  }
  if(is.na(indices[j])==T){
    indices<-c(indices,"Not found")
  }

  
}


indices<-list()
###2 for() is too slow. Maybe use sapply? => Changed the first for() into sapply. No significant difference

##Should only serach for the name column, not all the synonyms columns
##EcoCYC smart table
##parse genes.dat instead of genes.col
for(j in 1:3979){
  index<-grep(id_ECKs_CorrectedECKs_AssociatedGeneNames$associated_gene_names[j],all_gene_synonyms)
  
  if(length(index)>=2){  ## this is to deal situations like: "pps" find the indices of 1479: ppsR ydiA & 1717: ppsA pps
    for(i in 1:length(index)){
      isFound<-grep(paste("^",id_ECKs_CorrectedECKs_AssociatedGeneNames$associated_gene_names[j],"$"),
                    strsplit(all_gene_synonyms[i],split=",") %>% unlist )
      if(length(isFound)==1){
        ##index<-"???"
        index<-index[i]
      }else{ } 
    }
    
  }
    
  indices[[j]]<-index ##double brackets help. Otherwise there will be: Warning message: number of items to replace is not a multiple of replacement length
}


