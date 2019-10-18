

## For positive fitness
#########
no_of_genePairs_positive=choose(colSums(positive_phenotype_matrix),2)
co_annot_table=coannotation_positive 

result=list()
for(annot in c("Pwy","pcomplex","operon","regulator","kegg_modules")){
  
  
  result[[annot]]$avg_no_of_co_annotation=co_annot_table[annot,]
  names(result[[annot]]$avg_no_of_co_annotation)=colnames(All_Data_NAimputed)
  
  result[[annot]]$enrichment=co_annot_table[annot,]/random_coannotation[annot]
  names(result[[annot]]$enrichment)=colnames(All_Data_NAimputed)
  
  
    
  p_vals=c()
  for(i in 1:324){
    
    
    if(no_of_genes_positive[i] %in% c(0,1)){
      p_vals=c(p_vals,NA)
    }else{
      
      no_of_co_annotations=co_annot_table[annot,i]*no_of_genePairs_positive[i] #no. of co-annotations = avg co-annotations * no. of strains
      no_of_no_co_annotations=choose(no_of_genes_positive[i],2)-no_of_co_annotations
      population_co_annotations_exceptSubset=sum(strain1strain2_allAnnotations_allDistances[[annot]])-no_of_co_annotations
      population_no_co_annotations_exceptSubset=7914231-sum(strain1strain2_allAnnotations_allDistances[[annot]])-population_co_annotations_exceptSubset
      
      
      
      p_val=fisher.test(matrix(c(
        no_of_co_annotations,
        no_of_no_co_annotations,
        
        population_co_annotations_exceptSubset,
        population_no_co_annotations_exceptSubset),
        ,ncol=2,byrow=F),
        ,alternative = "greater")$p.value
    
      
      p_vals=c(p_vals,p_val)
      
      } 
    
  }
  
  result[[annot]]$p_vals=p_vals
  
}
#########