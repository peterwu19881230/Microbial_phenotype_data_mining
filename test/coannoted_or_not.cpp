#include <Rcpp.h>
using namespace Rcpp;
#include <omp.h> 
//On Mac I have to do this first: brew install libomp. Ref: https://stackoverflow.com/questions/25990296/how-to-include-omp-h-in-os-x
//And then I might have to do this (haven't tried): http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/

// [[Rcpp::export]]
IntegerVector coannoted_or_not(DataFrame distance_table ,List attribute_list) { 
  
  CharacterVector firstCol = distance_table[0];
  CharacterVector secondCol = distance_table[1];
  
  
  int dim1 = firstCol.size(); //Use the first column to get the dimension 1 of the distance_table
  
  IntegerVector result=rep(0,dim1); //Ref: http://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
  for(int i = 0; i <= (dim1-1); i++){
    
    //ref 1: https://stackoverflow.com/questions/23282849/rcpp-charactervector-why-arent-individual-elements-stdstrings
    //ref 2: https://www.geeksforgeeks.org/stdstring-class-in-c/
    std::string id1= as<std::string>(firstCol[i]); 
    std::string id2= as<std::string>(secondCol[i]);
    
    CharacterVector attr1=attribute_list[id1];
    CharacterVector attr2=attribute_list[id2];
    
    //Note: comparing c("A","B","C") with NA gives character(0), but NA with NA returns NA (logical NA)
    //ref: https://stackoverflow.com/questions/40515817/using-intersect-in-r-with-rcpp
    CharacterVector intersection=intersect(attr1,attr2);
    
    //ref: http://gallery.rcpp.org/articles/working-with-missing-values/
      if(intersection.size()>=1 && any( is_na(intersection) )==false){ //ref: https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-sugar.pdf
          result[i]=1;
      }
        
        
  }
  return result;
}



//Parallel version of the above function
//Ref: Udemy Advanced R: parallel-rcpp
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
IntegerVector coannoted_or_not_parallel(DataFrame distance_table ,List attribute_list,int threads) { 
  
  CharacterVector firstCol = distance_table[0];
  CharacterVector secondCol = distance_table[1];
  
  
  int dim1 = firstCol.size(); //Use the first column to get the dimension 1 of the distance_table
  
  IntegerVector result=rep(0,dim1); //Ref: http://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
  
  int nThreads = 0; 
  omp_set_num_threads(threads);
{
    nThreads = omp_get_num_threads();
  
  for(int i = 0; i <= (dim1-1); i++){
    
    //ref 1: https://stackoverflow.com/questions/23282849/rcpp-charactervector-why-arent-individual-elements-stdstrings
    //ref 2: https://www.geeksforgeeks.org/stdstring-class-in-c/
    std::string id1= as<std::string>(firstCol[i]); 
    std::string id2= as<std::string>(secondCol[i]);
    
    CharacterVector attr1=attribute_list[id1];
    CharacterVector attr2=attribute_list[id2]; 
    
    //Note: comparing c("A","B","C") with NA gives character(0), but NA with NA returns NA (logical NA)
    //ref: https://stackoverflow.com/questions/40515817/using-intersect-in-r-with-rcpp
    CharacterVector intersection=intersect(attr1,attr2);
    
    //ref: http://gallery.rcpp.org/articles/working-with-missing-values/
    if(intersection.size()>=1 && any( is_na(intersection) )==false){ //ref: https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-sugar.pdf
      result[i]=1;
    }
    
    
  }
  
}
  return result;
}
