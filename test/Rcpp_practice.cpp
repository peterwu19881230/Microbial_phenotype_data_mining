#include <Rcpp.h>
using namespace Rcpp;




//Ref: https://stackoverflow.com/questions/40515817/using-intersect-in-r-with-rcpp
//Note: I changed the IntegerVector to CharacterVector and it works
//Note: comparing c("A","B","C") with NA gives character(0), but NA with NA returns NA (logical NA)
// [[Rcpp::export]]
CharacterVector intersect_char( CharacterVector x, CharacterVector y){
  
  CharacterVector result=intersect( x, y );
  
  return result;
}



//From this function I know that the length of character(0) is 0, length of NA is 1
// [[Rcpp::export]]
IntegerVector length_cpp( CharacterVector x, IntegerVector y){
  
  int length_x=x.size();
  int length_y=y.size();
  
  IntegerVector result=IntegerVector::create(length_x,length_y);
  
  return result;
}


// [[Rcpp::export]]
LogicalVector is_naC2(NumericVector x) {
  
  //this is the sugar version of is_na() (ref: http://adv-r.had.co.nz/Rcpp.html#rcpp-na)
  return is_na(x); //is_na() can detect logical NA
}



//This function tells me: comparing c("a","b","c") with NA is different from comparing NA with NA
// [[Rcpp::export]]
bool length_x(CharacterVector x, CharacterVector y) {
  CharacterVector intersection=intersect(x,y);
  int leng=intersection.size();
  bool result=(leng>=1);
  return result;
}
//test code in R:
//length_x(c("a","b"),NA)
//length_x(c("a","b","c"),c("a","b","c"))



//This seems to be a good way to deal with NA
//Ref: http://dirk.eddelbuettel.com/code/rcpp/Rcpp-sugar.pdf
// [[Rcpp::export]]
bool my_anyNA(CharacterVector x) {
  bool result=any( is_na(x) );
  return result;
}


//Test when 2 char vectors are the inputs
// [[Rcpp::export]]
bool length_y(CharacterVector x, CharacterVector y) {
  CharacterVector intersection=intersect(x,y);
  
  bool result=any( is_na(intersection) );
  return result;
}


//Test if rep() function from sugar does the samething like in R: => yes it does
//Ref: http://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
// [[Rcpp::export]]
NumericVector rep_test(double x,int n){
  NumericVector result = rep( x, n );
  return result;
}






