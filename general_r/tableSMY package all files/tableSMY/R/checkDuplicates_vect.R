checkDuplicates_vect <-
function(vect){
  if(sum(duplicated(vect))>=1){
    print("Some duplicates are found:")
    table(vect)                
  }else{return("Everything in this vector is unique")}
}
