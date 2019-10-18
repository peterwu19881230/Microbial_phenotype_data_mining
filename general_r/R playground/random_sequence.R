leng=1000

generate_seq=function(leng=5,seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  paste( sample(c("A","C","G","T"),leng,replace = T),collapse="")
}

#test the function
generate_seq()
generate_seq(seed=101)
generate_seq(leng=100,seed=101)
