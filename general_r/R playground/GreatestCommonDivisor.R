GreatestCommonDivisor=function(a,b){
  
  small=min(a,b)
  big=max(a,b)
  
  
  while(small!=0){
    
    new_small=big %% small
    new_big=small
    
    small=new_small
    big=new_big
    
  }

  return(big)
}
