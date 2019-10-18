#compute Fibonacci numbers

#Note: n cannot be 0 or negative
Feb=function(n){
  Febs=numeric(n)
  Febs[1]=0
  Febs[2]=1
  
  if(n==1){
    return(0)
  }else if(n==2){
    return(1)
  }else{
    for(i in 3:n){
      Febs[i]=Febs[i-1]+Febs[i-2]
    }
    
    return(Febs[n])
  }
  
}


