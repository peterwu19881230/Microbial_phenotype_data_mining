## A very fast way of computing hamming distance matrix.Modified from: https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/
## I verified that it works. However it doesn't deal with data with NAs. => I wasn't able to modify it to deals with NA.
##(!)In the original function it calculates the distance between columns. I changed it to calculates the rows
hamming_distance <- function(X) {
  uniqs <- unique(as.vector(X))
  U <- X == uniqs[1]
  H <- U %*% t(U)
  for ( uniq in uniqs[-1] ) {
    U <- X == uniq
    H <- H + U %*% t(U)
  }
  ncol(X) - H
}

#This function assumes that there are only 1,0,-1 in the input matrix
##(!)1 caveat of this is that the diagnal distance (self-to-self relationship) is not always the same, which is fine
##because 0 are treated more like an unknown value
modified_hamming_distance=function(mat){
  X=list(ifelse(mat==-1,1,0),
         ifelse(mat==0,1,0),
         ifelse(mat==1,1,0)
  )
  
  S=list() #S stands for similarity
  count=0
  for(i in 1:3){  #Index: 1->-1,  2->0, 3->1
    for(j in 1:3){
      count=count+1
      S[[count]]=X[[i]]%*%t(X[[j]])
    }
  }
  
  
  sum_S=S[[1]]*1+S[[9]]*1 ##(-1,-1) (1,1) +1
  ncol(mat) - sum_S #Distance= total number - similarity 
}


#This function is the same as the previous one except it gives punishment
#This function assumes that there are only 1,0,-1 in the input matrix
##(!)1 caveat of this is that the diagnal distance (self-to-self relationship) is not always the same, which is fine
##because 0 are treated more like an unknown value
modified_hamming_distance_2=function(mat){
  X=list(ifelse(mat==-1,1,0),
         ifelse(mat==0,1,0),
         ifelse(mat==1,1,0)
  )
  
  S=list() #S stands for similarity
  count=0
  for(i in 1:3){  #Index: 1->-1,  2->0, 3->1
    for(j in 1:3){
      count=count+1
      S[[count]]=X[[i]]%*%t(X[[j]])
    }
  }
  
  #4 types of relationships I want to deal with:
  sum_S=S[[1]]*1+S[[9]]*1+ ##(-1,-1) (1,1) +1
    S[[5]]*0+ ##(0,0) +0
    S[[3]]*(-1)+S[[7]]*(-1)+ ##(-1,1) (1,-1) -1
    S[[4]]*(-0.5)+S[[2]]*(-0.5)+S[[6]]*(-0.5)+S[[8]]*(-0.5) ##(0,-1) (-1,0) (0,1) (1,0) -0.5
  
  ncol(mat) - sum_S #Distance= total number - similarity 
}

#This function is the same as modified_hamming_distance_2() except it doesn't give that -0.5 punishment
#This function assumes that there are only 1,0,-1 in the input matrix
##(!)1 caveat of this is that the diagnal distance (self-to-self relationship) is not always the same, which is fine
##because 0 are treated more like an unknown value
modified_hamming_distance_3=function(mat){
  X=list(ifelse(mat==-1,1,0),
         ifelse(mat==0,1,0),
         ifelse(mat==1,1,0)
  )
  
  S=list() #S stands for similarity
  count=0
  for(i in 1:3){  #Index: 1->-1,  2->0, 3->1
    for(j in 1:3){
      count=count+1
      S[[count]]=X[[i]]%*%t(X[[j]])
    }
  }
  
  #3 types of relationships I want to deal with:
  sum_S=S[[1]]*1+S[[9]]*1+ ##(-1,-1) (1,1) +1
    S[[5]]*0+ ##(0,0) +0
    S[[3]]*(-1)+S[[7]]*(-1) ##(-1,1) (1,-1) -1
  
  ncol(mat) - sum_S #Distance= total number - similarity 
}



#The problem of this version is that total No. of conditions are not used (various available 0s are not used in each pair)
#strains that are all 0 are removed because this function would give the distance=NaN
modified_hamming_distance_4=function(mat){
  X=list(ifelse(mat==-1,1,0),
         ifelse(mat==0,1,0),
         ifelse(mat==1,1,0)
  )
  
  S=list() #S stands for similarity
  count=0
  for(i in 1:3){  #Index: 1->-1,  2->0, 3->1
    for(j in 1:3){
      count=count+1
      S[[count]]=X[[i]]%*%t(X[[j]])
    }
  }
  
  sum_S=S[[1]]*1+S[[9]]*1+ ##(-1,-1) (1,1) +1
    S[[5]]*0+ ##(0,0) +0
    S[[3]]*(-1)+S[[7]]*(-1) ##(-1,1) (1,-1) -1
  
  totalCondUsed=(S[[1]]+S[[9]]+S[[3]]+S[[7]]) #Total conditions used to calculate scores for all pairwise combinations
  
  1-abs(sum_S/totalCondUsed) #Distance=1 - similarity /total number of condition used for calculating similarity
}


#The problem of this version is that total No. of conditions are not used (various available 0s are not used in each pair)
#strains that are all 0 are removed because this function would give the distance=NaN
##The same as modified_hamming_distance_4 but:   
##1. (-1,1) is not punished 
##2. (1,1),(-1,-1) have 1 point, (-1,1) has 0.5 point (Debby pointed out that even if it's different direction, they are somehow correlated <- biological process)
##3. (0,1) pairs are included in the total condition used 
modified_hamming_distance_5=function(mat){
  X=list(ifelse(mat==-1,1,0),
         ifelse(mat==0,1,0),
         ifelse(mat==1,1,0)
  )
  
  S=list() #S stands for similarity
  count=0
  for(i in 1:3){  #Index: 1->-1,  2->0, 3->1
    for(j in 1:3){
      count=count+1
      S[[count]]=X[[i]]%*%t(X[[j]])
    }
  }
  
  sum_S=S[[1]]*1+S[[9]]*1+ ##(-1,-1) (1,1) +1
    S[[5]]*0+ ##(0,0) +0
    S[[3]]*0.5+S[[7]]*0.5 ##(-1,1) (1,-1) +0.5
  
  totalCondUsed=(S[[1]]+S[[2]]+S[[3]]+S[[4]]+S[[6]]+S[[7]]+S[[8]]+S[[9]]) #Total conditions used to calculate scores for all pairwise combinations
  
  1-sum_S/totalCondUsed #Distance=1 - similarity /total number of condition used for calculating similarity
}


##Test script
##modified_hamming_distance_4(Ternary_Data[1:100,1:100])
##temp=modified_hamming_distance_4(Ternary_Data)
##melted=meltANDsort_dist(temp)
##melted2=meltANDsort_dist(MHD3_ternary)

##Test script
##set.seed(101)
##mat=matrix(sample(c(-1,0,1),12,replace=T),ncol=3)
##modified_hamming_distance_5(mat)
##mat=matrix(c(0,1,0,1,-1,1,0,1,0,-1,0,1,-1,0,1,1),nrow=2,byrow = T)
##modified_hamming_distance_5(mat)

