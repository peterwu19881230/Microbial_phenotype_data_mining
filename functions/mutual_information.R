#Mutual information
#mutualInfo_continouous=function(){}
##(Mutual information for the quantitative data. Not sure this is important for the paper)



##self-dfined mutual info for binary data. Ref: https://www.researchgate.net/post/How_can_i_calculate_Mutual_Information_theory_from_a_simple_dataset
##An article about a solution to when p(x,y) is 0: https://stats.stackexchange.com/questions/73502/conditional-mutual-information-and-how-to-deal-with-zero-probabilities
##(!) I think p(x,y) calculated in this ref were wrong
## Input: X, Y are paired data
mutualInfo_binary_old=function(X,Y,base=2){
  
  n_x=length(unique(X))
  n_y=length(unique(Y))
  x_bin=unique(X)
  y_bin=unique(Y)
  p_x=numeric(n_x)
  p_y=numeric(n_y)
  p_xy=numeric(n_x*n_y) 
  mutual_info_each=numeric(n_x*n_y)
  
  count=0
  for(i in 1:n_x){
    p_x[i]=sum(x_bin[i]==X)/length(X)
    
    for(j in 1:n_y){
      count=count+1
      p_y[j]=sum(y_bin[j]==Y)/length(Y)
      
      p_xy[count]=sum( (x_bin[i]==X) & (y_bin[j]==Y)  )/length(X) #Doesn't matter whether length(X) / length(Y) are used
      
      if( p_xy[count]==0 ){
        mutual_info_each[count]=0 #This is to prevent NaN generation when p(x,y) = 0 ( log(0)=-Inf )
      }else{mutual_info_each[count]=p_xy[count]*log(p_xy[count]/(p_x[i]*p_y[j]),base=base)} 
      
    }
  }
  mutual_info=sum(mutual_info_each)
  return(mutual_info)
}


#A shorter version of the above one. It gives the correct results but I haven't tested about the speed
mutualInfo_binary=function(X,Y,base=2){
pTabX=table(X)/length(X)
pTabY=table(Y)/length(Y)
p_xy=numeric(length(names(pTabX))*length(names(pTabY))) 
mutual_info_each=numeric(length(names(pTabX))*length(names(pTabY))) 
count=0
for(i in 1:length(pTabX)){
  
  for(j in 1:length(pTabY)){
    
    count=count+1
    p_xy[count]=sum( (as.numeric(names(pTabX[i]))==X) & (as.numeric(names(pTabY[j]))==Y)  )/length(X) #Doesn't matter whether length(X) / length(Y) are used
    
    if( p_xy[count]==0 ){
      mutual_info_each[count]=0 #This is to prevent NaN generation when p(x,y) = 0 ( log(0)=-Inf )
    }else{mutual_info_each[count]=p_xy[count]*log(p_xy[count]/(pTabX[i]*pTabY[j]),base=base)} 
  }
}
mutual_info=sum(mutual_info_each)
return(mutual_info)
}

#Weighted mutual information: this version doesn't really work (eg. compare X=Y=c(1,1,0,1,1) to X=Y=c(0,0,1,0,0))
#Ref: https://en.wikipedia.org/wiki/Mutual_information
mutualInfo_binary_weightedForNichols=function(X,Y,base=2){
  pTabX=table(X)/length(X)
  pTabY=table(Y)/length(Y)
  p_xy=numeric(length(names(pTabX))*length(names(pTabY))) 
  mutual_info_each=numeric(length(names(pTabX))*length(names(pTabY))) 
  count=0
  for(i in 1:length(pTabX)){
    
    for(j in 1:length(pTabY)){
      
      count=count+1
      p_xy[count]=sum( (as.numeric(names(pTabX[i]))==X) & (as.numeric(names(pTabY[j]))==Y)  )/length(X) #Doesn't matter whether length(X) / length(Y) are used
      
      if( p_xy[count]==0 ){
        mutual_info_each[count]=0 #This is to prevent NaN generation when p(x,y) = 0 ( log(0)=-Inf )
      }else{
        if( (as.numeric(names(pTabX[i]))==1 & as.numeric(names(pTabY[j]))==1)|
            (as.numeric(names(pTabX[i]))==-1 & as.numeric(names(pTabY[j]))==-1)
           ){ weight_xy=1
        } else if( (as.numeric(names(pTabX[i]))==1 & as.numeric(names(pTabY[j]))==-1)|
                   (as.numeric(names(pTabX[i]))==-1 & as.numeric(names(pTabY[j]))==1)
        ){weight_xy=0.5
        } else { weight_xy=0 }
        
        
        mutual_info_each[count]=weight_xy*p_xy[count]*log(p_xy[count]/(pTabX[i]*pTabY[j]),base=base)
        
        } 
    }
  }
  mutual_info=sum(mutual_info_each)
  return(mutual_info)
}


#To perform on a matrix like obj (pairwise MI), use this:
#mutinformation(t(Ternary_Data_324cutff_NAremoved[1:3,]) %>% as.data.frame) %>% natstobits

##Verify that it is consistent with my own written function
##mutualInfo_binary(Ternary_Data_324cutff_NAremoved[1,],Ternary_Data_324cutff_NAremoved[2,]) 
##mutualInfo_binary(Ternary_Data_324cutff_NAremoved[1,],Ternary_Data_324cutff_NAremoved[3,]) 
##mutualInfo_binary(Ternary_Data_324cutff_NAremoved[2,],Ternary_Data_324cutff_NAremoved[3,]) 




