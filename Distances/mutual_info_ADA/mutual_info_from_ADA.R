#Get mutual information for the ternary data from ADA

##ssh peterwu19881230@ada.tamu.edu
##module load R/3.5.0-iomkl-2017b-recommended-mt

##getwd()
##save(Ternary_Data_324cutff,file="Ternary_Data_324cutff.RData")
##scp Ternary_Data_324cutff.RData peterwu19881230@ada.tamu.edu:/home/peterwu19881230/
##scp Distances/mutual_info_ADA/mutual_info_from_ADA.R peterwu19881230@ada.tamu.edu:/home/peterwu19881230/
##scp Distances/mutual_info_ADA/mutual_info_Ternary.lsf peterwu19881230@ada.tamu.edu:/home/peterwu19881230/
##bsub < mutual_info_Ternary.lsf
##bjobs #This can look up the running (or pending) jobs

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



load("Data/sourced/Ternary_Data_324cutff.RData")





strainID=1:3979 #change this to a small number for testing purpose
combs=t(combn(strainID,2)) #total No. of combination: 7914231 


start.time=Sys.time()


library(parallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl,c("mutualInfo_binary","Ternary_Data_324cutff_NAremoved","combs")) 
mutual_info_ternaryData=parApply(cl=cl,X=combs,MARGIN=1,FUN=function(comb){
  mutualInfo_binary(Ternary_Data_324cutff_NAremoved[comb[1],],Ternary_Data_324cutff_NAremoved[comb[2],])
})

mutual_info_ternaryData_table=cbind(combs,mutual_info_ternaryData)
save(mutual_info_ternaryData_table,file="mutual_info_ternaryData_table.RData")
##(from my mac):
##scp peterwu19881230@ada.tamu.edu:/home/peterwu19881230/mutual_info_ternaryData_table.RData /Users/peterwu/Dropbox/Nichols_Data_mining/Data/sourced/ 

end.time=Sys.time()
end.time-start.time



