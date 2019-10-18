if(!exists("Binary_Data")){load("PrimaryDataSets.RData")}
load("combs_separate.RData")

start.time=Sys.time()

part=16 # Change this to process different portions


  binary_coef_hamming_No0_corrected=cbind(combs_separate[[part]],NA) #create the distance obj as a molten matrix first
for(i in 1:dim(binary_coef_hamming_No0_corrected)[1]){
  strain1=binary_coef_hamming_No0_corrected[i,1]
  strain2=binary_coef_hamming_No0_corrected[i,2]
  binary_coef_hamming_No0_corrected[i,3]=sum(!(Binary_Data[strain1,]==0&Binary_Data[strain2,]==0)&Binary_Data[strain1,]==Binary_Data[strain2,] ,na.rm=T)-sum((Binary_Data[strain1,]==1&Binary_Data[strain1,]==-1)|(Binary_Data[strain1,]==-1&Binary_Data[strain1,]==1),na.rm=T) #Punish the score if they are (1,-1) or (-1,1)
} ##If both are 0, don't count it. Count if 2 are the same.



end.time=Sys.time()
end.time-start.time
colnames(binary_coef_hamming_No0_corrected)=c("strain1","strain2","Hamming Distance_no0_corrected")
binary_coef_hamming_No0_corrected=as.data.frame(binary_coef_hamming_No0_corrected)
save(binary_coef_hamming_No0_corrected,file=paste("binary_coef_hamming_No0_corrected_",part,".RData",sep=""))


