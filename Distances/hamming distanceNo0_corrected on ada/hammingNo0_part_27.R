if(!exists("Binary_Data")){load("PrimaryDataSets.RData")}
load("combs_separate.RData")

start.time=Sys.time()

part=27 # Change this to process different portions


binary_coef_hamming_No0=cbind(combs_separate[[part]],NA) #create the distance obj as a molten matrix first
for(i in 1:dim(binary_coef_hamming_No0)[1]){
  strain1=binary_coef_hamming_No0[i,1]
  strain2=binary_coef_hamming_No0[i,2]
  binary_coef_hamming_No0[i,3]=sum(!(Binary_Data[strain1,]==0&Binary_Data[strain2,]==0)
    & Binary_Data[strain1,]==Binary_Data[strain2,] ,na.rm=T)
} ##If both are 0, don't count it. Count if 2 are the same.



end.time=Sys.time()
end.time-start.time
colnames(binary_coef_hamming_No0)=c("strain1","strain2","Hamming Distance")
binary_coef_hamming_No0=as.data.frame(binary_coef_hamming_No0)
save(binary_coef_hamming_No0,file=paste("binary_coef_hamming_No0_",part,".RData",sep=""))


