if(!exists("Binary_Data")){load("PrimaryDataSets.RData")}
load("combs_separate.RData")

start.time=Sys.time()

part=30 # Change this to process different portions


binary_dist_hamming=cbind(combs_separate[[part]],NA) #create the distance obj as a molten matrix first
for(i in 1:dim(binary_dist_hamming)[1]){
  strain1=binary_dist_hamming[i,1]
  strain2=binary_dist_hamming[i,2]
  binary_dist_hamming[i,3]=sum(Binary_Data[strain1,]!=Binary_Data[strain2,],na.rm=T)
}



end.time=Sys.time()
end.time-start.time
colnames(binary_dist_hamming)=c("strain1","strain2","Hamming Distance")
binary_dist_hamming=as.data.frame(binary_dist_hamming)
save(binary_dist_hamming,file=paste("binary_dist_hamming_",part,".RData",sep=""))
