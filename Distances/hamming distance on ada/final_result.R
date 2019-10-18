#Read the binary_dist_hamming_x.RData and concatnate them

binary.dist.hammnig=list()
for(i in 1:39){
  load(paste("Data/binary_dist_hamming_",i,".RData",sep=""))
  binary.dist.hammnig[[i]]=binary_dist_hamming
}
library(data.table)
#Same as do.call("rbind", l) on data.frames, but much faster
##Ref: https://www.rdocumentation.org/packages/data.table/versions/1.10.4-2/topics/rbindlist
binary_dist_hamming=rbindlist(binary.dist.hammnig)
save(binary_dist_hamming,file="Data/binary_dist_hamming.RData")