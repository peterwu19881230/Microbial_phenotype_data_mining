#Read the binary_dist_hamming_x.RData and concatnate them

binary.coef.hamming_No0_corrected=list()
for(i in 1:39){
  load(paste("Data/binary_coef_hamming_No0_corrected_",i,".RData",sep=""))
  binary.coef.hamming_No0_corrected[[i]]=binary_coef_hamming_No0_corrected
}
library(data.table)
#Same as do.call("rbind", l) on data.frames, but much faster
##Ref: https://www.rdocumentation.org/packages/data.table/versions/1.10.4-2/topics/rbindlist
binary_coef_hamming_No0_corrected=rbindlist(binary.coef.hamming_No0_corrected)
save(binary_coef_hamming_No0_corrected,file="Data/binary_coef_hamming_No0_corrected.RData")