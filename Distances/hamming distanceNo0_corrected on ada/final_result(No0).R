#Read the binary_dist_hamming_x.RData and concatnate them

binary.coef.hamming_No0=list()
for(i in 1:39){
  load(paste("Data/binary_coef_hamming_No0_",i,".RData",sep=""))
  binary.coef.hamming_No0[[i]]=binary_coef_hamming_No0
}
library(data.table)
#Same as do.call("rbind", l) on data.frames, but much faster
##Ref: https://www.rdocumentation.org/packages/data.table/versions/1.10.4-2/topics/rbindlist
binary_coef_hamming_No0=rbindlist(binary.coef.hamming_No0)

#I summed up the same element, so I have to correct the distance by: total numbers(324) - similarity
binary_coef_hamming_No0$`Hamming Distance`=324-binary_coef_hamming_No0$`Hamming Distance`

save(binary_coef_hamming_No0,file="Data/binary_coef_hamming_No0.RData")
