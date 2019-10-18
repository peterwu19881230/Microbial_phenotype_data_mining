#scp clustering/hamming\ distance\ on\ ada/*.* peterwu19881230@ada.tamu.edu:/home/peterwu19881230/

if(!exists("Binary_Data")){load("PrimaryDataSets.RData")}
strainID=1:3979
combs=t(combn(strainID,2)) #total No. of combination: 7914231 = 39*202929 (can be divided into 39 jobs)

combs_separate=list()
for(i in 1:39){
  combs_separate[[i]]=combs[(1+202929*(i-1)):(202929*i),]
}
save(combs_separate,file="clustering/hamming distance on ada/combs_separate.RData")










