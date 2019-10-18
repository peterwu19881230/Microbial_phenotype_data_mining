pos=read.table("Data/genomwide_metabolomic/zscore_pos.tsv")
neg=read.table("Data/genomwide_metabolomic/zscore_neg.tsv")


all=cbind(t(pos),t(neg)) #id from top to bottom matches those in the excel?


##compute pairwise MI matrix from Fuhrers (Failed to install on ada. Email sent to: ADA.coordinator@tamu.edu to ask for help on 10/8/2019)
if(!require(mpmi)){
  install.packages("mpmi")
  library(mpmi)
}

start.time = Sys.time()
mi_all=cminjk(t(all)) #running only the 1:10 features takes 7 sec, 1:20 takes 9 sec, 1:100 takes 2 min. Takes forever to run the fullset on my PC
end.time = Sys.time()
end.time - start.time 
fuhrer_mi=mi_all
save(fuhrer_mi,file="Data/fuhrer_mi.RData")