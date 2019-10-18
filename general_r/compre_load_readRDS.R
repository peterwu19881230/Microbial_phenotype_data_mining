# Compare load("x.RData") and readRDS("x.rds")


setwd("/Users/peterwu/Dropbox/Nichols-Data-mining")
load("Data/strain1.annot.allPCC.RData")

saveRDS(strain1.annot.allPCC,file="strain1.annot.allPCC.rds")

library(microbenchmark)
# Compare the two functions
microbenchmark(load("Data/strain1.annot.allPCC.RData"), 
                          readRDS("strain1.annot.allPCC.rds"), 
                          times = 10)

#Result:
#Unit: seconds
#                   expr                    min       lq     mean   median       uq      max    neval cld
#load("Data/strain1.annot.allPCC.RData") 4.014131 4.110251 4.537528 4.613306 4.889169 5.088521    10   a
#readRDS("strain1.annot.allPCC.rds")     4.147907 4.223073 4.720218 4.839004 5.045901 5.590443    10   a

#Conclusion: No significant difference. readRDS() is only slightly faster


